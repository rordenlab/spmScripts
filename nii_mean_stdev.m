function nii_mean_stdev (fnms, normBrightness, outname)
%Given multiple volumes, generate mean, standard deviation, SNR maps and report unusual images. Useful for quality assurance
% fnms : filenames to average (optional)
% normBrightness : if false (0) raw intensity is used. if true (1) image intensity
%                  scaled so minimum is 0 and median of non-zero voxels is 1
% outname : optional prefix appended to output image name
%Chris Rorden, 2014
% http://opensource.org/licenses/BSD-2-Clause
%SNR is mean/stdev http://www.ncbi.nlm.nih.gov/pubmed/17126038
%Example
% nii_mean_stdev;
% nii_mean_stdev(strvcat('a.nii','b.nii'),0);

if ~exist('fnms','var') %no files
 %fnms = spm_select(inf,'image','Select images to average');
 fnms = spm_select(inf,'^.*\.(gz|voi|img|nii)$','Select images to average');
end;
if ~exist('normBrightness','var') %no files
    normBrightness = str2double(cell2mat(inputdlg('Scale image intensity to normalize brightness? (0=no,1=yes):', 'Adjust brightness', 1,{'0'})));
end;
if ~exist('outname', 'var')
    outname = '';
end
if size(fnms,1) < 1
    return;
elseif size(fnms,1) == 1 %4D: compute mean and stDev using inbuilt Matlab functions
    [hdr, meanImg, sdImg, n] = statSub4D(normBrightness, fnms);
else %3D: compute mean and stDev using Welford one-pass http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    [hdr, meanImg, sdImg, n] = statSub3D(normBrightness, fnms);
end;
if (n < 2)
    error('Not enough images to estimate mean and stdev');
end
%save statistical images (optional)
saveImgSub([outname 'mean_of_' num2str(n) '.nii'],hdr,meanImg);
saveImgSub([outname 'stdev_of_' num2str(n) '.nii'],hdr,sdImg);
saveImgSub([outname 'snr_of_' num2str(n) '.nii'],hdr,meanImg ./ sdImg);
%report mean z-score for each image
if size(fnms,1) == 1
    reportZSub4D(meanImg, sdImg, normBrightness, fnms)
else
    reportZSub3D(meanImg, sdImg, normBrightness, fnms)
end
%end nii_mean_stdev() - local sub-functions follow

function reportZSub4D(meanImg, sdImg, normBrightness, fnms)
%report statistics for a single 4D image
[pth, nam, ext, vol] = spm_fileparts(fnms); %#ok<NASGU>
fnm = fullfile(pth, [nam ext]); %remove volume index 'vol' 
[hdr4d, img4d] = read_volsSub (fnm); %#ok<ASGLU>
mx = 0; %max value
mxi = 1; %index of max Z-score
sm = 0; %sum
for i = 1 : size(img4d,4)
    img = img4d(:,:,:,i);
    if normBrightness, img = normBrightnessSub(img); end; 
    img = (img-meanImg)./sdImg; %transform to z-score
    mn = mean(abs(img(isfinite(img(:))))); %remove not-a-number values, compute mean
    %saveImgSub(['z_' num2str(i) '.nii'],hdr4d(1),img);
    fprintf('Volume %d has a mean absolute Z-score of %g\n',i, mn);
    if (mn > mx)
        mx = mn;
        mxi = i;
    end
    sm = sm + mn;
end
fprintf('Most extreme volume was %d with a mean absolute Z-score of %g, mean Z-score for all volumes is %g\n\n',mxi, mx, sm/size(img4d,4));
%end reportZSub4D()

function reportZSub3D(meanImg, sdImg, normBrightness, fnms)
%compute mean and stDev using Welford one-pass http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
if (numel(fnms) < 1), return; end;
mx = 0; %max
mxi = 1; %index of maximum
sm = 0; %sum
for i=1:size(fnms,1)
	fnm = deblank(fnms(i,:));
    [hdr, img] = read_volsSub (fnm); %#ok<ASGLU>
    if normBrightness, img = normBrightnessSub(img); end; 
    img = (img-meanImg)./sdImg; %transform to z-score
    mn = mean(abs(img(isfinite(img(:))))); %remove not-a-number values, compute mean
    fprintf('%s has a mean absolute Z-score of %g\n',fnm, mn);
    if (mn > mx)
        mx = mn;
        mxi = i;
    end 
    sm = sm + mn;
end
fprintf('Most extreme image was %s with a mean absolute Z-score of %g, mean Z-score for all volumes is %g\n',deblank(fnms(mxi,:)), mx, sm/size(fnms,1));
%end reportZSub3D()

function [hdr, meanImg, sdImg, n] = statSub4D(normBrightness, fnms)
%compute statistics for a single 4D image
[pth, nam, ext, vol] = spm_fileparts(fnms); %#ok<NASGU>
fnm = fullfile(pth, [nam ext]); %remove volume index 'vol' 
if (exist(fnm,'file') == 0); error('%s unable to find image %s\n',which(mfilename),fnm); end;
[hdr4d, img4d] = read_volsSub (fnm);
hdr = hdr4d(1);
n = size(img4d,4);
if normBrightness 
    for v = 1 : n 
        img4d(:,:,:,v) = normBrightnessSub(img4d(:,:,:,v)); 
    end
end;
meanImg = mean(img4d,4);
sdImg = std(img4d,0,4);
%end statSub4D()

function [hdr, meanImg, sdImg, n] = statSub3D(normBrightness, fnms)
%compute statistics for a series of 3D images
hdr = [];
meanImg = [];
sdImg = [];
n = 0;
[hdr, meanImg, sdImg, n] = sumSub(hdr, meanImg, sdImg, n, normBrightness, fnms);
if n < 2, return; end;
sdImg = sqrt(sdImg/(n-1)); %convert to standard deviation
%end statSub3D

function [hdr, meanImg, m2Img, n] = sumSub(hdr, meanImg, m2Img, n, normBrightness, fnms)
%compute mean and stDev using Welford one-pass http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
if (numel(fnms) < 1), return; end; 
for i=1:size(fnms,1)
	fnm = deblank(fnms(i,:));
    [hdr, img] = read_volsSub (fnm);
    if normBrightness, img = normBrightnessSub(img); end;
    if isempty(meanImg), meanImg = zeros(size(img)); end;
    if isempty(m2Img), m2Img = zeros(size(img));  end;
    n = n + 1;
    delta = img - meanImg;
    meanImg = meanImg + (delta / n);
    m2Img = m2Img + delta.*(img-meanImg);
    %fprintf('%f\t%f\t%f\n',img(1),meanImg(1),m2Img(1));
end
%end sumSub()

function [hdr, img] = read_volsSub (fnm)
[fnm, isGz] = unGzSub (fnm); %convert FSL .nii.gz to .nii
hdr = spm_vol(fnm); %load header data
img = spm_read_vols(hdr); %load image data
if (isGz), delete(fnm); end; %remove .nii if we have .nii.gz
%end read_volsSub()

function [fnm, isGz] = unGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm);
isGz = false;
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));  
    isGz = true;
elseif strcmpi(ext,'.voi') %.voi -> 
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
    isGz = true;
end;  
%end unGzSub()

function saveImgSub(fname, hdr, img)
%save data as a NIfTI format image
hdr.fname = fname;
hdr.dt    =[16,0]; %set data type uint8=2; int16=4; int32=8; float32=16; float64=64
hdr.pinfo = [1;0;0]; %reset scale slope and intercept
spm_write_vol(hdr,img);
%end saveImgSub()

function img = normBrightnessSub(img)
img = img - min(img(:));
mdn = median(img(img > 0));
if mdn == 0, return; end;
img = img / mdn;
%end normBrightnessSub()