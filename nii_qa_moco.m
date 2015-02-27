function nii_qa_moco (imgName, rpName)
%Quality assurance: identify odd images based on 4D fMRI series and motion correction parameters 
% imgName : filename(s) of images to average, either single 4D volume or a series of 3D volumes
% rpName : name of motion corection file created by SPM (e.g. rp_img.txt)
%Chris Rorden, 2014
% BSG License http://opensource.org/licenses/BSD-2-Clause
%SNR is mean/stdev http://www.ncbi.nlm.nih.gov/pubmed/17126038
%DVAR from Powers et al. http://www.ncbi.nlm.nih.gov/pubmed/22019881 
%Example
% nii_qa_moco; %use GUI
% nii_qa_moco('4d.nii','rp_4d.txt'); %4d data
% nii_qa_moco('REST_LM1001.nii', '') %no moco file
% nii_qa_moco(strvcat('3d1.nii','3d2.nii','3d3.nii'),'rp_3d.txt');

% 0: check inputs
if ~exist('spm','file')
    error('%s requires SPM', mfilename);
end
if ~exist('imgName','var') %no image specified
 imgName = spm_select(inf,'^.*\.(gz|voi|img|nii)$','Select images to inspect');
end;
if ~exist('rpName','var') %no realignment parameter file specified
    %rpName = spm_select(inf,'^rp.*\.txt$','Select realignment parameter file');
    rpName = spm_select(inf,'^rp.*\.txt|^.*\.par$','Select realignment parameter file');
end
if (size(imgName,1) < 1) ||( numel(imgName) < 1), return; end;

% 1: Load image data
[hdr, img4d] = loadImg (imgName);
if size(img4d,4) < 3, error('Image data must have many volumes'); end;
% 2: Computemean, standard deviation and SNR measures
meanImg = mean(img4d,4);
sdImg = std(img4d,0,4);
% 3: Save images of mean, stanard deviation and SNR
saveImgSub(['mean_of_' num2str(size(img4d,4)) '.nii'],hdr,meanImg);
saveImgSub(['stdev_of_' num2str(size(img4d,4)) '.nii'],hdr,sdImg);
saveImgSub(['snr_of_' num2str(size(img4d,4)) '.nii'],hdr,meanImg ./ sdImg);
% 4: compute z-scores
zscore1d = zScoreSub(meanImg, sdImg, img4d);
if isempty(rpName)
    %report variability (without motion parameters)
    report1Text(zscore1d);
    report1Plot(zscore1d);
else
    % 5: estimate temporal derivative of timecourse variability
    dvars1d = deltaIntensitySub (meanImg, img4d, hdr);
    % 6: estimate head motion in mm
    fd1d = framewiseDisplacementSub (rpName);
    % 7: report variability
    reportText(zscore1d, fd1d, dvars1d);
    reportPlot(zscore1d, fd1d, dvars1d);
end
%end nii_qa_moco

function report1Plot (zscore1d)
zscore1d = zscore1d - min(zscore1d);
x = 1:numel(zscore1d);
figure;
plot(x,zscore1d/max(zscore1d),'r')
legend('Intensity Z');
xlabel('Time') 
ylabel('Amplitude (Normalized)') 
%end report1Plot()

function report1Text(zscore1d)
%report statistics for a single 4D image
diary qa.tab
fprintf('Volume\tZ-Score\n');
for i = 1 : numel(zscore1d)
    fprintf('%d\t%g\n',i, zscore1d(i));
end
[mx, mxi] = max(zscore1d);
fprintf('Most extreme Z-score %g (volume number %d)\n',mx, mxi);
diary off
%end report1Text()

function reportPlot (zscore1d, fd1d, dvars1d)
zscore1d = zscore1d - min(zscore1d);
fd1d = fd1d - min(fd1d);
dvars1d = dvars1d - min(dvars1d);
x = 1:numel(zscore1d);
figure;
plot(x,zscore1d/max(zscore1d),'r', x,fd1d/max(fd1d),'b', x, dvars1d/max(dvars1d),'g')
legend('Intensity Z','Position Change','Intensity Change');
xlabel('Time') 
ylabel('Amplitude (Normalized)') 
%end reportPlot()

function reportText(zscore1d, fd1d, dvars1d)
%report statistics for a single 4D image
diary qa_moco.tab
fprintf('Volume\tZ-Score\tFramewiseDisplacement\tFramewiseVariability\n');
for i = 1 : numel(zscore1d)
    fprintf('%d\t%g\t%g\t%g\n',i, zscore1d(i), fd1d(i), dvars1d(i));
end
M = corrcoef([zscore1d(:) fd1d(:)]);
[mx, mxi] = max(zscore1d);
fprintf('Most extreme Z-score %g (volume number %d): Correl PositionChange->IntensityZ R=%g\n',mx, mxi, M(2));
[mx, mxi] = max(fd1d);
fprintf('Most extreme Position Change (Framewise Displacement) %g mm (volume number %d)\n',mx, mxi);
[mx, mxi] = max(dvars1d);
M = corrcoef([fd1d(:) dvars1d(:)]);
fprintf('Most extreme Intensity Change (Framewise Variability) %g (volume number %d): Correl PositionChange->IntensityChange R=%g\n',mx, mxi,M(2));
diary off
%end reportText()

function zscore1d = zScoreSub (meanImg, sdImg, img4d)
%compute how unusual each volume is
zscore1d = zeros(1, size(img4d,4)); %pre-allocate array
for i = 1 : size(img4d,4)
    img = img4d(:,:,:,i);
    img = (img-meanImg)./sdImg; %transform to z-score
    zscore1d(i) = mean(abs(img(isfinite(img(:))))); %remove not-a-number values, compute mean
end

function dvars1d = deltaIntensitySub (meanImg, img4d, hdr)
%see DVARS, Powers et al. http://www.ncbi.nlm.nih.gov/pubmed/22019881 
% 0 : create a brain mask, here we use 1/8 mean 
ave = max(meanImg(:)); %mean(meanImg(:));
mn = min(meanImg(:));
thresh = mn + (ave-mn)/8;
mask = (meanImg >= thresh);%region of interest is all voxels that except the very darkest
dvars1d = zeros(1, size(img4d,4)); %pre-allocate array
for i = 2 : size(img4d,4)
    img = img4d(:,:,:,i) - img4d(:,:,:,i-1);
    img = img(mask);
    dvars1d(i) = var(img(:));
end
saveImgSub( 'mask.nii',hdr,mask);
%end deltaIntensitySub()

function fd = framewiseDisplacementSub (rpName)
%see Powers et al. http://www.ncbi.nlm.nih.gov/pubmed/22019881 
rp = load(rpName); %spm_load(rpName); %select the rp*.txt file
rpdx = rp;%preallocate, set 
for i = 1:6 
    rpdx(2:end,i) = rp(2:end,i) - rp(1:end-1,i);
end
[~, ~, ext] = fileparts(rpName);
if strcmpi(ext,'.par')
    % http://fsl.fmrib.ox.ac.uk/fsl/fsl-4.1.9/possum/input.html
    % https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;cda6e2ea.1112
    fprintf('Assuming motion parameters saved in format "Xrad Yrad Zrad Xmm Ymm Zmm" \n');
    for i = 1:3 %convert rotations from radians to mm
        rpdx(1:end,i) = rpdx(1:end,i) * 50; %brain assumed 50mm radius: circle circumference=2*pi*r http://en.wikipedia.org/wiki/Great-circle_distance
    end 
else
    fprintf('Assuming motion parameters saved in format "Xmm Ymm Zmm Xrad Yrad Zrad" \n');
    for i = 4:6 %convert rotations from radians to mm
        rpdx(1:end,i) = rpdx(1:end,i) * 50; %brain assumed 50mm radius: circle circumference=2*pi*r http://en.wikipedia.org/wiki/Great-circle_distance
    end
end
rpdx = abs(rpdx);
fd = sum(rpdx,2);
%end framewiseDisplacementSub()

function saveImgSub(fname, hdr, img)
%save data as a NIfTI format image
hdr.fname = fname;
hdr.dt    =[16,0]; %set data type uint8=2; int16=4; int32=8; float32=16; float64=64
hdr.pinfo = [1;0;0]; %reset scale slope and intercept
spm_write_vol(hdr,img);
%end saveImgSub()

function [hdr, img] = loadImg (fnm)
unm = unGzSub (fnm); %convert FSL .nii.gz to .nii
hdr = spm_vol(unm); %load header data
for i = 1: numel(hdr) %ignore 4D realignment .mat file that would provide different rotation to each image
    hdr(i).mat = hdr(1).mat;
end
img = spm_read_vols(hdr); %load image data
hdr = hdr(1);
%if input was nii.gz, remove .nii (otherwise FSL complains)
if strcmpi(fnm,unm), return; end;
for i = 1: size(unm,1)
   delete(unm(i,:));
end
%end deleteUnGzSub

function unfnm = unGzSub (fnm)
unfnm = [];
for i = 1: size(fnm,1)
    unfnm = [unfnm; uGzSub(fnm(i,:))]; %#ok<AGROW>
end
%end unGzSub

function fnm = uGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));  
elseif strcmpi(ext,'.voi') %.voi -> 
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
end;  
%end uGzSub()