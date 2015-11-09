function nii_mean_stdev_roi (fnm, roi)
%report descriptive statistics for portions of image named 'fnm' inside region of interest named 'roi'
% fnm : filename of 3D image with continuous brightness (.hdr/.img, .nii, or .nii.gz file)
% roi : filename of 3D masking image (.hdr/.img, .nii, .voi or .nii.gz file)
%Chris Rorden, 10/2015
% http://opensource.org/licenses/BSD-2-Clause
%SNR is mean/stdev http://www.ncbi.nlm.nih.gov/pubmed/17126038
%Example
% nii_mean_stdev_roi;
% nii_mean_stdev_roi('a.nii','b.voi');

if ~exist('fnm','var') %no files
 fnm = spm_select(1,'^.*\.(gz|voi|img|nii)$','Select image');
end;
if ~exist('roi','var') %no region of interest
 roi = spm_select(1,'^.*\.(gz|voi|img|nii)$','Select mask');
end;
[hdr, img] = read_volsSub (fnm); %#ok<ASGLU>
[rhdr, rimg] = read_volsSub (roi); %#ok<ASGLU>
if min(size(img) == size(rimg)) ~= 1
    error('Image and mask must have the same dimensions');
end
mimg = img(rimg > 0); %find voxels of image that are in mask
fprintf('Descriptives for image %s\n', fnm);
fprintf(' masked with image %s\n', roi);
fprintf(' number of voxels in mask : %d\n', numel(mimg(:)) );
fprintf(' intensity range %g..%g\n', min(mimg(:)), max(mimg(:)) );
fprintf(' mean intensity %g\n', mean(mimg(:)) );
fprintf(' standard deviation %g\n', std(mimg(:)) );
fprintf(' median %g\n', median(mimg(:)) );
%end nii_mean_stdev_roi() - local sub-functions follow

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
