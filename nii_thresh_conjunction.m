function nii_thresh_conjunction(fnms)
%Find voxels that survive threshold in ALL tests
% fnms : file name[s] of thresholded image[s] (optional)
%Notes
% The 'nii_thresh_conjunction' test is the valid conjunction test discussed in:
%    Nichols T, Brett M, Andersson J, Wager T, Poline JB. Valid conjunction
%    inference with the minimum statistic. Neuroimage. 2005 Apr 15;25(3):653-60.
%    This test gives the p-value of the z-values under the conjunction null,
%    i.e. the union of the null hypotheses for all terms
% http://nipy.sourceforge.net/nipy/devel/api/generated/nipy.modalities.fmri.glm.html#nipy.modalities.fmri.glm.Contrast
%Examples
% nii_thresh_conjunction(strvcat('a.nii','b.nii'));
if ~exist('fnms','var')
	fnms = spm_select(inf,'^thresh.*\.nii|.gz$','Select threshold images to conjunction');
end
if size(fnms,1) < 2
   error('You need to specify at least two images two conjoin');
end

fnm = deblank(fnms(1,:));
[hdr, minimg] = read_volsSub (fnm);
minimg(isnan(minimg)) = 0;
for i=2:size(fnms,1)
    fnm = deblank(fnms(i,:));
    [hdr, img] = read_volsSub (fnm);
    img(isnan(img)) = 0;
    if min(size(minimg) ~= size(img))
       error('image dimensions do not match : each image must be resliced to the same size');
    end
    minimg = min(minimg,img);
end;
%see if any voxels survive
minimg(minimg < 0) = 0;
minimg(minimg > 0) = 1;
numPos = sum((minimg(:)) > 0);
if  (numPos == 0)
   error('No positive voxels survive');
end
%save data
[pth nm ext] = spm_fileparts(hdr.fname);
hdr.fname = fullfile(pth, ['conjoin_' nm ext]);
spm_write_vol(hdr,minimg);


%local functions follow

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

