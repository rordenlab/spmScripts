function nii_mask(imgNames, maskName, thresh, val)
%Masks images imgName with image maskName: values less that thresh in mask are set to val in imgNames
%if a voxel in the mask is less than thresh corresponding voxel in input is set to val
% imgNames   : Filenames of image(s) to mask
% maskName   : Mask filename
% thresh     : Change voxels darker than this in the template image 
% val        : Value to set dark voxels to
%Output:
% masked images with 'm' prefix
%Example
%  nii_mask('T1.nii','lesion.nii',0.5, 0.0); %voxels darker than 0.5 in mask set to 0.0 in images

if nargin<1, imgNames = spm_select(Inf,'image','Select images for masking'); end
if nargin<2, maskName = spm_select(1,'image','Select mask'); end
if nargin<3, thresh = 0.025; end
%if nargin<4, val = 0; end
hdrM = spm_vol(maskName);
if numel(hdrM) > 1 
    error('Error: mask has multiple volumes, please explicitly specify masking volume, e.g. ''~/tpm.nii,1'' ');
end;
imgM = spm_read_vols(hdrM);
for j=1:size(imgNames,1)  
  hdr = spm_vol(deblank(imgNames(j,:)));
  if hdr.dim ~= hdrM.dim
    disp('mask and source must have same dimensions!');
    break; 
  end;
  img = spm_read_vols(hdr);
  clipped = sum(imgM(:) < thresh);
  if nargin<4, val = min(img(:)); end
  imgM(img > val) = 0;
  %img(imgM < thresh) = val;
  img(imgM > thresh) = val;
  
  [pth,nm,xt, ~] = spm_fileparts(hdr.fname);
  hdr.fname = fullfile(pth, ['m' nm xt]);  
  img(isnan(img)) = 0; % use ~isfinite instead of isnan to replace +/-inf with zero
  spm_write_vol(hdr,img);
  fprintf('Image %s had %d voxels set to %f because they were < %f in %s\n',hdr.fname, clipped, val, thresh, hdrM.fname); 
end
%end nii_mask()
