function imgD = dilate_erode(fname, doDilate)
%dilate or erode binary volume - either return image or save to disk
% fname : file name of a NIfTI image
% doDilate : if false then erode, else dilate
%Example:
%  dilate_erode; %use GUI
%  dilate_erode('ch2.nii'); %save dilated image to disk
%  dilate_erode('ch2.nii', false); %save eroded image to disk
%  im = dilate_erode('ch2.nii'); %save dilated image to disk
%  im = dilate_erode(im); %dilate memory a 2nd time
%Chris Rorden 2015 - BSD 2-clause license

if nargin < 1 %no files
 fname = spm_select(inf,'image','Select image to dilate');
end;
if ~isa(fname, 'char')
    img = fname;
else
    hdr = spm_vol(fname);
    img = spm_read_vols(hdr);
end
dims = size(img);
if size(img,3) < 3, error('Not a 3D volume'); end;
img = img(:); %convert 3D to 1D vector
thresh = (max(img)-min(img))/2 +  min(img);
img = (img > thresh); %binarize 0,1
if (sum(img == 1) < 1) || (sum(img == 0) < 1), error('No variability in image'); end;
imgD = img; %dilated image
imgD = imgD + [img(2:end); 0]; %dilate left
imgD = imgD + [0; img(1:end-1)]; %dilate right
imgD = imgD + [img(dims(1)+1:end); zeros(dims(1),1)]; %dilate anterior
imgD = imgD + [zeros(dims(1),1); img(1:end-dims(1))]; %dilate posterior
sliceVox = dims(1) * dims(2); %voxels per slice
imgD = imgD + [img(sliceVox+1:end); zeros(sliceVox,1)]; %dilate inferior
imgD = imgD + [zeros(sliceVox,1); img(1:end-sliceVox)]; %dilate superior
if (nargin > 1) && (~doDilate) %no files
    imgD = (imgD > 6); %binarize ERODED
else
    imgD = (imgD > 0); %binarize DILATED
end
imgD = reshape(imgD, dims); %convert from 1D to 3D
if (nargout > 0) return; end;%return image data
%Next lines only executed if no output specified
% save data to disk 
if ~exist('hdr','var'), error('Unable to save raw data'); end;
[pth,nam,ext] = spm_fileparts(fname);
hdr.fname = fullfile(pth,['x',  nam, ext]);
spm_write_vol(hdr,imgD);
%end dilate_erode()

  