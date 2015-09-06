function nii_tight_crop(fnm)
%discard exterior rows/columns/slices that are zeros, output has 'z' prefix
% fnms : file name of image (optional)
%Examples
% nii_tight_crop; %use GUI
% nii_tight_crop('img.nii');

if ~exist('fnm','var') %no filename specified
	fnm = spm_select(1,'image','Select image[s] for NaN removal'); 
end
%load image
hdr = spm_vol(deblank(fnm));
img = spm_read_vols(hdr);
if size(img,4) > 1
    error('%s designed for 3D images with only a single volume\n',mfilename);
end
img(isnan(img)) = 0;%zero nans
%collapse in x dim
for lo = 1: size(img,1)
    if max(max(img(lo,:,:))) > 0
        break
    end
end
if lo > 1
   img(1:(lo-1),:,:) = [] ;
end
for hi = size(img,1): -1 : 1
    if max(max(img(hi,:,:))) > 0
        break
    end
end
if hi < size(img,1)
   img((hi+1):end, :, :) = [] ;
end
%collapse in y dim
for lo = 1: size(img,2)
    if max(max(img(:, lo,:))) > 0
        break
    end
end
if lo > 1
   img(:,1:(lo-1),:) = [] ;
end
for hi = size(img,2): -1 : 1
    if max(max(img(:,hi,:))) > 0
        break
    end
end
if hi < size(img,2)
   img(:,(hi+1):end, :) = [] ;
end
%collapse z dim
for lo = 1: size(img,3)
    if max(max(img(:,:,lo))) > 0
        break
    end
end
if lo > 1
   img(:,:,1:(lo-1)) = [] ;
end
for hi = size(img,3): -1 : 1
    if max(max(img(:,:,hi))) > 0
        break
    end
end
if hi < size(img,3)
   img(:,:,(hi+1):end) = [] ;
end
%abort if there are no voxels to crop
if (hdr.dim(1) == size(img,1)) && (hdr.dim(2) == size(img,2)) && (hdr.dim(3) == size(img,3))
    fprintf('Unable to crop this image: positive intensity observed to all edges\n');
    return
end
%save cropped image
hdr.dim(1) = size(img,1);
hdr.dim(2) = size(img,2);
hdr.dim(3) = size(img,3);
[pth nm ext] = spm_fileparts(fnm);
hdr.fname = fullfile(pth, ['z' nm ext]); 
spm_write_vol(hdr,img);
%end nii_tight_crop()
