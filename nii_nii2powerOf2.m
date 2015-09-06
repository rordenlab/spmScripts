function nii_nii2powerOf2 (fnm)
%Pad NIfTI image so dimensions are now a power of two
%Some older graphics cards require power-of-two for textures
% fnm : image to modify
%Examples
% nii_nii2powerOf2; %use GUI
% nii_nii2powerOf2('vx256a.nii')

if ~exist('fnm','var')
	fnm = spm_select(1,'image','Select image to resize'); 
end
%load input image
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
dim = size(img); %input dimensions
%pad output image
odim = 2.^nextpow2(dim); %output dimensions
if all(odim == dim)
   error('Nothing to do: input image dimensions are power of 2');
end
dx = floor( (odim - dim)/2); %half of difference
oimg = zeros(odim);
oimg(dx(1)+1:dx(1)+dim(1), dx(2)+1:dx(2)+dim(2), dx(3)+1:dx(3)+dim(3)) = img;
%create output header
ohdr = hdr;
%since we have added slices we need to shift the origin...
v2m = hdr.mat; %voxel2mm transform
ohdr.mat(1:3,4) = v2m(1:3,4)' - dx*v2m(1:3,1:3)';
ohdr.dim = odim;
[pth nm ext] = spm_fileparts(fnm);
ohdr.fname = fullfile(pth, ['z' nm ext]);  
spm_write_vol(ohdr,oimg);
%nii_nii2powerOf2()        
