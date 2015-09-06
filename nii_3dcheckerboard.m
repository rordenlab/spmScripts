function nii_3dcheckerboard
%Create 3D object to test interpolation

dim = [15 15 1];
img = zeros(dim(1),dim(2),dim(3));
img(1:2:end) = 200;
img(1) = 50;

inname = fullfile(spm('Dir'),'canonical','avg152T1.nii');
hdr = spm_vol([inname,',1']); 
hdr.fname = ['check' num2str(dim(1)) 'x' num2str(dim(2)) 'x' num2str(dim(3)) '.nii'];
hdr.mat = [1 0 0 -(size(img,1)/2); 0 1 0 -(size(img,2)/2); 0 0 1 -(size(img,3)/2); 0 0 0 1]; 
hdr.dim(1) = size(img,1); hdr.dim(2) = size(img,2); hdr.dim(3) = size(img,3);
hdr.pinfo = [0;0;0];
spm_write_vol(hdr,img);