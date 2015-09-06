function nii_i3m2nii
%Convert i3m file to NIfTI
% http://fossies.org/dox/visit2.9.0/I3MConverter_8cpp_source.html
fnm = 'beetle.i3m';
fileID = fopen(fnm);
iMagic = fread(fileID,1,'uint32');
iVersion = fread(fileID,1,'uint32');
vSizeX  = fread(fileID,1,'uint32');
vSizeY  = fread(fileID,1,'uint32');
vSizeZ  = fread(fileID,1,'uint32');
if (iMagic ~= 69426942) || (iVersion ~= 1)
    error('Magic or version do not match');
    fclose(fileID);
end
img = fread(fileID,vSizeX*vSizeY*vSizeZ,'uint32');
img = bitshift(img, -24);
img = reshape(img, vSizeX,vSizeY,vSizeZ);
%optional: swizzle dimensions
img = permute(img,[1,3,2]);
img = flipdim(img,2);


inname = fullfile(spm('Dir'),'canonical','avg152T1.nii');
hdr = spm_vol([inname,',1']); 
hdr.fname = 'i3m.nii';
hdr.mat = [1 0 0 -(size(img,1)/2); 0 1 0 -(size(img,2)/2); 0 0 1 -(size(img,3)/2); 0 0 0 1]; 
hdr.dim(1) = size(img,1); hdr.dim(2) = size(img,2); hdr.dim(3) = size(img,3);
hdr.pinfo = [0;0;0];
spm_write_vol(hdr,img);
