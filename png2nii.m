function png2nii
%hint: use XnConvert to convert obscure formats (e.g. SGI bw) to PNG
nSlice = 98;
slicedim = 256; %256ows,256 columns, 16-bit
img = zeros(slicedim, slicedim, nSlice);
img = uint8(img);
for i = 1: nSlice
   fnm = sprintf('head-%02d_result.png', i-1);
   slice = imread(fnm);
   %fclose(fileID);
   img(:,:,i) = slice;
end
img = flipdim(img,1);
img = permute(img,[2,1,3]);

inname = fullfile(spm('Dir'),'canonical','avg152T1.nii');
hdr = spm_vol([inname,',1']); 
hdr.fname = 'png.nii';
hdr.mat = [1 0 0 -(size(img,1)/2); 0 1 0 -(size(img,2)/2); 0 0 1 -(size(img,3)/2); 0 0 0 1]; 
hdr.dim(1) = size(img,1); hdr.dim(2) = size(img,2); hdr.dim(3) = size(img,3);
hdr.pinfo = [0;0;0];
spm_write_vol(hdr,img);
%end png2nii()
