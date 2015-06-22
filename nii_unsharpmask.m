function nii_unsharpmask (fnm, smoothVox)
% Creates sharpened image with prefix 's'
%  fn : filename of image to sharpen
%  smoothVox : ammount of smoothing
% Example
%   nii_unsharpmask ('chris_t1.nii');

if nargin < 1 %no files
 fnm = spm_select(1,'image','Select image to sharpen');
end;
if nargin < 2 %smoothing not specified
    smoothVox = 2; 
end
%load image
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
%blur image
imgBlur = img+0; %+0 forces new matrix
spm_smooth(img,imgBlur,smoothVox,0); %create blurred image
dx = imgBlur - img;
img = img-dx;
img(img < 0) = 0;
%save image
[pth,nm,ext] = spm_fileparts(fnm);
hdr.fname = fullfile(pth, ['z' nm ext]);  
fprintf('Note sharpened images magnify noise. Saving as %s\n', hdr.fname);
spm_write_vol(hdr,img);