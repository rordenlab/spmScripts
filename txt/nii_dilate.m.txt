function nii_dilate(P, dilatemm, thresh);
%dilate binary volume
% FWHM: Full-width half maximum for Gaussian smooth in mm, thresh is threshold for boundary
% Thresh: Threshold for smoothed image
%   if thresh = 0.5: output will be approximately same size as input
%   if thresh > 0.5 then erosion, if thresh < 0.5 then dilation
%Example:
%  nii_dilate('T1_brain_mask.nii',8, 0.25); %dilate

if nargin <1 %no files
 P = spm_select(inf,'image','Select images to smooth');
end;
if nargin < 2 %no FWHM specified
 dilatemm = 8;
end;
if nargin < 3 %no FWHM specified
 thresh = 0.5;
end;
thresh = thresh *255;

for i=1:size(P,1)
  ref = deblank(P(i,:));
  [pth,nam,ext] = spm_fileparts(ref);
  src = fullfile(pth,[nam ext]);
  %smooth data
  %spm_smooth(src,smth,dilatemm); 
  %threshold smoothed data...
  hdr = spm_vol(src);
  in = spm_read_vols(hdr);
  in = in*255;
  smth = in;
  spm_smooth(in,smth,dilatemm);
  hdr.fname = fullfile(pth,['x',  nam, ext]);
  mask = smth;
  smth= zeros(size(mask));
  smth((mask >= thresh)) = 1; 

  hdr.dt(1)=2;%save binary data as bytes: uint8=2; int16=4; int32=8; float32=16; float64=64
  spm_write_vol(hdr,smth);
  %delete (hdr.fname); %this is 32-bit, but we will next save as 8-bit: SPM overwrites rather than erases images
  %hdr.dt(1)=2;%save binary data as bytes: uint8=2; int16=4; int32=8; float32=16; float64=64
  %spm_write_vol(hdr,w);
end

%%the the code below works, but grows diamonds not spheres.... 
%   also, dilation in voxels not mm
% if nargin <1 %no gray
%  wi = spm_select(1,'image','Select volume to dilate');
% end;
% if nargin <2 %no gray
%     ndilate = 10; %1=basic, 2=thorough
% end;
% 
% if ischar(wi), wi = spm_vol(wi); end;
% 
% mask = spm_read_vols(wi)*255;
% 
% kx=[0.25 0.5 0.25];
% ky=[0.25 0.5 0.25];
% kz=[0.25 0.5 0.25];
% mn = min(mask(:));
% for j=1:ndilate,
%     spm_conv_vol(mask,mask,kx,ky,kz,-[1 1 1]);
% end;
% w= zeros(size(mask));
% w((mask > mn)) = 1;
% 
% [pth,nam,ext]=fileparts(wi.fname);
% wi.fname = fullfile(pth,['d',  nam, ext]);
% spm_write_vol(wi,w);
