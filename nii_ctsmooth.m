function nii_ctsmooth
fnm = 'CT_Philips.nii'
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
oimg = img + 0.0;
mn = -20;
mx = 80;
img(img < mn) = mx;
img(img >= mx) = NaN;
simg = img + 0.0;
spm_smooth(img,simg, 3);
size(simg)
size(oimg)
size(img)
oimg(isfinite(img)) = simg(isfinite(img));

[pth nm ext] = spm_fileparts(fnm);
hdr.fname = fullfile(pth, ['s' nm ext]);  
%hdr.pinfo = [1; 0; 352];
spm_write_vol(hdr,oimg);
return
imgBin = 1.0 - (img >= mx);

[pth nm ext] = spm_fileparts(fnm);
hdr.fname = fullfile(pth, ['z' nm ext]);  
hdr.pinfo = [1; 0; 352];
spm_write_vol(hdr,imgBin);