function [cropName4D] = nii_cropVolumes (imgName4D, skipVol, nVol);
%given 4D NIfTI image creates NIfTI image with only selected input volumes
% imgName4D : name of source image (4D)
% skipVol   : number of initial volumes deleted from output
% nVol      : number of input volumes
%Examples
% nii_cropVolumes('img.nii',1,inf); %drop first volume
% nii_cropVolumes('img.nii',0,2); %preserve only first two volumes
% nii_cropVolumes('img.nii',1,2); %output has volumes 2,3 from input
%See also
% nii_numVolumes reports number of volumes in an image
if ~exist('imgName4D')
 imgName4D = spm_select(1,'image','Select 4D image to crop');
end
[pth,nam,ext,vol] = spm_fileparts( deblank (imgName4D));
imgName4D = fullfile(pth,[nam, ext]); %remove volume is 'img.nii,1' -> 'img.nii'
if ~exist(imgName4D), return; end;
if ~exist('skipVol','var') || ~exist('nVol','var')
    prompt = {'Skip first N volumes:','Retain N volumes'};
    dlg_title = 'Values for cropping';
    num_lines = 1;
    def = {'0','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    skipVol = str2num(answer{1});
    nVol = str2num(answer{2});
end;
if (nVol < 1), fprintf('%s quitting: you must retain at least one volume\n',mfilename); return; end;
cropName4D = fullfile(pth,['c', nam, ext]); %remove volume is 'img.nii,1' -> 'img.nii'
hdr = spm_vol([imgName4D]); 
[img] = spm_read_vols(hdr);
[nX nY nZ nV] = size(img);
if (skipVol < 0), skipVol = 0; end;
if ( (skipVol+nVol) > nV), nVol = nV - skipVol; end;
if ((skipVol == 0) & (nVol == nV)), fprintf('%s quitting: image only has %d volumes\n',mfilename,nV); return; end;
fprintf('%s has volumes %d..%d volumes from %s\n',cropName4D,skipVol+1,skipVol+nVol,imgName4D); 
hdr = hdr(1);
hdr.fname   = cropName4D;
nslices = hdr.dim(3);
hdr.dim(3) = nslices;
for vol=1:nVol
    hdr.n(1)=vol;
    imgMod = img(:, :, (1:nslices), skipVol+vol);
    spm_write_vol(hdr,imgMod(:, :, :, 1));
end;
