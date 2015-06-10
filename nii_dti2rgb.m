function nii_dti2rgb (fnmFA, fnmV1, old_RGB)
%Convert Fractional Anisotropy and Primary Vector maps to RGB
% fnmFA (optional) name of FA image
% fnmV1 (optional) name of V1 image - must be 4d volume with 3 directions
% old_RGB (optional) support Analyze planar RGB instead of NIfTI RGB
%For comments on old RGB see
% http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image/content/save_nii.m
%Example
% dti2rgb('dti_FA.nii.gz');

fnmFA = 'FA.nii.gz';
if ~exist('fnmFA','var')  %fnmFA not specified
   [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select FA image');
   fnmFA = [Apth, A];
end;
if ~exist('fnmV1','var') || ~exist(fnmV1,'file') %fnmFA not specified
   fnmV1 = regexprep(fnmFA,'FA.','V1.');
   if ~exist(fnmV1,'file')
        [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select V1 image');
        fnmV1 = [Apth, A];
   end
end;
if ~exist('old_RGB','var')  %fnmFA not specified
   old_RGB = true;%false;
end;
[hdrFA, imgFA] = nii_loadhdrimg(fnmFA);
[hdrV1, imgV1] = nii_loadhdrimg(fnmV1); %#ok<ASGLU>
if (size(imgV1,4) < 3) 
	error('V1 image must have at least 3 volumes');
end
%normalize images from 0..1
imgFA = (imgFA - min(imgFA(:)))/max(imgFA(:));
imgV1 = abs(imgV1);
%imgV1 = (imgV1 - min(imgV1(:)))/max(imgV1(:));
%create images
hdr = hdrFA;
i = 1;
imgFA(:) = imgFA(:) * 255;
if old_RGB %planar...
    imgFA = reshape(imgFA,hdr.dim(1)*hdr.dim(2),hdr.dim(3));
    imgV1 = reshape(imgV1,hdr.dim(1)*hdr.dim(2),hdr.dim(3),3);    
    img = zeros(hdr.dim(1)*hdr.dim(2),hdr.dim(3)*3);
    for s = 1: hdr.dim(3) %for each slice
     for rgb = 1:3 % for red,green,blue
        img(:,i) =imgV1(:,s,rgb) .* imgFA(:,s);
        i = i + 1;
     end
    end
else %triples rgbrgb...
    imgFA = repmat(imgFA(:),1,3);
    imgFA = imgFA';
    imgV1 = reshape(imgV1,hdr.dim(1)*hdr.dim(2)*hdr.dim(3),3);
    imgV1 = permute(imgV1,[2,1]);
    img = imgFA(:) .* imgV1(:);
end
hdr.dt = 128;
img = reshape(img,hdr.dim(1),hdr.dim(2),3*hdr.dim(3));
[pth nam ext] = fileparts(fnmFA);
fnmRgb = fullfile(pth, ['rgb' nam ext]);
nii_savehdrimg(fnmRgb, hdr,img);
%end nii_dti2rgb()
