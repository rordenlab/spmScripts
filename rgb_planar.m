function rgb_planar(fnm, ana2nii)
%convert between planar (Analyze) and triplet (NIfTI) image storage
%Analyze stores data rrrr..rggggg...gbbbb...b, NIfTI rgbrgbrgb
% fnm : name of image to convert
% ana2nii : if true (default) planar to triplet, if false triplet to planar
%Analyze (ana) requires planar data, NIfTI (nii) specificies triplet
%Examples
% rgb_planar('analyzeRGB.nii');
% rgb_planar('niftiRGB.nii', false);

%fnm = 'old.nii.gz'
if ~exist('fnm','var')
    [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select Red-Green-Blue image to convert');
    fnm = [Apth, A];
end
if ~exist('ana2nii','var')
    ana2nii = true;
end
[hdr, imgRGB] = nii_loadhdrimg(fnm);

%hdr = spm_vol(fnm);
if (hdr.dt(1) ~= 128), error('This script is for RGB images'); end;
nX = hdr.dim(1); 
nY = hdr.dim(2); %pixels per slice
nZ = hdr.dim(3); %number of slices
[pth nam ext] = fileparts(fnm);
if ana2nii
    %convert rrrr...rgggg...gbbbbb...b to rgbrgbrgb
    imgRGB = reshape(imgRGB,nX,nY,nZ*3);
    imgR = zeros(nX,nY,nZ);
    imgG = zeros(nX,nY,nZ);
    imgB = zeros(nX,nY,nZ);   
    i = 1;   
    for z = 1: nZ %for each slice
        imgR(:,:,z) = imgRGB(:,:,i);
        i = i + 1;
        imgG(:,:,z) = imgRGB(:,:,i);
        i = i + 1;
        imgB(:,:,z) = imgRGB(:,:,i);
        i = i + 1;
    end
    imgRGB = [imgR(:)'; imgG(:)'; imgB(:)'];
    imgRGB = reshape(imgRGB,nX,nY,nZ,3);
    fnm = fullfile(pth, ['a2n_' nam ext]);
    nii_savehdrimg(fnm, hdr,imgRGB);
else
    %convert rgbrgbrgb to rrrr...rgggg...gbbbbb...b
    imgR = imgRGB(1:3:end);
    imgG = imgRGB(2:3:end);
    imgB = imgRGB(3:3:end);
    imgR = reshape(imgR,nX,nY,nZ);
    imgG = reshape(imgG,nX,nY,nZ);
    imgB = reshape(imgB,nX,nY,nZ);
    imgRGB = reshape(imgRGB,nX,nY,nZ*3);
    sOut = 1;
    for sIn = 1: nZ
       imgRGB(:,:,sOut) = imgR(:,:,sIn);
       sOut = sOut + 1;
       imgRGB(:,:,sOut) = imgG(:,:,sIn);
       sOut = sOut + 1;
       imgRGB(:,:,sOut) = imgB(:,:,sIn);
       sOut = sOut + 1;
    end    
    fnm = fullfile(pth, ['n2a_' nam ext]);
    nii_savehdrimg(fnm, hdr,imgRGB);
end
