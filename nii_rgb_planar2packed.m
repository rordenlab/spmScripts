function nii_rgb_planar2packed(fnm)
%for 24-bit (Red/Green/Blue) images, convert between planar (Analyze) and triplet (NIfTI) image storage
%Analyze stores data as 2D planes rrrr..rggggg...gbbbb...b, NIfTI as packed triplets rgbrgbrgb...
%This script detects input images format and converts to the other style
% fnm : name of image to convert
%Examples
% nii_rgb_planar2packed('analyzeRGB.nii');

if ~exist('fnm','var') %file name not specified
    [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select Red-Green-Blue image to convert');
    fnm = [Apth, A];
end
[hdr, imgRGB] = nii_loadhdrimg(fnm);
if (hdr.dt(1) ~= 128), error('This script is for RGB images'); end;
[pth nam ext] = fileparts(fnm);
nX = hdr.dim(1); 
nY = hdr.dim(2); %pixels per slice
nZ = hdr.dim(3); %number of slices
isplanar = packedDx(imgRGB, nX, nY, nZ) > planarDx(imgRGB, nX, nY, nZ);
if isplanar
    fprintf('Converting from planar (Analyze) to packed triplets (NIfTI)\n');
    imgPacked = toPacked(imgRGB, nX, nY, nZ);
    fnm = fullfile(pth, ['a2n_' nam ext]);
    nii_savehdrimg(fnm, hdr,imgPacked);
else
    fprintf('Converting from packed triplets (NIfTI) to planar (Analyze)\n');
    imgPlanar = toPlanar(imgRGB, nX, nY, nZ); 
    fnm = fullfile(pth, ['n2a_' nam ext]);
    nii_savehdrimg(fnm, hdr,imgPlanar);
end
%nii_rgb_planar2packed()

function dx = packedDx(imgRGB, nX, nY, nZ)
%how similar are voxels to their neighbor assuming rgbrgbrgb...
imgRGB = reshape(imgRGB,nX*3,nY, nZ);
imgRGB2 = imgRGB([4:end, 1:3],:,:);
dx = sum(abs(imgRGB(:) - imgRGB2(:)));
%packedDx()

function dx = planarDx(imgRGB, nX, nY, nZ)
%how similar are voxels to their neighbor  rrr...rggg...gbbbb...b
imgRGB = reshape(imgRGB,nX*nY,3, nZ);
imgRGB2 = imgRGB([2:end, 1],:,:);
dx = sum(abs(imgRGB(:) - imgRGB2(:)));
%planarDx()

function imgRGB = toPacked(imgRGB, nX, nY, nZ)
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
%end toPacked()

function imgRGB = toPlanar(imgRGB, nX, nY, nZ)
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
%end toPlanar()