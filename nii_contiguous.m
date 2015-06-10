function nii_contiguous(fnm, conn)
%Given a region of interest map, generate map where each region is contiguous
% fnm : (optional) region of interest map
% conn : connectivity criterion: 6 (face), 18 (edge) or 26 (corner)
%Examples
% nii_contiguous; %use GUI
% nii_contiguous('img.nii');

if ~exist('fnm','var')
	fnm = spm_select(1,'image','Select region of interest map'); 
end
if ~exist('conn','var')
	conn = 6; %require voxels share face to be considered contiguous
end
tic;
hdr = spm_vol(fnm);
hdr = hdr(1); %just in case a 4D image is selected
imgIn = spm_read_vols(hdr);
nroiIn = 0;
nroiOut = 0;
imgOut = zeros(size(imgIn));
for i = 1: max(imgIn(:))
    imgi = double(imgIn > (i-0.5) & imgIn < (i+0.5));
    nvox = sum(imgi(:));
    if nvox < 1, continue; end;
    [L,num] = spm_bwlabel(imgi, conn); %use default 18-neighbor definition of cluster
    fprintf('Input region %d has %d voxels and %d contiguous regions\n',i, nvox, num);
    for r = 1 : num
       imgOut(L == r) = nroiOut + r; 
    end
    nroiIn = nroiIn + 1;
    nroiOut = nroiOut + num; 
end
fprintf('Input had %d regions, output has %d contiguous regions (required %gs)\n', nroiIn, nroiOut, toc);
%save output
[pth nm ext] = spm_fileparts(fnm);
hdr.fname = fullfile(pth, ['z' nm ext]);  
spm_write_vol(hdr,imgOut);
%end nii_contiguous()
