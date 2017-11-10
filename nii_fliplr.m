function nii_fliplr (V)
% Generates left-right mirror image
%   nii_mirror('C:\dir\img.nii');
% You can also pass multiple images - each will be mirrored
% Example - scans from participant 1 and 2
%  nii_mirror(strvcat('C:\dir\p1.nii','C:\dir\p2.nii'));

if nargin <1 %no files specified
 V = spm_select(inf,'image','Select image to left-right flip');
end
fprintf('%s: Flipping order of rows - please make sure this is the left-right axis\n', mfilename);
for i=1:size(V,1)
    ref = deblank(V(i,:));
    
    [pth,nam,ext] = spm_fileparts(ref); 
    hdr  = spm_vol(fullfile(pth,[nam ext]));
    img = spm_read_vols(hdr);
    imgflip = flipdim(img,1);
    nvol = numel(hdr);
    hdr = hdr(1);
    for vol = 1 : nvol
        hdr.n(1)=vol;
        hdr.fname = fullfile(pth,['RL' nam ext]);
        spm_write_vol(hdr,imgflip(:,:,:,vol));
    end;
end; %for each image

% for i=1:size(V,1)
%     ref = deblank(V(i,:));
%     hdr.n(1)=vol;
%     [pth,nam,ext] = spm_fileparts(ref); 
%     hdr  = spm_vol(fullfile(pth,[nam ext ',1']));
%     img = spm_read_vols(hdr);
%     imgflip = flipdim(img,1);
%     hdr.fname = fullfile(pth,['RL' nam ext]);
%     spm_write_vol(hdr,imgflip);
% end; %for each image
%end nii_fliplr()
