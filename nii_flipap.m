function nii_flipap (V)
% Generates anterior-posterior mirrored image (flip on Y dimension)
%  V : image(s) to flip
% Examples
%  nii_flipap('mri.nii');
%  nii_flipap(strvcat('C:\dir\p1.nii','C:\dir\p2.nii')); %scans from participant 1 and 2

if nargin <1 %no files specified
 V = spm_select(inf,'image','Select image to left-right flip');
end
fprintf('%s: Flipping order of rows - please make sure this is the left-right axis\n', mfilename);
for i=1:size(V,1)
    ref = deblank(V(i,:));

    [pth,nam,ext] = spm_fileparts(ref);
    hdr  = spm_vol(fullfile(pth,[nam ext]));
    img = spm_read_vols(hdr);
    imgflip = flipdim(img,2);
    nvol = numel(hdr);
    hdr = hdr(1);
    for vol = 1 : nvol
        hdr.n(1)=vol;
        hdr.fname = fullfile(pth,['AP' nam ext]);
        spm_write_vol(hdr,imgflip(:,:,:,vol));
    end;
end; %for each image
%end nii_flipap()