function nii_xor_mask(fnms, filt) 
%Exclusive or masking, Input: multiple images, output image: voxels are unique to each input
%  fnms: pair of input images
%  filt: if zero all non-zero voxels, if -1 then negative voxels, if +1 then all positive voxels
% Outputs: image with name1+'not'+name2 and name2+'not'+name1
%Examples
%  nii_or_mask;
%  nii_or_mask(strvcat('A.nii','B.nii'),0);
%Chris Rorden (2014)

if ~exist('fnms','var')
    fnms = spm_select(2,'image','Select source image that will be resliced');
end
if size(fnms,1) ~= 2
    error('You need to select 2 images, you selected %d',size(fnms,1));
end
if  ~exist('filt','var')  
    prompt = {'Voxel threshold (0=nonzero, -1=negative,+1=positive'};
    dlg_title = 'Preserve which intensities?';
    num_lines = 1;
    def = {'0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    filt = str2double(answer{1});
end;
%read first image
hdr1 = spm_vol(deblank(fnms(1,:)));
img1 = spm_read_vols(hdr1);
img1 = filtSub(img1, filt); %binarize
[pth1,nam1,ext1] = fileparts(hdr1.fname);
%read second image
hdr2 = spm_vol(deblank(fnms(2,:)));
img2 = spm_read_vols(hdr2);
img2 = filtSub(img2, filt); %binarize
[pth2,nam2,ext2] = fileparts(hdr2.fname);
%compute voxels exclusive to img1
oimg = img1;
oimg(img2==1) = 0;
if (sum(oimg(:)>0) == 0)
    fprintf('No voxels survive only in %s\n',hdr1.fname);
else
    fprintf('%d voxels survive only in %s\n',sum(oimg(:)>0),hdr1.fname);
    hdr1.fname  = fullfile(pth1,[nam1 '_not_' nam2 ext1]);
    hdr1 = spm_create_vol(hdr1); %save header to disk
    spm_write_vol(hdr1,oimg); %save image to disk
end
%compute voxels exclusive to img2
oimg = img2;
oimg(img1==1) = 0;
if (sum(oimg(:)>0) == 0)
    fprintf('No voxels survive only in %s\n', hdr2.fname);
else
    fprintf('%d voxels survive only in %s\n',sum(oimg(:)>0),hdr2.fname);
    hdr2.fname  = fullfile(pth2,[nam2 '_not_' nam1 ext2]);
    hdr2 = spm_create_vol(hdr2); %save header to disk
    spm_write_vol(hdr2,oimg); %save image to disk
end
%end nii_xor_mask(), local functions follow

function i = filtSub(i, fil)
if (fil < 0) 
    i(i > 0) = 0;
    i(i < 0) = 1;
elseif (fil > 0)
    i(i > 0) = 1;
    i(i < 0) = 0;    
else
   i(i ~= 0) = 1; %make all non-zero voxels = 1
end
%end filtSub()
    