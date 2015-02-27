function nii_combine_mask(fnms, filt, isAndMask) 
%Input: multiple images combined using AND or OR logical functions
%  fnms: filenames of input images
%  filt: if zero all non-zero voxels, if -1 then negative voxels, if +1 then all positive voxels
%  isAndMask : if 1 AND mask (voxel must survive ALL), else OR mask (voxel must survive ANY) 
% Outputs: image with 'm' appended to first filename
%Examples
%  nii_or_mask;
%  nii_or_mask(strvcat('A.nii','B.nii'),0, 1);
%Chris Rorden (2014)

if ~exist('fnms','var')
    fnms = spm_select(inf,'image','Select images to be compared');
end
if  ~exist('isAndMask','var')  
    prompt = {'Voxel threshold (0=nonzero, -1=negative,+1=positive','Find voxels common to all (1 AND) or any (0 OR)'};
    dlg_title = 'Preserve which intensities?';
    num_lines = 1;
    def = {'0','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    filt = str2double(answer{1});
    isAndMask = str2double(answer{2});
end;
%read first image
hdr1 = spm_vol(deblank(fnms(1,:)));
oimg = spm_read_vols(hdr1);
oimg = filtSub(oimg, filt); %binarize
for i=2:size(fnms,1)
    hdr = spm_vol(deblank(fnms(i,:)));
 	img = spm_read_vols(hdr);
    img = filtSub(img, filt); %binarize
    if isAndMask
        oimg = oimg & img; %combine with logical AND        
    else
        oimg = oimg | img; %combine with logical OR
    end
end;
if (sum(oimg(:)>0) == 0)
    fprintf('Error: no voxels in any image survives threshold\n');
    return
end
%save image
[pth,nam,ext] = fileparts(hdr1.fname);
hdr1.fname  = fullfile(pth,['m' nam ext]);
fprintf('%d voxels survive in %s\n',sum(oimg(:)>0),hdr1.fname);
hdr1 = spm_create_vol(hdr1); %save header to disk
spm_write_vol(hdr1,oimg); %save image to disk
%end nii_combine_mask()

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
    