function nii_fixslices(fnm, slices)
%interpolate slices corrupted by motion, for interleaved DWI
% fnm: image to repair
% slices : slices for fix
%Examples
% nii_fixslices; %use GUI
% nii_fixslices('YCN0572x.nii', [3, 5, 9, 11, 13, 15, 17]); 

if ~exist('fnm','var')
	fnm = spm_select(1,'image','Select image[s] for slice interpolation'); 
end
if ~exist('slices','var') 
    prompt = {'Slices:'};
    dlg_title = ['Options'];
    num_lines = [1 50];
    def = {'3 5 9 11 13 15'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer), return; end;
    slices = str2num(answer{1});
end;
if isempty(slices) || (numel(slices) < 1), return; end;
hdr = spm_vol(fnm);
imgIn = spm_read_vols(hdr);
img = imgIn;
for s = 1 : numel(slices)
    out = slices(s);
    above = out + 1;
    if above > size(img,3), above = size(img,3) - 1; end;
    below = out - 1; 
    if below < 1, below = 2; end;
    img(:,:,out) = 0.5 * (img(:,:,above) + img(:,:,below));
end
%save images
[pth nm ext] = spm_fileparts(fnm);
hdr.fname = fullfile(pth, ['i' nm ext]);  
spm_write_vol(hdr,img);
%end nii_fixslices()