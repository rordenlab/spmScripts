function nii_segment4d (image, n)
%convert a 4D image with order 1,2..n,1,2..n,1,2,--n into N separate images
% image : 4D image to segment
% n : number of interleaved timepoints, e.g. echoes
%Example
% nii_segment4d('img.nii', 3); %e.g. make 3 images 1st will have volumes 1,4,7...

if ~exist('image','var') || isempty(image) %no files specified
    image = spm_select(1,'image','Select 4D volume');
end
if ~exist('n','var')
	answer = inputdlg('Scaling factor ("Number of interleaved images' ,'Settings',1,{'2'});
	n=str2num(answer{1}); %#ok<ST2NM>
end
n = round(n);
if (n < 2), error('n must be 2 or larger'); end
[pth,nam,ext] = spm_fileparts(image);
image = fullfile(pth,[nam,ext]); %'img.nii,1' -> 'img.nii'
hdr = spm_vol(image);
img = spm_read_vols(hdr);
nvol = numel(hdr);
if rem(nvol, n) ~= 0, error('Number of volumes (%d) must be evenly divisible by n (%d)', nvol, n); end
hdr = hdr(1);
for i = 1 : n 
    hdr.fname = fullfile(pth, ['n' num2str(i) '_' nam ext]);
    k = i;
    for vol=1: (nvol / n)
        hdr.n(1)=vol;
        spm_write_vol(hdr,img(:, :, :, k));
        k = i + n;
    end;
end
%end nii_segment4d()