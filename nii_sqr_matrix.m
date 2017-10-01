function [hd,im] = nii_sqr_matrix(fnm)
%convert an image to have a square matrix in plane, e.g. 80x128 saved as 128x128
% fnm : NIfTI image with rectangular matrix (number of columns ~= number of rows)
%Examples
% nii_sqr_matrix('DTI_dir42_AP_M2029_POLAR1012_Session3.nii');
% nii_sqr_matrix; %

if ~exist('fnm','var') || ~exist(fnm,'file') %no files
    fnm = spm_select(1,'^.*\.(gz|nii)$','Select rectangular NIfTI');
end;
[hd,im] = loadSub(fnm);
if hd.dim(1) == hd.dim(2)
   fprintf('Matrix already square %dx%d',  hd.dim(1), hd.dim(2));
   return;
end
mx = max(hd.dim(1),hd.dim(2));
mn = min(hd.dim(1),hd.dim(2));
margin1 = round(0.5* (mx-mn));
margin2 = (mx-mn) - margin1;
if hd.dim(1) > hd.dim(2) %add rows
    m = zeros(size(im,1), margin1, size(im,3), size(im,4));
    im = [im m];
    m = zeros(size(im,1), margin2, size(im,3), size(im,4));
    im = [m im];
else
    m = zeros(margin1, size(im,2), size(im,3), size(im,4));
    im = [im; m];
    m = zeros(margin2, size(im,2), size(im,3), size(im,4));
    im = [m; im];    
end
warning('The origin will have shifted: coregistration will be off');
[p,n,x] = spm_fileparts(fnm);
movefile(fnm, fullfile(p,['rect_' n, x]));
%[p,n,x] = spm_fileparts(hd.fname);
%hd.fname = fullfile(p,['c' n, x]);
hd.fname
hd.dim(1) = size(im,1);
hd.dim(2) = size(im,2);
for vol=1:size(im,4)
    hd.n(1)=vol;
    spm_write_vol(hd,im(:, :, :, vol));
end;
%end nii_sqr_matrix()

function [hd,im, bve,bva] = loadSub(fnm)
im = []; bve = []; bva = [];
[p,n,x] = spm_fileparts(fnm);
if (strcmpi(x,'.bvec')) || (strcmpi(x,'.bval'))
    fnm = fullfile(p,[n,'.nii']);
    if ~exist(fnm,'file')
        fnm = fullfile(p,[n,'.nii.gz']);
    end
    [p,n,x] = spm_fileparts(fnm);
end;

if (length(x)==3)  && min((x=='.gz')==1) 
    fnm = char(gunzip(fnm));
    delnam = fnm;
    [p,n,x] = spm_fileparts(char(fnm));
else
    delnam = '';
end
hd = spm_vol(fnm); %input header
if hd(1).dt(1) == 128
    fprintf('Warning: skipping RGB image %s\n', fnm);
    if ~isempty(delnam), delete(delnam); end;
    return;
end
im = spm_read_vols(hd);%Input image
hd = hd(1);
if ~isempty(delnam), delete(delnam); end;
%end loadSub()
