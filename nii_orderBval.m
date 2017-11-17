function nii_orderBval (bval)
%DKE requires DWI volumes sorted by bvalue
%  bval: name of bval file to reorder (assumes file.bval, file.bvec,file.nii)
% Examples
%   nii_orderBval('DWI_dir42_AP_27_EP.bval');

if ~exist('bval','var')
 bval = spm_select(1,'^.*\.(bval)$','Select b-value file to re-order');
end;
if isempty(bval) || ~exist(bval,'file'), error('Unable to find bval file'); end
%load data
[hd,im, bve,bva] = loadSub(bval);
%sort data
if issorted(bva), warning('Already sorted "%s"', bval); return; end;
[sbva, idx] = sort(bva);
sbve = bve(idx, :);
dims = size(im);
im = reshape(im, prod(dims(1:3)), dims(4));
sim = im(:, idx);
sim = reshape(sim, dims);
%save data
[p,n,x] = fileparts(bval);
n = ['o', n]; %ordered
dlmwrite(fullfile(p,[ n, '.bval']), sbva','delimiter','\t');
dlmwrite(fullfile(p,[n, '.bvec']), sbve','delimiter','\t');
hd.fname = fullfile(p,[n,'.nii']);
for vol=1:dims(4)
    hd.n(1)=vol;
    spm_write_vol(hd,sim(:, :, :, vol));
end;
%end nii_orderBval()

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
bve = loadTxt(fullfile(p,[n,'.bvec']), 3);
if isempty(bve), bve = zeros( size(im,4),3); end;
bva = loadTxt(fullfile(p,[n,'.bval']));
if isempty(bva), bva = zeros( size(im,4),1); end;
%end loadSub()

function txt = loadTxt(fnm, nRow)
txt = [];
%fprintf('%s\n', fnm);
if exist(fnm,'file') == 0, return; end; 
fileID = fopen(fnm);
txt = cell2mat( textscan(fileID,'%f'));
fclose(fileID);
if exist('nRow', 'var')
    txt = reshape(txt, [numel(txt)/nRow nRow]);
end
%end loadTxt()