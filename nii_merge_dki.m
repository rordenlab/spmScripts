function nii_merge_dki (V)
%merge a series of DTI/DKI images for analysis
%  V: image(s) to merge
% Examples
%   nii_merge_dki('brain.nii.gz');

minVol = 3; %1 to merge all images, 3 to merge images with at least 3 volumes
if nargin <1 %no files
 V = spm_select(inf,'^.*\.(gz|nii)$','Select gz files to decompress');
end;
if ischar(V), V = cellstr(V); end;
num = numel(V);
hdr = []; img = []; bvec = []; bval = []; %vol = [];
if num < 2, error('You must specify multiple files'); end;
for n = 1: num
    [hd,im, bve,bva] = loadSub(char(V(n)));
    if (size(im,4) < minVol), continue; end;
    if isempty(im), continue; end;
    if ~isempty(img)
        if (size(im,1) ~= size(img,1)) || (size(im,2) ~= size(img,2)) || (size(im,3) ~= size(img,3))
            fprintf('WARNING: you are attempting to merge images with different dimensions\n');
        end
    end
    bvec = [bvec; bve];
    bval = [bval; bva];
    %vol = [vol; size(im,4)];
    hdr = [hdr; hd]; %#ok<*AGROW>
    img = cat(4,img, im);
    
    
end

if numel(hdr) < 2, error('Unable to load more than one image'); end;
hdrOut = hdr(1);
[p,n,x] = spm_fileparts(hdrOut.fname);
hdrOut.fname = fullfile(p,['m', n, x]);
for vol=1:size(img,4)
    hdrOut.n(1)=vol;
    spm_write_vol(hdrOut,img(:, :, :, vol));
end;
if abs(max(bvec(:))) == 0, return; end;
dlmwrite(fullfile(p,['m', n, '.bval']), bval','delimiter','\t');
dlmwrite(fullfile(p,['m', n, '.bvec']), bvec','delimiter','\t');
plotBvecSub(fullfile(p,['m', n, '.bvec']));
%end nii_merge_dki()

function plotBvecSub(fnm)
bvecs = load(fnm); % Assuming your filename is bvecs
figure('position',[100 100 500 500]);
plot3(bvecs(1,:),bvecs(2,:),bvecs(3,:),'*r');
axis([-1 1 -1 1 -1 1]);
axis vis3d;
rotate3d
%end plotBvecSub()

function [hd,im, bve,bva] = loadSub(fnm);
[p,n,x] = spm_fileparts(fnm);
if (length(x)==3)  && min((x=='.gz')==1) 
    fnm = char(gunzip(fnm));
    delnam = fnm;
    [p,n,x] = spm_fileparts(char(fnm));
else
    delnam = '';
end
hd = spm_vol(fnm); %input header
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