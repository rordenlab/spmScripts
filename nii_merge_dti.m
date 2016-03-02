function nii_merge_dti(fnms, minVol)
%Merge a set of DTI scans. Assumes img.nii has img.bvec/img.bval
%If no bvec/bval file is found it is assumed that this is a B0 series
% fnms: filenames to merge
% minVol: only add images with this many volumes (exclude pre-computed MD/ADC/trace)
%Examples
% nii_merge_dti; %use GUI;
% nii_merge_dti(strvcat('98_AP_3.nii.gz','99_AP_7.nii.gz'));
% nii_merge_dti(strvcat('98_AP_3.nii.gz','99_AP_7.nii.gz'), 2);%at least 2 vols
if nargin < 1,
   [files,pth] = uigetfile({'*.gz;*.nii;*.hdr;';'*.*'},'Choose positive DTI[s]', 'MultiSelect', 'on');
   fnms = strcat(pth,char(files));
end
if nargin < 2,
   minVol = 4;
end
if size(fnms,1) < 2
   error('%s requires more than one image', mfilename);
end
bvalCat = [];
bvecCat = [];
imgCat = [];
hdrOK = [];
for i=1:size(fnms,1)
    fnm = deblank(fnms(i,:));
    [hdr, img] = readNiftiSub(fnm);
    if isempty(hdr), continue; end; %unable to read file
    nVol = size(img,4);
    if nVol < minVol
        fprintf('Excluding image (only %d volumes): %s\n', nVol, fnm);
    else
        fprintf('Including %s %dx %dx %dx %d\n', fnm, size(img,1), size(img,2), size(img,3), size(img,4));
        imgCat = cat(4,imgCat, img);
        [bval, bvec] = getBvalBvecSub (fnm);
        if isempty(bval), bval = zeros(nVol,1); end;
        if isempty(bvec), bvec = zeros(nVol,3); end;
        bvecCat = [bvecCat; bvec];
        bvalCat = [bvalCat; bval];
        hdrOK = hdr;
    end
end
%save data
nV = size(imgCat,4);
if (nV < 2) || isempty(hdrOK), error('No volumes to merge'); end;
hdr = hdrOK;
[p,n] = fsl_filepartsSub(fnm);
n = [int2str(nV),'_', n];
hdr = hdr(1);
hdr.fname = fullfile(p, ['DTI_', n,'.nii']);
for vol=1:nV
    hdr.n(1)=vol;
    spm_write_vol(hdr,imgCat(:, :, :, vol));
end;
if (max(bvecCat) == min(bvecCat))
    fprintf('Warning: bvec/bval files will not be created (no variability in b-values)');
    return;
end;
dlmwrite(fullfile(p,['DTI_',n, '.bvec']),bvecCat','delimiter','\t');
dlmwrite(fullfile(p,['DTI_',n, '.bval']),bvalCat','delimiter','\t');
%end nii_merge_dti

function [hdr, img] = readNiftiSub(filename)
%load NIfTI (.nii, .nii.gz, .hdr/.img) image and header
% filename: image to open
% open4d: if true all volumes are loaded
%To do:
%  endian: rare, currently detected and reported but not handled
%Examples
% hdr = nii_loadhdrimg('myimg.nii');
% [hdr, img] = nii_loadhdrimg('myimg.nii');
% [hdr, img] = nii_loadhdrimg('img4d.nii');
if ~exist('filename','var')  %fnmFA not specified
   [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select image');
   filename = [Apth, A];
end
[fpth, fnam,fext] = fileparts(filename);
if strcmpi(fext,'.img') %hdr/img pair
    filename = fullfile(fpth, [fnam, '.hdr']);
end
if ~exist(filename, 'file')
    error('Unable to find file %s', filename);
end
%load data
if strcmpi(fext,'.gz') %unzip compressed data
	filename = gunzip(filename);
    filename = deblank(char(filename));
end;
hdr = spm_vol(filename);
if hdr(1).dt(1) == 128
   fprintf('Skipping RGB image %s\n', filename);
   hdr = [];
   img = [];
   return;
end
img = spm_read_vols(hdr);
if strcmpi(fext,'.gz') %fsl can not abide with coexisting img.nii and img.nii.gz
	delete(filename);
end;
%end nii_loadhdrimg()

function [bval, bvec] = getBvalBvecSub (imgName)
%read .bval file and return indices for B0 volumes
bval = [];
bvec = [];
[pth, nam] = fsl_filepartsSub(imgName);
nameVal = fullfile(pth,[nam '.bval']); %name for b-values
nameVec = fullfile(pth,[nam  '.bvec']); %name for b-vectors
if (exist(imgName, 'file') == 0) , fprintf('Unable to find required image %s\n',imgName); return; end;
if ( (exist(nameVal, 'file') == 0) || (exist(nameVec, 'file') == 0) ), fprintf('Unable to find required DTI files %s and %s\n',nameVec,nameVal); return; end;
%bval = importdata(nameVal); %<- does not work with Matlab 2014b on Linux
%read b-values
fileID = fopen(nameVal);
bval = cell2mat( textscan(fileID,'%d'));
fclose(fileID);
%read b-vectors
fileID = fopen(nameVec);
bvec = cell2mat( textscan(fileID,'%f'));
fclose(fileID);
if mod(numel(bvec),3) ~= 0
    error('Error: number of bvecs must be divisible by three. Found %d bvecs in %s', numel(bvec), nameVec);
end
bvec = reshape(bvec,numel(bvec)/3,3);
%end getB0vols()

function [pth nam ext] = fsl_filepartsSub(fileName)
% a.nii.gz has the extension ".nii.gz" not ".nii"
[pth nam ext] = fileparts(fileName);
if (length(ext)==3)  && min((ext=='.gz')==1)
	[pth nam ext2] = fileparts( fullfile(pth, nam)); %remove .nii for .nii.gz
    ext = [ext2 ext];
end;
%end fsl_filepartsSub()




