function nii_harvest

% /Root/subj1/study1/T1.nii
% /Root/subj1/study2/T2.nii
% /Root/subj2/study1/T1.nii
% /Root/subj2/study2/T2.nii
% /Root/subj2/study2/ASL.nii

outDir = '/Users/rorden/Desktop/pp';
%baseDir = '/Users/rorden/Desktop/harvest'; %'/Root'
if ~exist('baseDir','var') || isempty(baseDir)
    baseDir = pwd; %uigetdir('','Pick folder that contains all subjects');
end
subjDirs = subFolderSub(baseDir);
subjDirs = sort(subjDirs);
modalityKeys = {'T1', 'T2', 'Lesion', 'ASL', 'DTIrev','DTI', 'Rest', 'fMRI'}; %DTIREV before DTI!!! both "DTIREV.nii" and "DTI.nii" have prefix "DTI"
xperimentKeys = {'LIME', 'CT', 'R01'}; %order specifies priority: 1st item checked first!
%create empty structure
blank = [];
blank.subjName = [];
for i = 1: numel(modalityKeys)
    blank.nii.(modalityKeys{i}) =[];
end;
%1st: acquire data
nSubj = 0;
for s = 1: size(subjDirs,1)%1:nSubjDir2 %(nSubjDir2+1):nSubjDir  
    subjName = deblank(subjDirs{s});
    if subjName(1) == '.', continue; end;
    if (numel(subjName) > 1) && (subjName(2) == '4'), fprintf('SKIPPING %s\n', subjName); continue; end; %ignore folders with underscore, "M2015_needsmatfile"
    
    if isStringInKeySub (subjName,'_'), continue; end; %ignore folders with underscore, "M2015_needsmatfile"
    subjDir = [baseDir,filesep, subjName]; %no filesep
    %fprintf('%s\n', subjDir);
    nSubj = nSubj + 1;
    imgs(nSubj) = blank;
    imgs(nSubj).subjName = subjName;
    for x = 1: numel(xperimentKeys)
        xLabel = deblank(xperimentKeys{x}); %e.g. "R01"
        xDir = [subjDir,filesep, xLabel]; %no filesep
        if ~exist(xDir, 'file'), continue; end;
        %fprintf('%s\n', xDir);
        imgs(nSubj) = findImgsSub(imgs(nSubj), xDir, xLabel);
        
    end
    %imgs(nSubj).nii

end
fprintf('Found %d subjects in %s\n', nSubj, baseDir);
if nSubj < 1, return; end;
%report results
% 1st row: describe values
f = fieldnames(imgs(1).nii);
str = 'subj';
for i = 1: numel(f)
   str = sprintf('%s\t%s',str, f{i} );
end
fprintf('%s\n', str);
% subsequent rows: source of images
for s = 1: nSubj 
    str = imgs(s).subjName;
    for i = 1: numel(f)
        x = '';
        if ~isempty(imgs(s).nii.(f{i}))
           x = imgs(s).nii.(f{i}).x;
        end
        str = sprintf('%s\t%s',str, x );
    end
    fprintf('%s\n', str);
end
%copy core files to new folder
return
if exist(outDir, 'file') ~= 7, error('Unable to find folder %s', outDir); end;
for s = 1: nSubj 
    subj = deblank(imgs(s).subjName);
    subjDir = fullfile(outDir, subj);
    mat = [];
    if exist(subjDir,'file') == 0, mkdir(subjDir); end;
    for i = 1: numel(f)
        if ~isempty(imgs(s).nii.(f{i}))
            m = f{i}; % modality: T1, T2..
            x = imgs(s).nii.(f{i}).x; %e.g. experiment name "LIME", "CT" 
            imgin = imgs(s).nii.(f{i}).img; %e.g. '~/dir/m2000/CT/T1.nii'
            imgout = fullfile(subjDir, sprintf('%s_%s_%s.nii',m, subj, x));
            fprintf('%s -> %s\n',imgin, imgout);
            moveImgUnGz(imgin, imgout);
            mat.(m) = imgout;
        end
    end
    if ~isempty(mat)
        matName = fullfile(subjDir, [subj, '_limegui.mat']);
        fprintf('Creating %s\n',matName);
        save(matName,'-struct', 'mat');
        setAcpcSub(matName);
    end
    
end
%end nii_harvest

function setAcpcSub (matname)
m = load(matname);
if isfield(m,'T1')
    nii_setOrigin12(m.T1, 1, true); %T1 - crop
end
if isfield(m,'T2') && isfield(m,'Lesion')
    nii_setOrigin12({m.T2,m.Lesion}, 2, true); %T2
end
if isfield(m,'DTI') && isfield(m,'DTIrev')
    nii_setOrigin12({m.DTI,m.DTIrev}, 3, false); %DTI
elseif isfield(m,'DTI')
    nii_setOrigin12(m.DTI, 3, false); %DTI
end
%end setAcpcSub();

function moveImgUnGz(inname, outname)
[ipth, inam,iext] = fileparts(inname);
[opth, onam,oext] = fileparts(outname);
%load data
if strcmpi(iext,'.gz') %unzip compressed data
	inname = gunzip(inname);
    inname = deblank(char(inname));
    [ipth, inam,iext] = fileparts(inname);
end;
copyfile(inname, outname);
if strcmpi(iext,'.gz') %fsl can not abide with coexisting img.nii and img.nii.gz
	delete(filename);
end;
%copy bvec
ibvec = fullfile(ipth, [inam, '.bvec']);
if exist(ibvec, 'file'),
    obvec = fullfile(opth, [onam, '.bvec']);
    copyfile(ibvec, obvec);
end;
%copy bval
ibval = fullfile(ipth, [inam, '.bval']);
if exist(ibval, 'file'),
    obval = fullfile(opth, [onam, '.bval']);
    copyfile(ibval, obval);
end;
%end moveImgUnGz()


function imgs = findImgsSub(imgs, xDir, xLabel)
nameFiles = subImgSub(xDir);
if isempty(nameFiles), return; end;
%nameFiles = sort(nameFiles); %take first file for multiimage sequences, e.g. ASL
%nameFiles
%nameFiles(1)
f = fieldnames(imgs.nii);
for i = 1: numel(f)
    if ~isempty(imgs.nii.(f{i})), continue; end;
    for j = 1: numel(nameFiles) 
        if strncmpi(f{i},nameFiles(j),numel(f{i}))
           fname = fullfile(xDir, char(nameFiles(j)) );
           imgs.nii.(f{i}).x = xLabel;
           imgs.nii.(f{i}).img = fname;
           %if strncmpi(f{i},'DTI',numel(f{i}))
           % fprintf('%d\t%s\n', nVolSub(fname), fname);   
           %end
           %fprintf('-%s %s\n', xLabel, fname); 
           break;
        end
    end
end;

function nVol = nVolSub(filename)
%report number of volumes
nVol = 0;
[hdr, img] = readNiftiSub(filename);
if isempty(hdr), return; end;
nVol = numel(hdr);
%end nVolSub

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

%end findImgsSub()

function nameFiles=subImgSub(pathFolder)
nameFiles=subFileSub(pathFolder);
if isempty(nameFiles), return; end;
n = nameFiles; nameFiles = [];
for i = 1: numel(n) 
    [~,~,x] = fileparts(char(deblank(n(i))));
    if ~strncmpi('.gz',x, 3) && ~strncmpi('.nii',x, 4), continue; end;
    nameFiles = [nameFiles; n(i)]; %#ok<AGROW>
end
%end subFileSub()


function nameFiles=subFileSub(pathFolder)
d = dir(pathFolder);
isub = ~[d(:).isdir];
nameFiles = {d(isub).name}';
%end subFileSub()

function isKey = isStringInKeySub (str, imgKey)
isKey = true;
for k = 1 : size(imgKey,1)
    key = deblank(imgKey(k,:));
    pos = strfind(lower(char(str)),lower(key));
    if ~isempty(pos), isKey = pos(1); return; end;
end
isKey = false;
%isStringInKey()

function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()