function dicm_sort(src, outFolder)
%Sort DICOMs from source folder and each series in unique folder of outFolder
% src : input folder with many DICOM images
% outFolder : output folder where images will be saved
%
%Requirements
% Xiangrui Li's dicm2nii  https://github.com/xiangruili/dicm2nii
%Example
% dicm_sort() %use graphical interface
% dicm_sort(pwd,pwd)

src = '~/Downloads/disk/DICOM';
outFolder = '~/Downloads/a2';

if isempty(which('dicm_hdr')), error('dicm2nii required'); end;
if ~exist('src','var') %file not specified
   src = uigetdir(pwd,'Select input folder');
end;
if ~exist('outFolder','var') %file not specified
   outFolder = uigetdir(pwd,'Select output folder');
end;
dicm_sortSub(src, outFolder);
%end dicm_sort()

function dicm_sortSub(src, outFolder)
d = dir(src);
isub = [d(:).isdir];
dirs = {d(isub).name}';
for i = 1 : numel(dirs)
    if isempty(dirs{i}) || (dirs{i}(1) == '.'), continue; end;
    dicm_sortSub(fullfile(src, dirs{i}), outFolder);
end
isub = ~[d(:).isdir];
fnms = {d(isub).name}';
%end subFileSub()
for i = 1 : numel(fnms)
    inName = fullfile(src, fnms{i});
    try
        s = dicm_hdr(inName);
    catch
        fprintf('dicm_hdr unable to read %s\n', inName);
        continue;
    end
    if isempty(s), continue; end;
    if ~isfield(s, 'SeriesNumber') || ~isfinite(s.SeriesNumber), continue; end;
    if ~isfield(s, 'AcquisitionNumber') || ~isfinite(s.AcquisitionNumber), continue; end;
    if ~isfield(s, 'InstanceNumber') || ~isfinite(s.InstanceNumber), continue; end;
    if ~isfinite(s.SeriesNumber), continue; end;
    %outDir = fullfile(outFolder, sprintf('%d_%d',s.SeriesNumber, s.AcquisitionNumber));
    outDir = fullfile(outFolder, sprintf('%d', s.AcquisitionNumber));
    
    if ~exist(outDir,'dir'), mkdir(outDir), end;
    %outName = fullfile(outDir, fnms{i});
    baseName = sprintf('%04d.dcm', s.InstanceNumber);
    if isfield(s, 'MRImageDiffBValueNumber') && isfinite(s.MRImageDiffBValueNumber)
        baseName = sprintf('b%02d_%s', s.MRImageDiffBValueNumber, baseName);
    
    end;
    if isfield(s, 'MRImageGradientOrientationNumber') && isfinite(s.InstanceNumber)
        baseName = sprintf('g%03d_%s', s.MRImageGradientOrientationNumber, baseName);
    end;
    outName = fullfile(outDir, baseName);
    if exist(outName,'file'), 
        fprintf('File already exists %s\n', outName); 
        while exist(outName,'file')
            outName = fullfile(outDir, sprintf('%04d_%04ddup.dcm', s.InstanceNumber, round(rand()*9999)));
        end;
    end;
    copyfile(inName, outName);
end
%end dicm_sortSub()