function fiberQA (baseDir)

if exist('baseDir','var') 
    cd(baseDir);
end


m = dir('*.mat');
if isempty(m), error('Unable to find mat files'); end;
fprintf('Found %d subjects in %s\n',numel(m), pwd);
label = jhuLabelSub;
isGM_R = isGMr_Sub(label);
nROI = sum(isGM_R(:));
msk = triu(ones(nROI,nROI),1);
%fieldname = 'dtifc_jhu'; 
%fieldname = 'dti_jhu'; 
%fieldname = 'dtimx_jhu';
fieldname = 'dtimn_jhu';
mat = [];
nam = [];
matOK = 0;
tic
for s = 1:numel(m)
    snam = m(s).name;
    mt = getMatSub(snam, fieldname, msk, isGM_R);
    if isempty(mt)
        fprintf('Skipping %s\n', snam);
        continue;
    end
    matOK = matOK + 1;
    mat(matOK,:) = mt; 
    nam = [nam, {snam}]; %#ok<*AGROW>
end %for each subject - generate population statistics
if matOK == 0, return; end;
%fprintf('Found %d files with field %s\n',  matOK, fieldname);
mn = mean(mat);
rOK = zeros(matOK,1);
fid = fopen([fieldname '.txt'], 'w');
for s = 1:matOK
    [r] = corrcoef (mn, mat(s,:));
    fprintf('%s\tr=\t%g\n',  nam{s}, r(2,1) );
    fprintf(fid, '%s\tr=\t%g\n',  nam{s}, r(2,1) );
    rOK(s) = r(2,1);
end
fclose(fid);
fprintf('For %d files with %s, r min=\t%g\tmean=\t%g\tmax=\t%g\n',  matOK, fieldname, min(rOK), mean(rOK), max(rOK));
toc
%end fiberQA()


function mt = getMatSub(matname, fieldname, msk, isGM_R)
fld = fieldSub(matname, fieldname);
mt = [];
if  isempty(fld), return; end;
%fprintf('Processing %s\n',  matName);
mt =  fld.r; %189x189 matrix
mt = mt(isGM_R,isGM_R); %53x53 matrix
mt = mt(msk ~=0);
%end 

function label = jhuLabelSub 
pth = which('nii_stat');
[pth] = fileparts (pth);
pth = [pth filesep 'roi' filesep 'jhu.txt'];
if ~exist(pth,'file'), error('Unable to find %s\n',pth); end;
fid = fopen(pth);  % Open file
label=[];
tline = fgetl(fid);
while ischar(tline)
    label=strvcat(label,tline); %#ok<DSTRVCT,REMFF1>
    tline = fgetl(fid);
end
fclose(fid); 
%end labelSub()

function fld = fieldSub(matname, fieldname)
m = load(matname);
if isfield(m, fieldname)
    fld = m.(fieldname);
else
    fld = [];
end
%end fieldSub()

function isGM_R = isGMr_Sub(label)
n = size(label,1);
isGM_R = zeros(n,1);
for i = 1: n
    if ~isempty(strfind(label(i,:), '_L|')) && ~isempty(strfind(label(i,:), '|1'))
        isGM_R(i) = 1;
    end
    
end
isGM_R = logical(isGM_R);
%end isGMr_Sub()
