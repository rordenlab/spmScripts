function nii_fixnames(dir)
%MUSC uses non-standard names for sequences, fix them!
% dir: folder with NIfTI images to rename

dir = pwd;
modalityKeysOld = {'CRH_ASL', 'mb_diff_PA_','mb_diff_', 'RESTING_STATE_PA_', 'RESTING_STATE_'};
modalityKeysNew = {'ASL', 'DTIrev','DTI_', 'RestRev_', 'Rest_'}; 
nameFiles=subImgSub(dir);

for i = 1: numel(modalityKeysOld)
    keyOld = char(modalityKeysOld(i));
    keyNew = char(modalityKeysNew(i));
    for j = 1: numel(nameFiles)
        inname = char(nameFiles(j));
        if ~exist(inname, 'file'), continue; end; %renamed by previous key
        if strncmpi(keyOld, inname, numel(keyOld))
           %fprintf('%d %s -> %s %s\n', i, (keyOld), (keyNew), char(nameFiles(j)) );
           [p,n,x] = fileparts(inname);
           n(1:numel(keyOld)) = '';
           n = [keyNew, n];
           outname = fullfile(p,[n,x]);
           
           moveImgUnGz(inname, outname);
           %break;
        end
    end
end;
%nii_fixnames()

function moveImgUnGz(inname, outname)
[ipth, inam,iext] = fileparts(inname);
[opth, onam,oext] = fileparts(outname);
if strcmpi(oext,'.gz') %unzip compressed data
    [~, onam,oext] = fileparts(onam);
    outname = fullfile(opth, [onam, oext]);
end
%load data
if strcmpi(iext,'.gz') %unzip compressed data
    gzname = inname;
	inname = gunzip(inname);
    inname = deblank(char(inname));
    [ipth, inam,iext] = fileparts(inname);
    movefile(inname, outname);
    delete(gzname);
else
    movefile(inname, outname);
end;
fprintf('%s -> %s\n', [inam, iext], [onam, oext] );

%copy bvec
ibvec = fullfile(ipth, [inam, '.bvec']);
if exist(ibvec, 'file'),
    obvec = fullfile(opth, [onam, '.bvec']);
    movefile(ibvec, obvec);
end;
%copy bval
ibval = fullfile(ipth, [inam, '.bval']);
if exist(ibval, 'file'),
    obval = fullfile(opth, [onam, '.bval']);
    movefile(ibval, obval);
end;
%end moveImgUnGz()

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