function nii_merge(vols, outnam)
%merge multiple images together
% vols: name of images to weld together
% outnam: optional output name, e.g. 'DTI_P195.nii'
%Examples
% nii_merge({'DTIA_P195.nii','DTIB_P195.nii'})
% nii_merge({'DTIA_P199.nii','DTIB_P199.nii'},'DTI_P199.nii')
% nii_merge %use GUI
if ~exist('vols','var') || isempty(vols) %no files specified
 vols = spm_select(inf,'image','Merge images');
end
if ischar(vols), vols = cellstr(vols); end
if numel(vols) < 2, error('You must specify at least 2 images!'); end;
for v = 1 : numel(vols) %apply parameters from first session to others
    [pth,nam,ext] = spm_fileparts(deblank(vols{v})); %remove ',1'
    fname = fullfile(pth,[nam ext]); %strip volume label
    hdr = spm_vol(fname); %read header - this time 4D if specified
    img = spm_read_vols(hdr); %load image
    hdr = spm_vol([fname,',1']); %read header of 1st volume
    hdr = hdr(1);
    nameVal = fullfile(pth,[nam '.bval']); %name for b-values
    nameVec = fullfile(pth,[nam  '.bvec']); %name for b-vectors

    val = vReadSub(nameVal, size(img,4), 1);
    vec = vReadSub(nameVec, size(img,4), 3);
    if v == 1
        hdrAll = hdr;
        imgAll = img;
        valAll = val;
        vecAll = vec;
        
    else
        if min(hdrAll.dim(1:3) == hdr.dim(1:3)) == 0
            error('Images must have the same number of columns, rows and slices');
        end
        imgAll = cat(4, img, imgAll);
        valAll = [valAll; val]; %#ok<AGROW>
        vecAll = [vecAll; vec]; %#ok<AGROW>
        
    end
    
end
if exist('outnam', 'var')
    [pth,nam,ext] = spm_fileparts(outnam);     
else
    [pth,nam,ext] = spm_fileparts(hdrAll.fname);
    nam = ['x' nam];
end
hdrAll.fname = fullfile(pth,[ nam ext]); %strip volume label
for vol=1:size(imgAll,4)
    hdrAll.n(1)=vol;
    spm_write_vol(hdrAll,imgAll(:, :, :, vol));
end;
vSaveSub(fullfile(pth,[nam '.bval']), valAll);
vSaveSub(fullfile(pth,[nam '.bvec']), vecAll);
%end nii_merge()

function vSaveSub(nam, v);
if ~any(v(:)), return; end; %all zeros...
dlmwrite(nam,v','delimiter','\t');
%end vSaveSub()

function v = vReadSub(nam, nCol, nRow)
zeros(nCol,nRow);
if ~exist(nam,'file'), return; end;
fileID = fopen(nam);
v = cell2mat( textscan(fileID,'%f'));
fclose(fileID);
v = reshape(v, [numel(v)/nRow nRow]);
%end vecReadSub()