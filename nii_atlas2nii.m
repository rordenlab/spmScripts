function nii_atlas2nii (atlasfnm)
%Convert NIfTI indexed atlas to series of NIfTI images, one per region
% atlasfnm = (optional) filenames to convert, e.g. 'aal.nii.gz'
%Examples
% nii_atlas2nii; %use GUI
% nii_atlas2nii('aal.nii.gz');

if ~exist('atlasfnm','var') %no files specified
 atlasfnm = spm_select(1,'^.*\.(gz|voi|img|nii)$','Select Atlas to convert');
end;
[niiname, isGz] = unGzSub (atlasfnm); %eg. aal.nii.gz -> aal.nii
[p,n,~] =spm_fileparts(niiname);
hdr = spm_vol(niiname);
if ~spm_type(hdr.dt,'intt'), error('Atlases must have integer datatype %s\n', niiname); end;
img = spm_read_vols(hdr);
for i = 1 : max(img(:)) %for each region
    imgI = img;
    imgI(imgI ~= i) = 0; %mask images not of index value
    imgI(imgI ~= 0) = 1; %binarize image
    numVx = sum(imgI(:));
    fprintf('region %d has %d voxels\n', i, numVx);
    if numVx < 1, continue; end; %skip empty regions
    hdr.fname = fullfile(p, [ n '_' num2str(i) '.nii']);  
    spm_write_vol(hdr,imgI);
end;
if isGz %delete uncompressed image: FSL does not allow 'img.nii' and 'img.nii.gz' to exist in same folder
    delete(niiname);
end
%end nii_atlas2nii()

function [fnm, isGz] = unGzSub (fnm)
isGz = false;
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz -> .nii
    isGz = true;
    fnm = char(gunzip(fnm));    
elseif strcmpi(ext,'.voi') %.voi -> .nii
    isGz = true;
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
end;  
%end unGzSub()