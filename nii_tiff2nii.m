function nii_tiff2nii (files)
%Convert TIFF stack to NIFTI
% files : (optional) cell string with filenames to convert
%Examples
% nii_tiff2nii %use GUI
% nii_tiff2nii({'DWI_STAC.tiff'})
% nii_tiff2nii({'DWI_STAC.tiff','test.TIF'})

if ~exist('files','var')
    [files,pth] = uigetfile({'*.tiff;';'*.*'},'Select TIFF Image[s]', 'MultiSelect', 'on'); 
    files = strcat(pth, files); %add path to all file names
end;
files = cellstr(files); %make cellstr nii_tiff2nii({'DWI_STAC.tiff'}) == nii_tiff2nii('DWI_STAC.tiff')
niiname = fullfile(spm('Dir'),'canonical','avg152T1.nii'); %sample nifti file for header
if ~exist(niiname,'file')
    error('Unable to find a NIfTI image named %s', niiname);
end;
hdr = spm_vol([niiname,',1']); 

for f=1:size(files,2)
    %determine image information
    nam = strvcat(deblank(files(:,f)));
    i = imfinfo(nam);
    if i(1).BitsPerSample ~= 8
       warning('Only tested on 8-bit data'); 
    end
    sz = [i(1).Width,  i(1).Height, numel(i)];
    fprintf('dimensions\t%s\t%d\t%d\t%d\n',nam, sz(1), sz(2), sz(3));
    %load all 2D slices as a single 3D volume
    img = zeros(sz(1), sz(2), sz(3));
    for z = 1 : sz(3)
        s = imread(nam, z)';
        if (size(s,1) ~= sz(1)) || (size(s,2) ~= sz(2)) 
            error('Slice %d has different dimensions %d x %d', z, size(s,1), size(s,2) );
        end
        img(:,:,z) = s;
    end
    %create NIfTI header
    [p,n] = spm_fileparts(which(nam));
    [~, pdir2] = fileparts(p); %get file's directory name
    hdr.fname = fullfile(p,[ pdir2, '_', n, '.nii']);
    hdr.mat = eye(4);
    %set voxel size in mm
    hdr.mat(1,1) = -256/sz(1);
    hdr.mat(2,2) = -256/sz(2);
    hdr.mat(3,3) = 130/sz(3);
    %set origin
    hdr.mat(1,4) = -0.5 * hdr.mat(1,1) * sz(1);
    hdr.mat(2,4) = -0.5 * hdr.mat(2,2) * sz(2);
    hdr.mat(3,4) = -0.5 * hdr.mat(3,3) * sz(3);
    hdr.dim = sz;
    hdr.pinfo = [0;0;0];
    spm_write_vol(hdr,img);
    nii_setOrigin12(hdr.fname,3,true);
end;
%end nii_tiff2nii()
