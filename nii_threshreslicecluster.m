function nii_threshreslicecluster(srcNam, tarNam, thresh, clusterMM3, binarize)
%Reslice image to isotropic 1mm resolution, zero dark voxels, zero small clusters
%  srcNam: (optional) name of NIfTI image to reslice
%  tarNam: (optional) image to match- either filename of NIfTI header or loaded NIfTI header structure
%  thresh: (optional) threshold (e.g. 3 then voxels >3 survive, if -3 then voxels<-3 survive)
%  clusterMM3: (optional) only clusters larger than this value survive (cubic millimeter)
%  binarize: (optional) if true: image is binarized as 0 or 1, if false unmasked voxel intensities not changed
% Outputs: thresholded mage
%Examples
% nii_threshreslicecluster; %use gui
% nii_threshreslicecluster('unthresh.nii','ch256.nii.gz',3,1000)
%For usage example see
% https://f1000research.com/articles/4-466/v1
isGz = 1; %1 to convert .nii -> .nii.gz, else output uncompressed
isSaveTemp = true; %save resliced image before filter
if ~exist('clusterMM3','var'), clusterMM3 = 32; end;
if ~exist('binarize','var'), binarize = false; end;
if ~exist('srcNam','var')
   [nam,pth] = uigetfile({'*.nii;*.hdr;*.nii.gz';'*.*'},'Select image to reslice and threshold');
   if isequal(nam,0), return; end;
   srcNam =fullfile (pth, nam);
end
if ~exist('tarNam','var')
   [nam,pth] = uigetfile({'*.nii;*.hdr;*.nii.gz';'*.*'},'Select target template image');
   if isequal(nam,0), return; end;
   tarNam =fullfile (pth, nam);
end
if ~exist('thresh','var')
    prompt = {'Enter threshold:','Enter minimum cluster size (mm^3)','Compress output (1=yes,0=no)', 'Binarize (1=yes,0=no)'};
    dlg_title = 'Values for thresholding';
    num_lines = 1;
    def = {'3','864','1','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    thresh = str2double(answer{1});
    clusterMM3 = str2double(answer{2});
    isGz = str2double(answer{3});
    binarize = str2double(answer{4});
end;
srcNam = unGzSub (srcNam);
inhdr = spm_vol(srcNam); %load input header
inimg = spm_read_vols(inhdr); %load input image
tarNam = unGzSub (tarNam);
tarhdr = spm_vol(tarNam); %load input header
outhdr            = inhdr;
[pth,nam,ext] = fileparts(inhdr.fname);
outhdr.fname      = fullfile(pth,['r' nam ext]);
outhdr.dim = tarhdr.dim;
outhdr.mat = tarhdr.mat;
outimg = zeros(outhdr.dim(1:3));
for i = 1:outhdr.dim(3)
    M = inv(spm_matrix([0 0 -i])*inv(outhdr.mat)*inhdr.mat);
    outimg(:,:,i) = spm_slice_vol(inimg, M, outhdr.dim(1:2), 2); % (1=linear interp, 2=spline)
end
%threshold
if thresh > 0
    outimg(outimg(:) < thresh) = 0;
else
    outimg(outimg(:) > thresh) = 0;
end
if isSaveTemp
    outimg(isnan(outimg)) = 0;
    outhdr.fname      = fullfile(pth,['temp' nam ext]);
    outhdr = spm_create_vol(outhdr); %save header to disk
    spm_write_vol(outhdr,outimg); %save image to disk
    outhdr.fname      = fullfile(pth,['r' nam ext]);
end
%use >0 not ~=0 so NaNs are not detected
fprintf('%s has %d voxels that exceed %g \n', outhdr.fname,sum(outimg(:)> 0), thresh);
%cluster threshold
if clusterMM3 > 0
    outmm3=prod(abs(outhdr.mat(1:3, 1:3)*[1;1;1]));
    bw=outimg;
    bw(~isfinite(bw)) = 0;
    bw(bw ~=0) = 1;
    [bw,nCluster] = spm_bwlabel(bw,18);
    nLargeCluster = 0;
    clustVox = ceil(clusterMM3/outmm3);
    for i = 1:nCluster
        %fprintf('Cluster %d has %d voxels\n', i,sum(bw(:)== i));
        if sum(bw(:)== i) < clustVox
            bw(bw == i) = 0;
        else
            nLargeCluster = nLargeCluster + 1;
            fprintf('Cluster %d has %d voxels\n', i,sum(bw(:)== i));
        end
    end
    outimg(bw(:) == 0) = 0;
    fprintf('%s has %d voxels (each %gmm^3) that exceed %g and are part of the %d clusters of at least %dmm^3 (%d voxels))\n', ...
        outhdr.fname,sum(outimg(:)~= 0), outmm3, thresh, nLargeCluster,clusterMM3, clustVox);
end %if clusterVox

%binarize
if binarize, 
    outimg(~isfinite(outimg)) = 0;
    outimg(outimg ~=0) = 1; 
end;
if sum(outimg(:)~= 0) < 1
    fprintf('No voxels survive thresholding\n');
    return
end
%write filtered image to disk
outhdr = spm_create_vol(outhdr); %save header to disk
spm_write_vol(outhdr,outimg); %save image to disk
if (isGz) && (strcmpi(ext,'.nii')) % compress .nii files, but not .hdr files
    gzip(outhdr.fname);
    delete(outhdr.fname);
end
%end nii_threshreslice()


function fnm = unGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));
elseif strcmpi(ext,'.voi') %.voi ->
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
end;
%end unGzSub()