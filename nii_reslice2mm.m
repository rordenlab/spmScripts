function fnm = nii_reslice2mm(fnm)
%reslice 3D or 4D data to be 2mm isotropic with a standard bounding box
if ~exist('fnm','var') || isempty(fnm)
    fnm = spm_select(1,'image','Select image to reslice');
end;
[p,n,x] = spm_fileparts(fnm); %strip volume, e.g. 'img.nii,1'->img.nii
fnm = fullfile(p, [ n, x]);
inhdr = spm_vol(fnm);
inimg = spm_read_vols(inhdr); %load input image
outhdr = inhdr(1);
outhdr.mat = [-2 0 0 80; 0 2 0 -114; 0 0 2 -52; 0 0 0 1]; %<- 2mm
outhdr.dim = [79 95 69]; %<-2mm,
%1mm-> outhdr.mat = [-1 0 0 79; 0 1 0 -113; 0 0 1 -71; 0 0 0 1];
%1mm-> outhdr.dim = [157 189 156];
outimg = zeros(outhdr.dim);
fnm = fullfile(p, ['x', n, x]);
outhdr.fname = fnm;
for vol = 1: numel(inhdr) %for every volume of 4D data
    for i = 1: outhdr.dim(3)
        M = inv(spm_matrix([0 0 -i])*inv(outhdr.mat)*inhdr(1).mat);
        outimg(:,:,i) = spm_slice_vol(inimg(:,:,:,vol), M, outhdr.dim(1:2), 1); %",0"=nearest,", 1"=linear interp)
    end
    outhdr.n(1)=vol;
    spm_write_vol(outhdr,outimg); %save image to disk
end;
%end reslice2mm()