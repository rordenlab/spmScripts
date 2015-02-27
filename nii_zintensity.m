function nii_zintensity (Subj, Maskname);
% Creates version of image where intensity is z-score of brightness in mask
%  Filename: Continuous image[s] used for median calculation
%  Mask: (optional). list of mask image[s] - region used for peak
%
%Example
%  mdn = nii_median('brain.nii','mask.nii');
%  mdn = nii_median('MNI152_T1_2mm_brain.nii.gz	','mask.voi');

%select files with a dialog
if nargin<1
    Subj = spm_select(inf,'image','Select image(s) for intensity normalization (Z scores)');
end;

if (nargin < 2) %user did not specify a mask - request one...
   [Maskname, pathname] = uigetfile({'*.nii;*.img;*.gz;*.voi';'*.*'},'Select the Mask'); 
   Maskname = fullfile(pathname, Maskname);
end;
%load image
for i=1:size(Subj,1)
    Filename = deblank(Subj(i,:));
    Filename = nii_ungz(Filename); %optional: convert .nii.gz to .nii
    vi = spm_vol(Filename);
    img = spm_read_vols(vi);
    %load mask
    Maskname = nii_ungz(Maskname);
    vm = spm_vol(Maskname);
    mask = spm_read_vols(vm);
    % make mask binary...
    mn = min(mask(:));
    mask = (mask ~= mn);
    imgmasked = img(mask);
    %return result
    format long;
    mn=mean(imgmasked);
    st= std(imgmasked);
    fprintf('%s has %d voxels, after masking with %s the remaining %d voxels have a mean intensity of %f\n',Filename,length(img(:)), Maskname,length(imgmasked(:)),mn); 
    if (st == 0) 
        fprintf('nii_zintensity error: can not compute Z-scores when standard error is zero!');
        return;
    end;
    % next part saves transformed data
    %z=(img-mn)./ st;
     z((img==0)) = 0; %<- use this to retain zeros....
    [pth,nam,ext]=fileparts(vi.fname);
    vi.dt(1)=16; %save as 32-bit float
    vi.fname = fullfile(pth,['z',  nam, ext]);
    spm_write_vol(vi,z);
end;