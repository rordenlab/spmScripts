function [MeanSignal] = nii_timecourse_roi(IMGname, ROIname);
%Signal in region of interest: average signal in each volume of Img within mask image ROI
%  IMGname : NIfTI image to examine, can be 4D  
%  ROIname: NIFTI region[s] of interest where bright voxels define region
%Example
%  meansig = nii_timecourse_roi('rest_filt.nii','lfrontal.nii');
%  meansig = nii_timecourse_roi('rest_filt.nii',strvcat('rfrontal.nii','lfrontal.nii'));
%             to calculate correlation between these two regions:  R = corrcoef(meansig)

if nargin<1, IMGname = spm_select(1,'image','Select 4D image you wish to analyze'); end;
if nargin<2, ROIname = spm_select(inf,'image','Select region[s] of interest'); end;

if (exist(IMGname) == 0); MeanSignal = 0;fprintf('%s unable to find image %s\n',which(mfilename),IMGname);  return; end;
IMGname = nii_ungz(IMGname); %optional: decompress .voi/nii.gz to .nii files
ROIname = nii_ungz(ROIname); %optional: decompress .voi/nii.gz to .nii files

tic; %start timing

%load 4D dataset only once!
hdr = spm_vol(IMGname);
[XYZV] = spm_read_vols(hdr);
[nX nY nZ nV] = size(XYZV);


MeanSignal = zeros(nV, size(ROIname,1));
for r=1:size(ROIname,1) %for each region
    ROI = deblank(ROIname(r,:));
    if (exist(ROI) == 0); fprintf('%s unable to find image %s\n',which(mfilename),ROI);  return; end;
    rhdr = spm_vol(ROI);
    if (rhdr(1).dim(1) ~= hdr(1).dim(1)); fprintf('image and mask must have same number of columns/rows/slices: %s ~= %s\n',Img, ROI); return; end;
    mask = spm_read_vols(rhdr);
    mn = min(mask(:));
    mask = (mask ~= mn);%region of interest is all voxels that except the very darkest
    for vol = 1 : nV %for each volume
        imgmasked = XYZV(:,:,:,vol);
        imgmasked = imgmasked(mask);
        MeanSignal(vol,r) = mean(imgmasked);
    end; %for each volume
end; %for each region
toc; %report elapsed time

