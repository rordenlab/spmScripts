function nii_diploe (T1, T2)
%Wide diploic space can disrupt segmentation-normalization
%This script darkens regions around brain based on T2 scan
% T1 : T1-weighted image with wide diploic space (marrow, cancellous bone)
% T2 : T2-weighted image
%Version
% Chris Rorden 20171109
%License
% This is a simple wrapper for SPM12, and carries the same license:
% GNU General Public Licence (either version 2, or at your option, any later version)
%Example
% nii_diploe('T1.nii', 'T2.nii');
% nii_diploe; %use GUI

if ~exist('T1','var')  %no T1 specified
 T1 = spm_select(1,'image','Select T1 scan');
end
if ~exist('T2','var') %no T2 specified
 T2 = spm_select(1,'image','Select T2 scan');
end
%init batch system https://en.wikibooks.org/wiki/SPM/Batch
spm('defaults','fmri');
spm_jobman('initcfg');
%coregister/reslice T2 to match T1
T2 = coregEstWriteSub(T1, T2, []); 
%segment T2 scan 
normSegSub(T2, [], false); %segment T2 scan
%brain extract T2 based on T2
T2bet = extractSub(0.25, T2, prefixSub('c1',T2), prefixSub('c2',T2), prefixSub('c3',T2));
%make T1 where we darken regions outside the brain
maskT1 = maskSub(T1, T2bet, 0.25, true);
%Segment the masked T1 scan
%normSegSub(maskT1, T2, true); %<- This uses both the T1 and T2 scans
normSegSub(maskT1, [], true); %<- This uses only the T1 scan
%end nii_diploe()

function [T2, lesion] = coregEstWriteSub(T1, T2, lesion)
%coregister T2 to match T1 image, apply to lesion
if isempty(T1) || isempty(T2), return; end;
fprintf('Coregistering %s to match %s\n',T2,T1);
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {T1};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {T2};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {lesion};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
T2 = prefixSub('r',T2);
if ~isempty(lesion), lesion = prefixSub('r',lesion); end;
%end coregEstSub()

function maskImgName = maskSub(imgName, maskName, modulation, rimMask)
%reduce voxels outside of mask by factor modulation
% imgName: image to mask
% mask: image with zeros in regions to be attenuated
% modulation : fraction of attenuation outside the mask
% rimMask : if true, only rim around bounary of mask will be attenuated
%Example
% maskSub('i.nii', 'msk.nii, 0.25, true); %decrease intensity outside mask to 25% input
%
hdr = spm_vol(maskName);
mask = spm_read_vols(hdr);
mask(mask ~= 0) = 1;
if rimMask
    inmask = mask + 0;
    spm_smooth(inmask, mask,4);
    mask(mask > 0.05) = modulation;
    mask(mask < modulation) = 1;
    mask(inmask==1) = 1;
end;
hdr = spm_vol(imgName);
img = spm_read_vols(hdr);
if ~isequal(size(mask), size(img)), error('Images size mismatch %s %s', imgName, maskName); end;
img = img .* mask;
maskImgName = prefixSub('m',imgName);
hdr.fname = maskImgName;
%hdr.dt(1)    =16;
spm_write_vol(hdr,img);
%end maskSub()

function t1Bet = extractSub(thresh, t1, c1, c2, c3)
%subroutine to extract brain from surrounding scalp
% t1: anatomical scan to be extracted
% c1: gray matter map
% c2: white matter map
% c3: [optional] spinal fluid map
fprintf('Brain extraction of %s\n', t1);
%[pth,nam,ext] = spm_fileparts(t1);
%load headers
mi = spm_vol(t1);%bias corrected T1
gi = spm_vol(c1);%Gray Matter map
wi = spm_vol(c2);%White Matter map
%load images
m = spm_read_vols(mi);
g = spm_read_vols(gi);
w = spm_read_vols(wi);
if nargin > 4 && ~isempty(c3)
   ci = spm_vol(c3);%CSF map
   c = spm_read_vols(ci);
   w = c+w;
end;
w = g+w;
if thresh <= 0
    m=m.*w;
else
    mask= zeros(size(m));
    mask(w >= thresh) = 255;
    spm_smooth(mask,mask,1); %feather the edges
    mask = mask / 255;
    m=m.*mask;
end;
mi.fname = prefixSub('b',t1);
%mi.fname = fullfile(pth,['b',  nam, ext]);
mi.dt(1) = 4; %16-bit precision more than sufficient uint8=2; int16=4; int32=8; float32=16; float64=64
spm_write_vol(mi,m);
t1Bet = mi.fname;
%end extractSub()

function normSegSub(img, img2, createDartel)
template = fullfile(spm('Dir'),'tpm','TPM.nii');
if ~exist(template,'file'), error('Unable to find template named %s',template); end
matlabbatch{1}.spm.spatial.preproc.channel(1).vols = {img};
matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];
if ~isempty(img2)
    matlabbatch{1}.spm.spatial.preproc.channel(2).vols = {img2};
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel(2).write = [0 0];
end
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[template ',1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 createDartel];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[template ',2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 createDartel];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[template ',3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 createDartel];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[template ',4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[template ',5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[template ',6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 2;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
spm_jobman('run',matlabbatch);
%end normSegSub()

function nam = prefixSub (pre, nam)
[p, n, x] = spm_fileparts(nam);
nam = fullfile(p, [pre, n, x]);
%end prefixSub()