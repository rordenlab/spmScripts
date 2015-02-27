function nii_render (P, PreserveCSF, Normalize, Thresh, DilateVox)
%use segmentation-normalization to generate scalp stripped image for rendering
% P : input images
% Thresh: Threshold, if 0.01 then tissue >1% gray/white will be preserved
%     Note: Thresh=0 modulates output (T1*(wm+gm))
% PreserveCSF : [optional] if true, subsequent rendering includes spinal fluid
%     Useful for showing brain injury.
%     If true, Threshold based on wm+gm+csf: consider ~0.5
% Normalize: [optional] if true, tissue warped to standard space
% DilateVox: [optional] Pad edges of the brain with a few voxels
%
%examples
% nii_render('T1.nii');
% nii_render('T1.nii',false,true,0.05);

if nargin <1 %no files
 P = spm_select(inf,'image','Select images to scalp strip');
end;
if nargin <2 %Keeping spinal fluid not specified
 PreserveCSF = false; %default: thorough cleanup
end;
if nargin <3 %Normalize not specified
    Normalize = false;
end;
if nargin <4 %Threshold not specified
    if PreserveCSF
        Thresh = 0.75; %at least 75% gm+wm+csf 
    else
        Thresh = 0.005; %a very low threshold effectively dilates brain, so feathering CSF not brain
    end;
end;
if nargin < 5 %Normalize not specified
    DilateVox = 0;
end;
cleanup = 2; %2= thorough cleanup; 1=light cleanup, 0= nocleanup
gm = fullfile(spm('Dir'),'tpm','grey.nii');
wm = fullfile(spm('Dir'),'tpm','white.nii');
csf = fullfile(spm('Dir'),'tpm','csf.nii');
spm('defaults','fmri');
spm_jobman('initcfg');
for i=1:size(P,1)
    ref = deblank(P(i,:));
    [pth,nam,ext] = spm_fileparts(ref);
    matlabbatch{1}.spm.spatial.preproc.data = {ref};    
    if Normalize
        matlabbatch{1}.spm.spatial.preproc.output.GM = [0 1 0];
        matlabbatch{1}.spm.spatial.preproc.output.WM = [0 1 0];
        if PreserveCSF 
            matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 1 0];
        else
            matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];            
        end;
        matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
    else
        matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
        matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 1];
        if PreserveCSF 
            matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 1];
        else
            matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];            
        end;
        matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
    end;
    matlabbatch{1}.spm.spatial.preproc.output.cleanup = cleanup;
    matlabbatch{1}.spm.spatial.preproc.opts.tpm = {gm; wm; csf};
    matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2; 2; 2; 4];
    matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
    matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
    matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
    matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};
    fprintf('Unified segmentation of %s with cleanup level %d threshold %f, job %d/%d\n', ref, cleanup, Thresh, i, size(P,1));
    fprintf('  If segmentation fails: use SPM''s DISPLAY tool to set the origin as the anterior commissure\n');    
    spm_jobman('run',matlabbatch);
    if Normalize
        writenormsub(fullfile(pth,['m',  nam, ext]), fullfile(pth,[  nam,'_seg_sn.mat'])) ;
        if PreserveCSF
            extractsub(Thresh,fullfile(pth,['wm',  nam, ext]),fullfile(pth,['wc1',  nam, ext]),fullfile(pth,['wc2',  nam, ext]),fullfile(pth,['wc3',  nam, ext]), DilateVox);
        else
            extractsub(Thresh,fullfile(pth,['wm',  nam, ext]),fullfile(pth,['wc1',  nam, ext]),fullfile(pth,['wc2',  nam, ext]),'', DilateVox);            
        end;
    else %if normalized else native space
        if PreserveCSF
            extractsub(Thresh,fullfile(pth,['m',  nam, ext]),fullfile(pth,['c1',  nam, ext]),fullfile(pth,['c2',  nam, ext]),fullfile(pth,['c3',  nam, ext]), DilateVox);
        else
            extractsub(Thresh,fullfile(pth,['m',  nam, ext]),fullfile(pth,['c1',  nam, ext]),fullfile(pth,['c2',  nam, ext]),'', DilateVox);            
        end;    
    end;
end; %for each image...
%end main function

function writenormsub (img, mat)
%subroutine to reslice image (img) using spatial transforms (mat);
matlabbatch{1}.spm.spatial.normalise.write.subj.matname = {mat};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {img};
matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [NaN NaN NaN];
matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
%end writenormsub()

function extractsub(thresh, t1, c1, c2, c3, DilateVox)   
%subroutine to extract brain from surrounding scalp
% t1: anatomical scan to be extracted
% c1: gray matter map
% c2: white matter map
% c3: [optional] spinal fluid map
[pth,nam,ext] = spm_fileparts(t1);
%load headers
mi = spm_vol(t1);%bias corrected T1
gi = spm_vol(c1);%Gray Matter map
wi = spm_vol(c2);%White Matter map
%load images
m = spm_read_vols(mi);
g = spm_read_vols(gi);
w = spm_read_vols(wi);
if ~isempty(c3)
   ci = spm_vol(c3);%CSF map
   c = spm_read_vols(ci);
   w = c+w; 
end;
w = g+w;
mi.fname = fullfile(pth,['render',  nam, ext]);
if thresh <= 0
    m=m.*w;
else
    mask= zeros(size(m));
    for px=1:length(w(:)),
      if w(px) >= thresh
        mask(px) = 255;
      end;
    end;
    if DilateVox > 0
        mask = dilatesub(mask, DilateVox);
    end;
    spm_smooth(mask,mask,1); %feather with 1mm FWHM
    mask = mask/255;
    m=m.*mask;
end;
mi.dt(1)=4;%save as 16-bit uint8=2; int16=4; int32=8; float32=16; float64=64
spm_write_vol(mi,m); 
%end extractsub()

function [result] = dilatesub(mask, ndilate)
     kx=[0.25 0.5 0.25];
     ky=[0.25 0.5 0.25];
     kz=[0.25 0.5 0.25];
     mn = min(mask(:));
     for j=1:ndilate,
         spm_conv_vol(mask,mask,kx,ky,kz,-[1 1 1]);
     end;
     result= zeros(size(mask));
     result((mask > mn)) = 255;
%end dilatesub()

