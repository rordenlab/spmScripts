function coivox = nii_setOrigin12(vols, modality, cropBB)
%Align images so that origin and alignment roughly match MNI space
%  vols : cell string of image name(s) - first image used for estimate, others yoked
%  modality : modality of first image 1=T1, 2=T2, 3=EPI
%  cropBB : (optional) crop resulting image to standard bounding box
%Example
% nii_setOrigin12('T1_P001.nii', 1, true); %T1
% nii_setOrigin12({'T2_P171.nii','LS_P171.nii'}, 2, true); %T2 with yoked Lesion
% nii_setOrigin12({'APDTI_LM1021.nii.gz','PADTI_LM1021.nii.gz'}, 3, true); %DTI
%Chris Rorden 12/2014 (now supports SPM12)

% if ~exist('cropBB','var') 
%     cropBB = true;
%     fprintf('%s will crop image!\n', mfilename);
% end
if ~exist('vols','var') || isempty(vols) %no files specified
 vols = spm_select(inf,'image','Reset origin for selected image(s) (estimated from 1st)');
end
if ischar(vols), vols = cellstr(vols); end
vols = nii_ungz (vols, true, true);
if ~exist('modality','var') || isempty(modality) %no files specified
 modality = 1;
 fprintf('%s Modality not specified, assuming %d (1=T1,2=T2,3=EPI)\n', mfilename, modality);
end
nii_isSPM12orNewer;
setCenterOfIntensitySub(vols);
coregEstTemplateSub(vols, modality);
deleteMatFilesSub(vols);
if exist('cropBB','var') && (cropBB) %only if requested
 nii_clip2bb(vols, [], modality < 3, true, modality > 2); %clipZ for T1 and T2
end
%end MAIN FUNCTION - LOCAL FUNCTIONS FOLLOW

function deleteMatFilesSub(vols)
for v = 1:   numel(vols) 
    [pth, nam, ~, ~] = spm_fileparts(deblank(vols{v}));
    fname = fullfile(pth,[nam '.mat']);
    if exist(fname,'file')
        delete(fname);
    end
end %for each volume
%end deleteMatFilesSub()

function setCenterOfIntensitySub(vols)
coivox = ones(4,1);
%extract filename 
[pth,nam,ext, ~] = spm_fileparts(deblank(vols{1}));
fname = fullfile(pth,[nam ext]); %strip volume label
%report if filename does not exist...
if (exist(fname, 'file') ~= 2) 
 	fprintf('%s error: unable to find image %s.\n',mfilename,fname);
	return;  
end;
hdr = spm_vol([fname,',1']); %load header 
img = spm_read_vols(hdr); %load image data
img = img - min(img(:));
img(isnan(img)) = 0;
%find center of mass in each dimension (total mass divided by weighted location of mass
% img = [1 2 1; 3 4 3];
sumTotal = sum(img(:));
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZ_mm = hdr.mat * coivox; %convert from voxels to millimeters
fprintf('%s center of brightness differs from current origin by %.0fx%.0fx%.0fmm in X Y Z dimensions\n',fname,XYZ_mm(1),XYZ_mm(2),XYZ_mm(3)); 
for v = 1:   numel(vols) 
    fname = deblank(vols{v});
    if ~isempty(fname)
        [pth,nam,ext, ~] = spm_fileparts(fname);
        fname = fullfile(pth,[nam ext]); 
        hdr = spm_vol([fname ',1']); %load header of first volume 
        fname = fullfile(pth,[nam '.mat']);
        if exist(fname,'file')
            destname = fullfile(pth,[nam '_old.mat']);
            copyfile(fname,destname);
            fprintf('%s is renaming %s to %s\n',mfilename,fname,destname);
        end
        hdr.mat(1,4) =  hdr.mat(1,4) - XYZ_mm(1);
        hdr.mat(2,4) =  hdr.mat(2,4) - XYZ_mm(2);
        hdr.mat(3,4) =  hdr.mat(3,4) - XYZ_mm(3);
        spm_create_vol(hdr);
        if exist(fname,'file')
            delete(fname);
        end
    end
end%for each volume
%end setCenterOfIntensitySub()

function coregEstTemplateSub(vols, modality)
%vols: images to coregister - first used for estimate
if modality == 2
   template = fullfile(spm('Dir'),'canonical','avg152T2.nii');
elseif modality == 1
    template = fullfile(spm('Dir'),'canonical','avg152T1.nii');
else
    template  = fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii');
end
if ~exist(template,'file')
    error('Unable to find template named %s\n', template);
end
if ischar(vols)
    vols = cellstr(vols);
end
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {template};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[deblank(vols{1}),',1']};%{'/Users/rorden/Desktop/3D.nii,1'};
if  numel(vols) > 1
    matlabbatch{1}.spm.spatial.coreg.estimate.other = vols(2:end)';%transpose: column vector required!
else
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
end
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
%end coregEstTemplateSub()

