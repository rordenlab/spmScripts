function nii_fslflirt (Src, Ref, Shadow,  OutDir, Prefix, IsInputMask,IsLinearInterp, InvRefWeight)
%Normalize P->Ref using FLIRT
% Src: image to normalize
% Ref: [optional] target of normalization
% Shadow: image[s] to be resliced using Src->Ref transform
% OutDir: [optional] output folder, defaults to Src's directory
% Prefix: [optional] prefix appended to warped images, defaults to 'w'
% IsInputMask : [optional] if 1 then 1st shadow used as 
% IsLinearInterp : [optional] if false then Shadow resliced using nearest neighbor
% InvRefWeight : [optional] Ref is masked with the inverse of this image
%Linear transform, mutual information cost function
%
%Note: requires FSL. One could use SPM's coregister function instead
%      I wrote this for data pre-processed with FSL
%Example: compute MNI->fMRI, then transform ROI from MNI to native space
%  nii_fslflirt('T1_LM1001.nii','DTIA_LM1001e.nii.gz','LS_LM1001.nii');
%  nii_fslflirt('T1_LM1001.nii','DTIA_LM1001e.nii.gz','LS_LM1001.nii','','',true);
%  nii_fslflirt('/Users/rorden/nii_rest/sct1_1mm','wT1_LM1001.nii','',pwd,'w',true,false,'LS_LM1001.nii');

if ~exist('Src','var') || isempty(Src)
    [nam,pth] = uigetfile({'*.nii;*.nii.gz;';'*.*'},'Select Source image'); 
    if isequal(nam,0), return; end;
    Src =fullfile (pth, nam);
end;
if ~exist('Ref','var') || isempty(Ref)
    [nam,pth] = uigetfile({'*.nii;*.nii.gz;';'*.*'},'Select Reference image'); 
    if isequal(nam,0), return; end;
    Ref =fullfile (pth, nam);
    %Ref = [fsldir '/data/standard/MNI152_T1_2mm_brain.nii.gz'] ;
end;
if ~exist('OutDir','var') || isempty(OutDir)
    OutDir = spm_fileparts(Src);
end;
Src = fullPathSub (Src);
Ref = fullPathSub (Ref);
if ~exist('Prefix','var') || isempty(Prefix), Prefix = 'w'; end;
if exist('IsInputMask','var') && IsInputMask && ~isempty(Shadow)
    ssrc = fullPathSub(deblank(Shadow(1,:)));
    [~,snam,sext] = spm_fileparts(ssrc);
    mask = fullfile(OutDir,['b' snam sext]);     
    fslCmdSub(sprintf('fslmaths %s -binv %s',ssrc,mask));
    mask = [' -inweight ' mask];
else
    mask = '';
end;
if exist('InvRefWeight','var') && ~isempty(InvRefWeight)
    InvRefWeight = fullPathSub (InvRefWeight);

    [~,snam,sext] = spm_fileparts(InvRefWeight);
    refmask = fullfile(OutDir,['b' snam sext]);     
    fslCmdSub(sprintf('fslmaths %s -binv %s',InvRefWeight,refmask));
    refmask = [' -refweight ' refmask];
else
    refmask = '';
end;

[~,nam,ext] = spm_fileparts(Src);
mat = fullfile(OutDir,[ nam '.mt']);
dest = fullfile(OutDir,[Prefix nam ext]);  
%command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/flirt -in %s -ref %s -out %s -omat %s -bins 256 -cost corratio -searchrx -45 45 -searchry -45 45 -searchrz -45 45 -dof 12  -interp trilinear"\n',Src,Ref,dest,mat);
%system(command);
fslCmdSub(sprintf('flirt -in %s -ref %s%s%s -out %s -omat %s -bins 256 -cost corratio -searchrx -45 45 -searchry -45 45 -searchrz -45 45 -dof 12  -interp trilinear',Src,Ref,mask,refmask,dest,mat));
if exist('Shadow','var') && ~isempty(Shadow) %Shadow Reg
    for i=1:size(Shadow,1)
        ssrc = fullPathSub(deblank(Shadow(i,:)));
        [~,snam,sext] = spm_fileparts(ssrc);
        dest = fullfile(OutDir,[Prefix snam sext]);  
        %command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/flirt -in %s -ref %s -out %s -applyxfm -init  %s   -interp trilinear"\n',ssrc,Ref,dest,mat);
        %system(command);
        if exist('IsLinearInterp','var') && ~IsLinearInterp 
            interp = 'nearestneighbour';
        else
            interp = 'trilinear';
        end
        fslCmdSub(sprintf('flirt -in %s -ref %s -out %s -applyxfm -init  %s   -interp %s',ssrc,Ref,dest,mat,interp));
    end;%for each shadow
end;%if shadows specified

%local functions follow

function fnm = fullPathSub (fnm)
%fsl requires a full path to files
if isempty(fnm), return; end;
[pth, nam, ext] = fileparts(fnm);
if ~isempty(pth), return; end;
fnm = fullfile(pwd, [nam ext]);
% fslNameSub()


function fslCmdSub (Cmd)
%execute a fsl command, e.g. fslCmd('fslinfo a.nii');
fsldir= '/usr/local/fsl/';
if ~exist(fsldir,'dir')
	error('%s: fsldir (%s) not found',mfilename, fsldir);
end
setenv('FSLDIR', fsldir);
flirt = [fsldir 'bin/flirt'];
if ~exist(flirt,'file')
	error('%s: flirt (%s) not found',mfilename,flirt);
end
command=sprintf('sh -c ". %setc/fslconf/fsl.sh; %sbin/%s"\n',fsldir,fsldir, Cmd);
fprintf(command);
system(command);
%end fslCmdSub()




