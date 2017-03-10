function nii_dtibatch (dtiBvecNames, isEddyCorrect)
%test all possible vectors and polarities
%assumes angulations have been correctly adjusted
% dtiNii: name of bvec file(s), e.g. img.bvec
% isEddyCorrect : if true simple undistortion applied, if false than quick and dirty
%Examples
% nii_dtibatch
% nii_dtibatch('DTI_axial_04.bvec');
% nii_dtibatch(strvcat('DTI_axial_04.bvec','DTI_30_deg_AP_05.bvec'));
%Chris Rorden 4 July 2014, BSD License
%Subsequently you can view the processed images from the command line
% fslview DTI_axial_04_FA.nii.gz DTI_axial_04_V1.nii.gz

if ~exist('dtiBvecNames','var') %file not specified
   [A,Apth] = uigetfile({'*.bvec';'*.*'},'Select b-vector file(s)', 'MultiSelect', 'on');
   dtiBvecNames = strcat(Apth,char(A));
end;
if ~exist('isEddyCorrect','var') %file not specified
   isEddyCorrect = true;
end;
fsldir= '/usr/local/fsl';
if ~exist(fsldir,'dir')
    error('Unable to find %s', fsldir); 
end
for i=1:size(dtiBvecNames,1)
    dtiBvec = deblank(dtiBvecNames(i,:)); %positive image
    [pth,nam,ext] = fileparts(dtiBvec);
    imgNam = fullfile(pth, [nam '.nii']); %img.nii
    if ~exist(imgNam,'file')
    	imgNam = fullfile(pth, [nam '.nii.gz']); %img.nii.gz
    end;
    if ~exist(imgNam,'file'), warning('Unable to find %s\n', imgNam); continue; end;
    refVol = refVolSub(dtiBvec);
    %next: permute all possible b-vector alternatives
    dtiSub(fsldir,imgNam,dtiBvec,refVol, isEddyCorrect);
    tracktionSub(dtiBvec);    
end
viewSub(fsldir,dtiBvec) %display results
%end main loop.... subroutines follow

function tracktionSub(dtiBvec)
if strcmpi(computer, 'GLNXA64')
    exeNam = 'tracktionLX'
else
   exeNam = 'tracktion'; 
end
exeNam = fullfile(fileparts(which(mfilename)), exeNam);
p = fileparts(exeNam);
if isempty(p)
   exeNam = fullfile(pwd, exeNam); 
end
if ~exist(exeNam,'file')
   fprintf('Skipped tractography: unable to find %s\n', exeNam);
   return; 
end
command=sprintf('%s "%s"\n',exeNam, dtiBvec);
system(command);
%end tracktionSub()

function ref = refVolSub(vNam) %find first B0 volume
[pth,nam,ext] = fileparts(vNam);
bNam = fullfile(pth, [ nam '.bval'] ); %Eddy corrected data
if ~exist(bNam, 'file'), error('Unable to find file %s',bNam); end;
b = textread(bNam);
ref = find(b==0, 1, 'first');
if isempty(ref), error('No b-zero volume in file %s', bNam); end;
ref = ref - 1;%fsl indexes volumes from 0
%end refVolSub()

function maskNam = betSub(fsldir,imgNam, refVol) %brain extract
setenv('FSLDIR', fsldir);
[pth,nam,ext] = fileparts(imgNam);
if strcmpi(ext,'.gz'), [~,nam] = fileparts(nam); end;
inNam = imgNam;
maskNam = fullfile(pth, nam ); %will generate image "dti_mask.nii.gz"
if (exist('refVol','var')) && (refVol > 0)
    %fslroi WIPDiffHR2.2SENSE4.2_ID123_12 msk 32 1
    refNam = fullfile(pth, [nam '_' num2str(refVol)] ); %will generate image "dti_mask.nii.gz"
    command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/fslroi %s %s %d 1"\n',imgNam,refNam, refVol);
    system(command);
    inNam = refNam;
end
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/bet %s %s -f 0.3 -g 0 -n -m"\n',inNam,maskNam);
maskNam = fullfile(pth, [nam '_mask.nii.gz']); %will generate image "dti_mask.nii.gz"
system(command);
%end betSub()

function eccNam = eddySub(fsldir,imgNam,refVol) %eddy current correct and brain extract
fprintf('Using eddy_correct: consider using eddy (and topup)\n');
[pth,nam,ext] = fileparts(imgNam);
eccNam = fullfile(pth, [nam, '_ecc']); %will generate image "dti_eddy.nii.gz"
setenv('FSLDIR', fsldir);
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin'])
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/eddy_correct %s %s %d"\n',imgNam,eccNam, refVol);
system(command);
eccNam = fullfile(pth, [nam, '_ecc.nii.gz']); %Eddy corrected data
%end eddySub()

function dtiSub(fsldir,imgNam,vNam, refVol, isEddyCorrect) %compute vectors
%%/usr/local/fsl/bin/dtifit --data=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti_eddy.nii.gz --out=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti --mask=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti_mask.nii.gz --bvecs=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/s004a001.bvec --bvals=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/s003a001.bval
[pth,nam,ext] = fileparts(vNam);
bNam = fullfile(pth, [ nam '.bval'] ); %Eddy corrected data
maskNam = betSub(fsldir,imgNam, refVol);
if ~exist(maskNam, 'file'), error('BET failed to create %s', maskNam); end;
if isEddyCorrect
    eccNam = eddySub(fsldir,imgNam,refVol);
else
    eccNam = imgNam;
end
[pth,nam,ext] = fileparts(vNam);
outNam = fullfile(pth, nam); %Eddy corrected data
setenv('FSLDIR', fsldir);
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin'])
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/dtifit --data=%s --out=%s --mask=%s --bvecs=%s --bvals=%s"\n',eccNam, outNam, maskNam,vNam,bNam);
system(command);
delete(fullfile(pth, [nam '_V2.nii.gz']));
delete(fullfile(pth, [nam '_V3.nii.gz']));
delete(fullfile(pth, [nam '_L1.nii.gz']));
delete(fullfile(pth, [nam '_L2.nii.gz']));
delete(fullfile(pth, [nam '_L3.nii.gz']));
delete(fullfile(pth, [nam '_MO.nii.gz']));
delete(fullfile(pth, [nam '_MD.nii.gz']));
delete(fullfile(pth, [nam '_S0.nii.gz']));
%end dtiSub

function viewSub(fsldir,vNam) %compute vectors
%%/usr/local/fsl/bin/dtifit --data=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti_eddy.nii.gz --out=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti --mask=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti_mask.nii.gz --bvecs=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/s004a001.bvec --bvals=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/s003a001.bval
[pth,nam,ext] = fileparts(vNam);
faNam = fullfile(pth, [nam '_FA.nii.gz']); %Eddy corrected data
v1Nam = fullfile(pth, [nam '_V1.nii.gz']); %Eddy corrected data
setenv('FSLDIR', fsldir);
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin'])
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/fslview %s %s &"\n',faNam,v1Nam);
system(command);
%end dtiSub