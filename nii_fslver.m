function [major, minor, rev] = nii_fslver(fsldir, majorReq, minorReq, revReq )
%Returns fsl version
% fsldir : (optional) location of fsl
% majorReq : (optional) minimum required major version
% minorReq : (optional) minimum required major version
% minorRevision : (optional) minimum required revision
%Examples
% [major, minor, rev] = nii_fslver; %read version
% nii_fslver('',5,0,10); %error if version prior to 5.0.10
if ~exist('fsldir','var') || isempty(fsldir)
    fsldir= '/usr/local/fsl';
end
fnm = fullfile(fsldir,'etc','fslversion');
major = 0;
minor = 0;
rev = 0;
if ~exist(fnm,'file')
    warning('Unable to determine FSL version: no file named %s',fnm);
    return; 
end;
[major, minor, rev] = textread(fnm, '%d %d %d', 1, 'delimiter', '.');
if ~exist('majorReq','var') || (major > majorReq), return; end;
if ~exist('minorReq','var') || (minor > minorReq), return; end;
if ~exist('revReq','var') || (rev >= revReq), return; end;
error('Require FSL version %d.%d.%d but found %d.%d.%d', majorReq, minorReq, revReq, major, minor, rev);
%end nii_fslver



