function [isFullSphere, meanLength] = nii_meanBvec(vNam) %does bvec sample whole sphere?
%given FSL style bvec, determine if sampling over full or half sphere
%useful for determining if one should run eddy (requires full sphere) or eddy_correct
%see nii_plotBvec for visualization
% vNam : name of bvec file
%Output
% isFullSphere: returns true if mean vector near zero (vectors cancel each other)
% meanLength: length of mean vector, near zero for full sphere
%Examples
% isFull = nii_meanBvec('a.bvec');
% isFull = nii_meanBvec(); %use GUI

if ~exist('vNam','var') %no files
 vNam = spm_select(1,'^.*\.(bvec)$','Select bvec file to view');
end;
if ~exist(vNam, 'file'), error('Unable to find file %s',vNam); end;
v = textread(vNam);
v( :, all( ~any( v ), 1 ) ) = []; %delete vectors with all zeros (e.g. B=0)
v = mean(v,2); %mean vector
meanLength = sqrt(sum((v) .^ 2)); %mean vector length
%if sperical then vectors cancel out and mean length near zero
%if half-sphere than the "center of mass" for vectors is biased
isFullSphere = (meanLength < 0.25);
%end isFullSphereSub()