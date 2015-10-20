function plotBvec(fnm)
%Show sampling of DTI bvec file
% fnm : name of FSL format BVec file
%Notes
% http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/EDDY
% http://www.mccauslandcenter.sc.edu/CRNL/tools/advanced-dti
% http://www.emmanuelcaruyer.com/q-space-sampling.php
%Examples
% plotBvec;
% plotBvec('myVec.bvec');

if nargin <1 %no files
 fnm = spm_select(1,'^.*\.(bvec)$','Select bvec file to view');
end;
bvecs = load(fnm); % Assuming your filename is bvecs
figure('position',[100 100 500 500]);
plot3(bvecs(1,:),bvecs(2,:),bvecs(3,:),'*r');
axis([-1 1 -1 1 -1 1]);
axis vis3d;
rotate3d
%end plotBvec()