function nii_nii2stl (fnm, sthresh)
%input: T1 file in NIfTI format, output: mesh image
% fnm : input image
% sthresh : threshold for deciding if tissue is brain or air
%Similar
% http://bartferguson.nl/blog/?p=5
% http://hackaday.com/2015/08/25/you-own-your-mri-brainscan-do-something-interesting-with-it/
% http://www.pegors.net/3d-printing.html
% http://www.instructables.com/id/3D-Printing-from-MRI-data-in-5-steps/
% https://www.shapeways.com/blog/archives/21040-how-to-make-a-3d-print-of-your-brain.html
%Example
% nii_nii2stl('cr.nii')

if ~exist('fnm','var')  %fnmFA not specified
   fnm = spm_select(1,'image','Select T1 scan'); 
end;
if ~exist('sthresh','var'), sthresh = 0.15; end; 
if isempty(which('spm')), error('Install SPM12'); end;
if isempty(which('nii_reslice_target')), error('Put nii_reslice_target.m in your Matlab path'); end;
if true
    %1 get approximately aligned to template
    coregEstTemplateSub(fnm);
    %2 crop so we have a rough bounding box
    fnm = cropImgSub(fnm);
    % segment to find tissue types
    newSegSub(fnm);
    %3: set non-brain tissue to zero
else
    [p,n,x] = spm_fileparts(fnm);
    fnm = fullfile(p,['r',n,x]);
end
extractSub(sthresh, fnm, false); %save scalp extracted image
bfnm = extractSub(sthresh, fnm, true); %save binarized image
%4: make watertight
ffnm = watertightSub(bfnm);
%5 make mesh
makeMeshSub(ffnm);
%6: make hollow (cheaper)
hfnm = makeHollowSub(ffnm);
%7 make mesh
makeMeshSub(hfnm);
%end nii_nii2stl()

function fnm = cropImgSub(fnm)
%tarhdr.mat = [-1 0 0 79; 0 1 0 -113; 0 0 1 -51; 0 0 0 1];
%tarhdr.dim = [157 189 136];
tarhdr.mat = [-0.75 0 0 78.3750; 0 0.75 0 -113; 0 0 0.75 -51; 0 0 0 1];
tarhdr.dim = [209 252 181];
hdr = nii_reslice_target(fnm,'',tarhdr);
fnm = hdr.fname;
%cropImgSub()

function coregEstTemplateSub(vols)
%vols: images to coregister - first used for estimate
template = fullfile(spm('Dir'),'canonical','avg152T1.nii');
if ~exist(template,'file')
    error('Unable to find template named %s\n', template);
end
if ischar(vols)
    vols = cellstr(vols);
end
mbatch{1}.spm.spatial.coreg.estimate.ref = {template};
mbatch{1}.spm.spatial.coreg.estimate.source = {[deblank(vols{1}),',1']};%{'/Users/rorden/Desktop/3D.nii,1'};
if  numel(vols) > 1
    mbatch{1}.spm.spatial.coreg.estimate.other = vols(2:end);% {''};
else
    mbatch{1}.spm.spatial.coreg.estimate.other = {''};
end
mbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
mbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
mbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
mbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',mbatch);
%end coregEstTemplateSub()

function hfnm = makeHollowSub(fnm)
cropBottom = 15;
fprintf('Making hollow %s, cropping bottom %d slices for an escape hole\n', fnm, cropBottom);
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
img = (img - min(img(:))) / ( max(img(:)) - min(img(:)) ); %normalize 0..1
smth = img + 0.0;
spm_smooth(img,smth,4);
img(smth > 0.995) = 0;
%smooth a second time: hopefully marching-cubes will provide get subvoxel precision 
orig = img + 0.0;
spm_smooth(orig,img,1); %just one voxel fwhm - very mild
img(:,:,1:cropBottom) = 0;
[pth,nam,ext, ~] = spm_fileparts(deblank(fnm));
hfnm = fullfile(pth,['h',  nam, ext]);
hdr.fname = hfnm;
hdr.pinfo(1) = 1/255;
spm_write_vol(hdr,img);
%end makeHollowSub()

function sfnm = makeMeshSub(fnm)
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
thresh = max(img(:)) / 2;
reduce = 0.1; %reduce mesh complexity to 10%, smaller values for smaller files
FV = isosurface(img,thresh);
FV = reducepatch(FV,reduce);
%BELOW: FAST vector for converting from slice indices to mm
vx = [ FV.vertices ones(size(FV.vertices,1),1)];
vx = mtimes(hdr.mat,vx')';
FV.vertices = vx(:,1:3);
[pth,nam] = spm_fileparts(fnm);
sfnm = fullfile(pth,[nam, '.stl']);
stlwriteSub(sfnm, FV);
mfnm = fullfile(pth,[nam, '.mat']);
matwriteSub(mfnm, FV);
%end makeMeshSub()

function matwriteSub(fnm,  FV)
%Save mesh in MATLAB format
s.vertices = FV.vertices;
s.faces = FV.faces; %#ok<STRNU>
%s.vertexColors = [];
%s.colorMap = [];
%s.colorMin = [];
%explicitly request v7 for compression and compatibility
% http://undocumentedmatlab.com/blog/improving-save-performance
save(fnm,'-v7','-struct', 's');
fprintf('If you have installed MRIcroS, you can view this mesh using the command:\n');
fprintf('  MRIcroS(''addLayer'', ''%s'') \n', fnm);
%end writeMat()

function newSegSub(t1)
%apply new segment - return name of warping matrix
template = fullfile(spm('Dir'),'tpm','TPM.nii');
if ~exist(template,'file')
    error('Unable to find template named %s',template);
end
fprintf('NewSegment of %s\n', t1);
mbatch{1}.spm.spatial.preproc.channel(1).vols = {t1};
mbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.001;
mbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;
mbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];
mbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[template ',1']};
mbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
mbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
mbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
mbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[template ',2']};
mbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
mbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
mbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
mbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[template ',3']};
mbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
mbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
mbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
mbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[template ',4']};
mbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
mbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
mbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
mbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[template ',5']};
mbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
mbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
mbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
mbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[template ',6']};
mbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
mbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
mbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
mbatch{1}.spm.spatial.preproc.warp.mrf = 1;
mbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
mbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
mbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
mbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
mbatch{1}.spm.spatial.preproc.warp.samp = 3;
mbatch{1}.spm.spatial.preproc.warp.write = [0 0];
spm_jobman('run',mbatch);
%end newSegSub()

function t1Bet = extractSub(thresh, t1, probMap)   
%subroutine to extract brain from surrounding scalp
% t1: image extracted, requires corresponding 'c1','c2' tissue maps
fprintf('Brain extraction of %s\n', t1);
[pth,nam,ext] = spm_fileparts(t1);
c1 = fullfile(pth,['c1',nam,ext]);
c2 = fullfile(pth,['c2',nam,ext]);
if ~exist(c1,'file') || ~exist(c2,'file')
    error('Unable to find tissue maps %s %s', c1,c2);
end
%load headers
mi = spm_vol(t1);%bias corrected T1
gi = spm_vol(c1);%Gray Matter map
wi = spm_vol(c2);%White Matter map
%load images
m = spm_read_vols(mi);
g = spm_read_vols(gi);
w = spm_read_vols(wi);
w = g+w;
if thresh <= 0
    if ~exist('probMap','var') || ~probMap
        m=m.*w;
    end;
else
    mask= zeros(size(m));
    mask(w >= thresh) = 255;
    maskIn = mask;
    spm_smooth(maskIn,mask,1); %feather the edges
    mask = mask / 255;
    if ~exist('probMap','var') || ~probMap
        m=m.*mask;
    else
        m = mask;
    end;
    
end;
if exist('probMap','var') && probMap
   %  m = m * 1000;
    mi.pinfo(1) = 1/255;
    mi.fname = fullfile(pth,['p',  nam, ext]);
else
    mi.fname = fullfile(pth,['b',  nam, ext]);

end
mi.dt(1) = 4; %16-bit precision more than sufficient uint8=2; int16=4; int32=8; float32=16; float64=64
spm_write_vol(mi,m);
t1Bet = mi.fname;
%end extractSub()

function ffnm = watertightSub(fnm)
%this is going to be very slow....
fprintf('Flood fill of %s\n', fnm);
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
%img = dilateSub(img);
fimg = floodSub(img, true); %blurry image with center max values
img(fimg == max(fimg(:))) = max(img(:)); %fill ventricles
img = floodSub(img, false); %blurry image with center max values
[pth,nam,ext] = spm_fileparts(fnm);
ffnm = fullfile(pth,['f',  nam, ext]);
hdr.fname = ffnm;
hdr.pinfo(1) = 1/255;
spm_write_vol(hdr,img);
%end watertightSub()

% function img = dilateSub(img)
% img = (img - min(img(:))) / ( max(img(:)) - min(img(:)) );
% thresh = 0.2;
% smth = img + 0.0;
% spm_smooth(img,smth,2);
% img= zeros(size(img));
% img((smth >= thresh)) = 1; 
%%end dilateSub()

function img = floodSub(img, isFirstPass)
%similar to http://www.mathworks.com/matlabcentral/fileexchange/12184-floodfill3d
if isFirstPass
    thresh = max(img(:))/ 1000; %e.g. 0.1%
else
    thresh = max(img(:))/ 2; %e.g. 0.50%
end
imgBin = (img < thresh);
imgBin = double(imgBin);      % In case a logical matrix comes in.
imgBin(1,:,:) = NaN; %pad LEFT 
imgBin(end,:) = NaN; %pad RIGHT
imgBin(:,1,:) = NaN; %pad ANTERIOR 
imgBin(:,end,:) = NaN; %pad POSTERIOR
imgBin(:,:,1) = NaN; %pad INFERIOR 
imgBin(:,:,end) = NaN; %pad SUPERIOR
imgBin = floodFill3DSub(imgBin, [2,2,2]);
%fill bubbles
mx = max(img(:));
img(isfinite(imgBin))= mx;
%next: optional - blur to feather edges - only useful if marching cubes uses subvoxel smoothing
maskIn = img + 0.0;
if isFirstPass
    spm_smooth(maskIn,img,2); %feather the edges a lot: only center is filled
else
    spm_smooth(maskIn,img,1); %feather the edges a lot: weak blur
end
%end floodSub()

%Copyrights for Matlab Central Code
% stlwriteSub Copyright (c) 2015, Sven Holcombe
% floodFill3DSub Copyright (c) 2006,  F Dinath
% 
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function [A] = floodFill3DSub(A, point)
% http://www.mathworks.com/matlabcentral/fileexchange/12184-floodfill3d
% [B] = FloodFill3D(A, slice);
% This program flood fills a 6-connected 3D region. The input matrix MUST
% be a binary image. The user will select a seed (point) in the matrix to
% initiate the flood fill. You must specify the matrix slice in which you
% wish to place the seed.
% 
% A = binary matrix
% slice = a chosen slice in the matrix where you wish to place the seed.
%
% Enjoy,
% F. Dinath
if A(point(1), point(2), point(3));
    A(point(1), point(2), point(3)) = NaN;
    a{1} = sub2ind(size(A), point(1), point(2), point(3));
    i = 1;
    while 1
        i = i+1;
        a{i} = []; %#ok<AGROW>
        [x, y, z] = ind2sub(size(A), a{i-1});
        ob = nonzeros((A(sub2ind(size(A), x, y, z-1)) == 1).*sub2ind(size(A), x, y, z-1));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        ob = nonzeros((A(sub2ind(size(A), x, y, z+1)) == 1).*sub2ind(size(A), x, y, z+1));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        ob = nonzeros((A(sub2ind(size(A), x-1, y, z)) == 1).*sub2ind(size(A), x-1, y, z));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        ob = nonzeros((A(sub2ind(size(A), x+1, y, z)) == 1).*sub2ind(size(A), x+1, y, z));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        ob = nonzeros((A(sub2ind(size(A), x, y-1, z)) == 1).*sub2ind(size(A), x, y-1, z));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        ob = nonzeros((A(sub2ind(size(A), x, y+1, z)) == 1).*sub2ind(size(A), x, y+1, z));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        if isempty(a{i});
            break;
        end
    end
end
%end floodFill3DSub()

function stlwriteSub(fnm, FV)
% http://www.mathworks.com/matlabcentral/fileexchange/36770-stlwrite-write-binary-or-ascii-stl-file
%   Original idea adapted from surf2stl by Bill McDonald. Huge speed
%   improvements implemented by Oliver Woodford. Non-Delaunay triangulation
%   of quadrilateral surface courtesy of Kevin Moerman. FaceColor
%   implementation by Grant Lohsen.
%
%   Author: Sven Holcombe, 11-24-11
% Create the facets
facets = single(FV.vertices');
%The following line is similar to MeshLab's Filters/Normals/InvertFacesOrientation

%facets([1,3],:)=facets([3,1],:);
facets = reshape(facets(:,FV.faces'), 3, 3, []);
% Compute their normals
V1 = squeeze(facets(:,2,:) - facets(:,1,:));
V2 = squeeze(facets(:,3,:) - facets(:,1,:));
normals = V1([2 3 1],:) .* V2([3 1 2],:) - V2([2 3 1],:) .* V1([3 1 2],:);
clear V1 V2
normals = bsxfun(@times, normals, 1 ./ sqrt(sum(normals .* normals, 1)));
facets = cat(2, reshape(normals, 3, 1, []), facets);
clear normals
fid = fopen(fnm, 'w');
fprintf(fid, '%-80s', 'stlwrite from Sven Holcombe');             % Title
fwrite(fid, size(facets, 3), 'uint32');           % Number of facets
% Add one uint16(0) to the end of each facet using a typecasting trick
facets = reshape(typecast(facets(:), 'uint16'), 12*2, []);
% Set the last bit to 0 (default) or supplied RGB
facets(end+1,:) = 0;
fwrite(fid, facets, 'uint16');
%end stlwriteSub()
% Close the file
fclose(fid);
fprintf('Wrote %d facets\n',size(facets, 2));
%end stlwriteSub()
