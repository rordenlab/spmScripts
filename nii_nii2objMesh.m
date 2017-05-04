function outnm = nii_nii2objMesh (fnm, thresh, clusterVox, isSmooth, reduce, outnm, floodfill)
%convert NIfTI image to mesh, all arguments are optional
% fnm : nifti image to meshify
% thresh : air/surface threshold, e.g. 2=voxels darker than 2 are air
%          n.b. for Altases, set thresh=0: one image will be created per region
% clusterVox : only clusters of with more voxels than this value survive
% isSmooth : if true, input is blurred, otherwise raw data
% reduce : reduction ratio, e.g. 0.2 will decimate 80% of triangles
% outnm : name for output image, e.g. "mesh.obj"
% floodfill : if true, remove bubbles inside objects (e.g. ventricles embedded below cortex)
%Examples
% nii_nii2objMesh %use GUI
% nii_nii2objMesh('motor.nii', 3, 30, false, 0.3);
% nii_nii2objMesh('MNI152_T1_1mm_brain.nii.gz', 5500, 1000000, true, 0.1);
% nii_nii2objMesh('HarvardOxford-cort-maxprob-thr0-1mm.nii.gz', 0, 1, true, 0.1); %convert smoothed atlas
% nii_nii2objMesh('natbrainlab.nii.gz', 0, 1, false, 0.1); %convert unsmoothed atlas
if ~exist('fnm','var')  %no files specified
    [files,pth] = uigetfile({'*.gz;*.nii;*.hdr;';'*.*'},'Choose images to merge', 'MultiSelect', 'off');
    fnm = strcat(pth,char(files));
end
%load image
[hdr,img, fnm] = loadNiiSub(fnm);
img(isnan(img)) = 0;
%set preferences
if ~exist('outnm','var') 
    [pth nm] = spm_fileparts(fnm);
    outnm = fullfile(pth, [nm '.obj']);
end;
if ~exist('thresh','var') thresh = (0.5 * max(img(:))-min(img(:))) + min(img(:)); end;
if ~exist('clusterVox','var') clusterVox = 32; end;
if ~exist('isSmooth','var') isSmooth = false; end;
if  ~exist('floodfill','var') floodfill = false; end;
fprintf('Image intensity range %g..%g\n', min(img(:)), max(img(:)));
if ~exist('reduce','var') %provide user interface if prefences not specified
    prompt = {'Enter threshold (2=only voxels brighter than 2, 0=image is atlas)','Enter minimum cluster size (voxels)','Smooth (blur) input(1=yes, 0=no)','Reduction ratio (0.5=remove 50% of vertices)'};
    dlg_title = 'Values for thresholding';
    num_lines = 1;
    def = {num2str(thresh),num2str(clusterVox), num2str(isSmooth), '0.3'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    thresh = str2double(answer{1});
    clusterVox = str2double(answer{2});
    isSmooth = str2double(answer{3});
    reduce = str2double(answer{4});
end;
if thresh == 0
    if clusterVox > 1, fprintf('nb: cluster size ignored for Atlases'); end;
    if floodfill, fprintf('nb: floodfill ignored for Atlases'); end;
    saveAtlasSub(hdr, img, isSmooth, reduce, outnm);
    return;
end
if isSmooth
    presmooth = img+0; %+0 forces new matrix
    spm_smooth(presmooth,img,1,0); %smooth data
end
isInvert = thresh < 0;
if isInvert
   img = -img;
   thresh = -thresh;
end
[img, nVoxPost] = clusterSub(hdr, img, thresh, clusterVox, isInvert);
if nVoxPost < 1, fprintf('No voxels exceed threshold\n'); return; end;
if floodfill
    img = floodSub(img, thresh, 2); %optional: fill bubbles inside
end;
img2meshSub(hdr, img, outnm, reduce, thresh);
%end nii_nii2obj()

function saveAtlasSub(hdr, img, isSmooth, reduce, outnm)
if ~(spm_type(hdr.dt,'intt')) %integer output
    error('Threshold=0 suggests Atlas, but image is not saved as discrete integers.');
end;
mn = min(img(:))+1;
mx = max(img(:));
if (mn >= mx), error('No variability in input image'); end;
if isSmooth, fprintf('n.b. each region will be smoothed'); end;
[pth nm] = spm_fileparts(outnm);
for i = mn : mx
    outnm = fullfile(pth, [nm '_' num2str(i) '.obj']);
    bw = (img == i)+0; %black and white image, +0 converts logical->double
    if isSmooth
        presmooth = bw+0; %+0 forces new matrix
        spm_smooth(presmooth,bw,1,0); %smooth data
    end
    vx = sum(bw(:) >= 0.5);
    if vx > 0
        fprintf('Region %d has %d voxels\n', i, vx);
        img2meshSub(hdr, bw, outnm, reduce, 0.5);
    end;
end
%end saveAtlasSub()

function img2meshSub(hdr, img, outnm, reduce, thresh)
if (max(img(:))< thresh) || (min(img(:)) > thresh) return; end;
FV = isosurface(img,thresh);
FV = reducepatch(FV,reduce);
FV.vertices = FV.vertices(:,[2,1,3]); %isosurface swaps X/Y
vx = [ FV.vertices ones(size(FV.vertices,1),1)];
vx = mtimes(hdr.mat,vx')';
FV.vertices = vx(:,1:3);
FV.faces = fliplr(FV.faces); %reverse winding
writeObjSub(FV.vertices,FV.faces, outnm);
%end img2meshSub()

function [img, nVoxPost] = clusterSub(hdr, img, thresh, clusterVox, isInvert)
bw=img;
bw = (bw >= abs(thresh)) * 1.0;
nVox = sum(bw(:));
nVoxPost = nVox;
if (clusterVox < 2) || (nVox < 1), return; end;
[bw,nCluster] = spm_bwlabel(bw,18);
nClusterPost = nCluster;
minClusterVox = nCluster;
maxClusterVox = 0;
for i = 1:nCluster
    vox = sum(bw(:)== i);
    %fprintf('Cluster %d has %d voxels\n', i,vox);
    minClusterVox = min(minClusterVox,vox);
    maxClusterVox = max(maxClusterVox,vox);
    if vox < clusterVox
        nVoxPost = nVoxPost - vox;
        nClusterPost = nClusterPost - 1;
        bw(bw == i) = 0;        
    end
end
%any voxels that do not survive threshold need to be set just below threshold
bw = (bw < 1); %mask now inverted and logical: TRUEs are voxels that can be reduced
maxsubthresh = thresh - (2 * eps);%max(img(img(:) < thresh));
mask =  (img >= thresh) .* bw; %mask voxels that are over thresh undersized
img(mask == 1) = maxsubthresh;
if isInvert, thresh = - thresh; end;
if minClusterVox ~= maxClusterVox
    fprintf('Cluster sizes range from %d to %d voxels\n',minClusterVox, maxClusterVox);
end
outmm3=prod(abs(hdr.mat(1:3, 1:3)*[1;1;1]));
fprintf('%d voxels (%d clusters) exceed %g, of which %d (%d) are in clusters larger than %d voxels (%gmm^3)\n', ...
    nVox, nCluster, thresh, nVoxPost, nClusterPost, clusterVox, clusterVox * outmm3);
%end clusterSub()

function [hdr,img, fnm] = loadNiiSub(fnm)
isGz = false;
if ~exist('spm_fileparts','file'), fprintf('Please install SPM12 or later'); end;
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));
    isGz = true;
elseif strcmpi(ext,'.voi') %.voi ->
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
    isGz = true;
end;
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
if isGz, delete(fnm); end;
%end loadNiiSub()
function writeObjSub(vertex,face,filename)
% --- save Face/Vertex data as WaveFront Object format file
%inputs:
%	vertex: vertices matrix where cols are xyz and each row a vertix
%	face: face matrix where cols are xyz and each row is face
%	fileName: the Wavefront Object file to create
%notes
% https://en.wikipedia.org/wiki/Wavefront_.obj_file
[nF nFd] =size(face);
[nV nVd] =size(vertex);
if (nF <1) || (nV <3 || (nFd ~=3) || (nVd ~=3)), warning('Problem with writeObj'); return; end;  %#ok<WNTAG>
fid = fopen(filename, 'wt');
fprintf(fid, '# WaveFront Object format image created with MRIcroS\n');
fprintf(fid, 'v %.12g %.12g %.12g\n', vertex');
fprintf(fid, 'f %d %d %d\n', (face)');
fclose(fid);
%end writeObjSub()

function img = floodSub(img, thresh, smoothVox)
%similar to http://www.mathworks.com/matlabcentral/fileexchange/12184-floodfill3d
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
if exist('smoothVox', 'var')
    spm_smooth(maskIn, img, smoothVox); %feather the edges a lot: weak blur
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