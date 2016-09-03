function nii_makeDTI
%Generates a simple NIfTI format DTI image
% demonstrates voxel offset
fnm = 'test';
pixDimMM = 3; %distance between voxel centers
%bvec = makeBVec7; %6 diffusion bvectors small file but somewhat uneven
bvec = makeBVec33; %32 diffusion bvectors large file 
writeBvecBval (bvec, fnm);
nVol = size(bvec,1);
dim = [9, 9, 9, nVol]; %image resolution in columns, rows, slices, volumes
dtype = 32; %precision of data
ofile    = [fnm, '.nii'];%spm_file(parfile,'path',opts.outdir,'ext',opts.ext);
scale = 1;
inter = 0;
switch dtype
    case 8
        dtype = spm_type('int8'); % uint8?
    case 16
        dtype = spm_type('int16'); % uint16?
    case 32
        dtype = spm_type('float32');
    case 64
        dtype = spm_type('float64');
    otherwise
        error('Unknown data type.');
end
mat = eye(4);
mat(1:3, 1:3) = mat(1:3, 1:3) * pixDimMM;
mat(1,4) = -pixDimMM * ((dim(1)+1)/2); %center, indexed from 1 not 0
mat(2,4) = -pixDimMM * ((dim(2)+1)/2);
mat(3,4) = -pixDimMM * ((dim(3)+1)/2);
%SPM's matrices are indexed from 1
% hdr = spm_vol('test.nii,1');
% mm = (hdr.mat * [0 0 1 1; 0 0 9 1]') %mm of 1st and 9th slice (z=-12..12)
%if you use fslhd you will see matrix indexed from 0
% >fslhd test.nii
% mat0 = [3 0 0 -12; 0 3 0 -12; 0 0 3 -12; 0 0 0 1]
% mm = (mat0 * [0 0 0 1; 0 0 8 1]') %mm of 1st and 9th slice (z=-12..12)
dato     = file_array(ofile,dim,[dtype 0],0,scale,inter);
N        = nifti;
N.dat    = dato;
N.mat    = mat;
N.mat0   = mat;
N.mat_intent  = 'Scanner';
N.mat0_intent = 'Scanner';
N.descrip     = 'made with Matlab';
%N.timing = struct('toffset',[],'tspace',[]); % store TR
create(N);
for i=1:nVol
    N.dat(:,:,:,i) = makeDtiSub(dim(1:3), bvec(i,:));
end
%end nii_mk()

function writeBvecBval (bvec, fnm)
dlmwrite([fnm '.bvec'],bvec', 'delimiter','\t');
bval = zeros(size(bvec,1),1);
for i = 1 : size(bvec,1) %B=0 images have zero-length vectors, otherwise B-weighted
    if norm(bvec(i,:)) >= 0.01
        bval(i) = 1000;
    end
end
dlmwrite([fnm '.bval'],bval', 'delimiter','\t');
dlmwrite([fnm '.trakvis.bvec'],bvec(bval > 0,:), 'delimiter','\t');
%end writeBvecBval()

function bvec = makeBVec33
%Returns a set of bvectors: 0,0,0 for B=0
% http://www.emmanuelcaruyer.com/q-space-sampling.php
bvec=[0 0 0;
    0.049	-0.919	-0.391
0.726	0.301	-0.618
-0.683	0.255	-0.684
0.845	-0.502	-0.186
-0.73	-0.619	-0.288
-0.051	0.039	0.998
-0.018	0.871	-0.491
-0.444	0.494	0.747
-0.989	-0.086	-0.116
-0.47	-0.855	0.221
0.412	0.4	0.819
-0.552	0.79	-0.267
-0.123	-0.477	0.871
-0.848	0.141	0.51
-0.341	-0.788	-0.512
0.361	-0.529	0.768
-0.472	0.85	0.234
-0.856	-0.481	0.189
0.797	0.162	0.582
0.467	-0.009	-0.884
0.013	0.998	-0.056
0.882	-0.387	0.267
0.017	-0.536	-0.844
-0.442	-0.651	0.617
0.365	-0.058	0.929
0.977	-0.004	-0.213
-0.406	-0.902	-0.145
-0.627	0.614	0.479
-0.354	0.772	-0.528
-0.658	-0.472	-0.586
0.423	0.322	-0.847
0.212	-0.754	-0.622];
%end makeBVec33()

function bvec = makeBVec7
%Returns a set of bvectors: 0,0,0 for B=0
% http://www.emmanuelcaruyer.com/q-space-sampling.php
bvec=[0 0 0;
    0.049	-0.919	-0.391;
    0.726	0.301	-0.618;
    -0.683	0.255	-0.684;
    0.845	-0.502	-0.186;
    -0.730	-0.619	-0.288;
    -0.051	0.039	0.998;
    -0.018	0.871	-0.491];
%end makeBVec7()

function img = makeDtiSub (dim, bvec)
isBZero = norm(bvec) < 0.01; %if B=0 then unweighted image
img = zeros(dim);
center = (dim+1)/2; %e.g. center of 9 voxels is the 5th
sz = norm(center);
for z = 1 : dim(3)
	vec = z/dim(3);
    if vec < 0.5
        vec = [0 0 1]; %head/foot
    else
        vec = [1 0 0]; %left/right
    end
    dwi = norm(cross(bvec,vec)); %is this orthogonal
    if isBZero
       dwi = 1; 
    end
    %fprintf('%g\n', dwi);
	for y = 1 : dim(2) 
		for x = 1 : dim(1)
            dx = (sz -norm(center-[x,y,z]))/sz;
            if (~isBZero) && (dx < 0.5)
                dx = 0; 
            end
			img(x,y,z)=dwi * dx;
		end;
	end;
end;
img = img * 4096; %TrackVis likes large integer values
%end makeDtiSub()