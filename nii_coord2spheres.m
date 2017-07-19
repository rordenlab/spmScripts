function nii_coord2spheres(coord, outname, exename)
%generate spheres at coordinates for display as Surf Ice nodes
% coord : Nx5 array of coordinates (X,Y,Z,Clr,Sz) see example
% outname: (optional) Name of output file, default is 'spheres.node'
% exename: (optional) Full path to Surf Ice
%Example
% nii_coord2spheres(); %demo
% c=[0,0,0,1,1; 40,0,0,2,2]; %1st at origin, 2nd 40mm right of origin
% nii_coord2spheres(c);

%check inputs - insert defaults if arguments not provided
if ~exist('coord','var')
   fprintf('Generating demo spheres\n'); 
   %3 spheres
   % first at origin
   % second 40mm right, twice as large
   % third 60mm posterior different color
   coord=[0,0,0,1,1; 40,0,0,1,2; 0,-60,0,2,1];
end
if ~exist('outname','var')
   outname='spheres.node';
end
if ~exist('exename','var') 
   exename = '/Users/rorden/Desktop/Surf_Ice/Surfice/surfice.app/Contents/MacOS/surfice';
end
if ~exist(exename,'file')
    fprintf('Hint: edit exename to automatically launch Surf Ice\n'); 
end
if size(coord,2) == 3
    fprintf('Assuming all spheres are of same size and color\n');
    coord = [coord, ones(size(coord,1),2)];
end
if size(coord,2) ~= 5 
   error('Each coordinate must have 5 properties [XYZCS] (XYZ=space, C=color, S=size)'); 
end
%now do all the work
fid = fopen(outname,'w');
for i = 1 : size(coord,1)
   fprintf(fid, '%g %g %g %g %g\n',coord(i,1), coord(i,2), coord(i,3), coord(i,4), coord(i,5)); 
end
fclose(fid);
fprintf('Use Surf Ice Nodes/AddNodes to open file %s\n',  outname);
%optional: show results
if ~exist(exename,'file'), return; end;
[p,n,x] = fileparts(outname);
if isempty(p), outname = fullfile(pwd,[n,x]); end;
cmd = sprintf('%s -S "begin meshload(''BrainMesh_ICBM152_smoothed.mz3''); nodeload(''%s''); SHADERXRAY(0.5, 0.5); end."', exename, outname);
system(cmd); 
