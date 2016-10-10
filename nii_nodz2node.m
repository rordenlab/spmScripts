function  nii_nodz2node (fnm)
%convert SurfIce .nodez file to BrainNet .node and .edge file
% fnm : name of nodz file to convert

if ~exist('fnm','var') || isempty(fnm)
    [nam,pth] = uigetfile({'*.nodz;';'*.*'},'Select Source image'); 
    if isequal(nam,0), return; end;
    fnm =fullfile (pth, nam);
end;
if ~exist(fnm,'file'), return; end;
str = fileread(fnm);
key = '#ENDNODE';
pos = strfind(str,key);
if isempty(pos), fprintf('Unable to find "%s" in %s\n', key, fnm); return; end;
%save node
[p,n] = fileparts(fnm);
fileID = fopen(fullfile(p, [n, '.node']),'w');
fprintf(fileID, str(1: (pos(1)-1)));
fclose(fileID);
%save edge
fileID = fopen(fullfile(p, [n, '.edge']),'w');
s = str((pos(1) + numel(key) + 1): end  ); %key plus EOLN character
s = strtrim(s); %just in case windows files have CR+LF EOLN 
fprintf(fileID, s);
fclose(fileID);
%nii_nodz2node
