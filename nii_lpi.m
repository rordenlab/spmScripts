function nii_lpi(fnms)
%Flip image to be in LPI orientation
% fnm: image to flip
%Rationale
% http://www.diedrichsenlab.org/imaging/suit_function.htm
%  The algorithm works best if the T1 image is brought into LPI-orientation, and the origin of the image is set to the anterior commissure.
%Examples
% nii_lpi()
% nii_lpi('T1.nii')
% nii_lpi({'s1T1.nii', 's22T1.nii'})

if ~exist('fnms', 'var') %user did not specify a mask - request one...
   [fnm, pathname] = uigetfile({'*.nii;*.img;*.gz;*.voi';'*.*'},'Select the Mask', 'MultiSelect', 'on'); 
   fnms = fullfile(pathname, fnm);
end
fnms = cellstr(fnms); 
if isempty(which('nii_tool'))
    error('Make sure nii_tool is in your Matlab path: https://github.com/xiangruili/dicm2nii');
end
for f=1:numel(fnms)
    nii2lpi(fnms{f});
end
%end nii_lpi()

function nii2lpi(fnm)
nii = nii_tool('load', fnm);
R0 = nii_xform_mat(nii.hdr); % original background R
[~, perm, flp] = reorient(R0, nii.hdr.dim(2:4), 0);
f = [0, 0, 0];
p = [1, 2, 3];
if isequal(p,perm) && isequal(f,flp) 
    fprintf('Already in LPI %s\n', fnm);
    return;
end
[nii, perm, flp] = nii_reorient(nii, false);
[p,n,x] = fileparts(fnm);
onm = fullfile(p,['lpi_',n,x]);
nii_tool('save', nii, onm);
fprintf('Converted to LPI %s\n', onm);
%end nii2lpi()

% normalize columns
function v = normc(M)
v = bsxfun(@rdivide, M, sqrt(sum(M .* M)));
% v = M ./ sqrt(sum(M .* M)); % since 2016b

% reorient nii to diagnal major
function [nii, perm, flp] = nii_reorient(nii, leftHand, ask_code)
if nargin<3, ask_code = []; end
[R, frm] = nii_xform_mat(nii.hdr, ask_code);
dim = nii.hdr.dim(2:4);
pixdim = nii.hdr.pixdim(2:4);
[R, perm, flp] = reorient(R, dim, leftHand);
fps = bitand(nii.hdr.dim_info, [3 12 48]) ./ [1 4 16];
if ~isequal(perm, 1:3)
    nii.hdr.dim(2:4) = dim(perm);
    nii.hdr.pixdim(2:4) = pixdim(perm);
    nii.hdr.dim_info = [1 4 16] * fps(perm)' + bitand(nii.hdr.dim_info, 192);
    nii.img = permute(nii.img, [perm 4:8]);
end
sc = nii.hdr.slice_code;
if sc>0 && flp(fps==3)
    nii.hdr.slice_code = sc+mod(sc,2)*2-1; % 1<->2, 3<->4, 5<->6
end
if isequal(perm, 1:3) && ~any(flp), return; end
if frm(1) == nii.hdr.sform_code % only update matching form
    nii.hdr.srow_x = R(1,:);
    nii.hdr.srow_y = R(2,:);
    nii.hdr.srow_z = R(3,:);
end
if frm(1) == nii.hdr.qform_code
    nii.hdr.qoffset_x = R(1,4);
    nii.hdr.qoffset_y = R(2,4);
    nii.hdr.qoffset_z = R(3,4);
    R0 = normc(R(1:3, 1:3));
    dcm2quat = dicm2nii('', 'dcm2quat', 'func_handle');
    [q, nii.hdr.pixdim(1)] = dcm2quat(R0);
    nii.hdr.quatern_b = q(2);
    nii.hdr.quatern_c = q(3);
    nii.hdr.quatern_d = q(4);
end
for i = find(flp), nii.img = flip(nii.img, i); end

% Return true if input is char or single string (R2016b+)
function tf = ischar(A)
tf = builtin('ischar', A);
if tf, return; end
if exist('strings', 'builtin'), tf = isstring(A) && numel(A)==1; end

% Return xform mat and form_code: form_code may have two if not to ask_code
function [R, frm] = nii_xform_mat(hdr, ask_code)
% [R, form] = nii_xform_mat(hdr, asked_code);
% Return the transformation matrix from a NIfTI hdr. By default, this returns
% the sform if available. If the optional second input, required form code, is
% provided, this will try to return matrix for that form code. The second
% optional output is the form code of the actually returned matrix.
fs = [hdr.sform_code hdr.qform_code]; % sform preferred
if fs(1)==fs(2), fs = fs(1); end % sform if both are the same
f = fs(fs>=1 & fs<=4); % 1/2/3/4 only
if isempty(f) || ~strncmp(hdr.magic, 'n', 1) % treat it as Analyze
    frm = 0;
    try % try spm style Analyze
        [pth, nam, ext] = fileparts(hdr.file_name);
        if strcmpi(ext, '.gz'), [~, nam] = fileparts(nam); end
        R = load(fullfile(pth, [nam '.mat']));
        R = R.M;
    catch % make up R for Analyze: suppose xyz order with left storage 
        R = [diag(hdr.pixdim(2:4)) -(hdr.dim(2:4).*hdr.pixdim(2:4)/2)'; 0 0 0 1];
        R(1,:) = -R(1,:); % use left handed
    end
    return;
end
if numel(f)==1 || nargin<2 || isempty(ask_code) % only 1 avail or no ask_code
    frm = f;
else % numel(f) is 2, numel(ask_code) can be 1 or 2
    frm = f(f == ask_code(1));
    if isempty(frm) && numel(ask_code)>1, frm = f(f == ask_code(2)); end
    if isempty(frm) && any(f==2), frm = 2; end % use confusing code 2
    if isempty(frm), frm = f(1); end % no match to ask_code, use sform
end

if frm(1) == fs(1) % match sform_code or no match
    R = [hdr.srow_x; hdr.srow_y; hdr.srow_z; 0 0 0 1];
else % match qform_code
    R = quat2R(hdr);
end

% Reorient 4x4 R
function [R, perm, flp] = reorient(R, dim, leftHand)
% [R, perm, flip] = reorient(R, dim, leftHand)
% Re-orient transformation matrix R (4x4), so it will be diagonal major and
% positive at diagonal, unless the optional third input is true, which requires
% left-handed matrix, where R(1,1) will be negative. 
% The second input is the img space dimension (1x3). 
% The perm output, like [1 2 3] or a permutation of it, indicates if input R was
% permuted for 3 axis. The third output, flip (1x3 logical), indicates an axis 
% (AFTER perm) is flipped if true.
a = abs(R(1:3,1:3));
[~, ixyz] = max(a);
if ixyz(2) == ixyz(1), a(ixyz(2),2) = 0; [~, ixyz(2)] = max(a(:,2)); end
if any(ixyz(3) == ixyz(1:2)), ixyz(3) = setdiff(1:3, ixyz(1:2)); end
[~, perm] = sort(ixyz);
R(:,1:3) = R(:,perm);
flp = R([1 6 11]) < 0; % diag(R(1:3, 1:3))
if nargin>2 && leftHand, flp(1) = ~flp(1); end
rotM = diag([1-flp*2 1]);
rotM(1:3, 4) = (dim(perm)-1) .* flp; % 0 or dim-1
R = R / rotM; % xform matrix after flip