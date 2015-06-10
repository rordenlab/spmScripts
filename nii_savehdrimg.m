function nii_savehdrimg(fname, hdr, img)
if ~isfield(hdr.private, 'hk')
    error('Incompatible header format (must be read by nii_loadhdrimg)');
end
%next: adjust private header to etch visible protions
hdr.private.dime.datatype = hdr.dt(1); %datatype
hdr.private.dime.scl_slope = hdr.pinfo(1); %slope
hdr.private.dime.scl_inter = hdr.pinfo(2); %intercept
nDim = numel(size(img));
if nDim < 3, nDim = 3; end;
hdr.private.dime.dim = ones(1,8);
hdr.private.dime.dim(1) = nDim;
for i = 1: nDim
    if size(img,i) > 0
        hdr.private.dime.dim(i+1) = size(img,i);
    end
end

if (hdr.dt(1) == 128) 
    hdr.private.dime.dim(1) = nDim; 
    hdr.private.dime.dim(4) = (numel(img)/size(img,1))/size(img,2)/3;
    hdr.private.dime.dim(5) = 1;
end
%now save the data
write_niiSub(fname, hdr, img);
%end nii_savehdrimg()

function write_niiSub(fname, hdrx, img)
% Subfunction: simplified from Jimmy Shen's NIfTI toolbox
hdr = hdrx.private;
[pth, nam, ext] = fileparts(fname);
isGz = false;
if strcmpi(ext,'.gz')
    isGz = true;
    [~, nam, ext] = fileparts(nam);
    if ~strcmpi(ext,'.nii')
        [pth, nam] = fileparts(fname);
        ext = '.nii';
    end
end
fileprefix = fullfile(pth, nam); %without extension
switch double(hdr.dime.datatype),
    case   1,    hdr.dime.bitpix = int16(1 );    precision = 'ubit1';
    case   2,    hdr.dime.bitpix = int16(8 );    precision = 'uint8';
    case   4,    hdr.dime.bitpix = int16(16);    precision = 'int16';
    case   8,    hdr.dime.bitpix = int16(32);    precision = 'int32';
    case  16,    hdr.dime.bitpix = int16(32);    precision = 'float32';
    case  32,    hdr.dime.bitpix = int16(64);    precision = 'float32';
    case  64,    hdr.dime.bitpix = int16(64);    precision = 'float64';
    case 128,    hdr.dime.bitpix = int16(24);    precision = 'uint8';
    case 256,    hdr.dime.bitpix = int16(8 );    precision = 'int8';
    case 511,    hdr.dime.bitpix = int16(96);    precision = 'float32';
    case 512,    hdr.dime.bitpix = int16(16);    precision = 'uint16';
    case 768,    hdr.dime.bitpix = int16(32);    precision = 'uint32';
    case 1024,   hdr.dime.bitpix = int16(64);    precision = 'int64';
    case 1280,   hdr.dime.bitpix = int16(64);    precision = 'uint64';
    case 1792,   hdr.dime.bitpix = int16(128);   precision = 'float64';
    otherwise
        error('This datatype is not supported');
end
hdr.dime.glmax = round(double(max(img(:))));
hdr.dime.glmin = round(double(min(img(:))));
if strcmpi(ext,'.nii')
    fid = fopen([fileprefix '.nii'], 'w');
    if fid < 0, error('Cannot open file %s.nii.',fileprefix); end
    hdr.dime.vox_offset = 352;        
    hdr.hist.magic = 'n+1';
    save_nii_hdrSub(hdr, fid);
    skip_bytes = double(hdr.dime.vox_offset) - 348;
    fwrite(fid, zeros(1, skip_bytes), 'uint8');
else
    fid = fopen([fileprefix '.hdr'], 'w');
    if fid < 0, error('Cannot open file %s.hdr.', fileprefix); end
    hdr.dime.vox_offset = 0;
    hdr.hist.magic = 'ni1';
    save_nii_hdrSub(hdr, fid);
    fclose(fid);
    fid = fopen([fileprefix '.img'], 'w');
end
fwrite(fid, img, precision);
fclose(fid);
if isGz
    gzip([fileprefix '.nii']);
    delete([fileprefix '.nii']);
end
%end write_nii()

% Subfunction: simplified from Jimmy Shen's NIfTI toolbox
function save_nii_hdrSub(hdr, fid)
fseek(fid, 0, -1);
fwrite(fid, hdr.hk.sizeof_hdr(1),    'int32');	% must be 348.
fwrite(fid, padChar(hdr.hk.data_type, 10), 'uchar');
fwrite(fid, padChar(hdr.hk.db_name, 18), 'uchar');
fwrite(fid, hdr.hk.extents(1),       'int32');
fwrite(fid, hdr.hk.session_error(1), 'int16');
fwrite(fid, hdr.hk.regular(1),       'uchar');	% might be uint8
fwrite(fid, hdr.hk.dim_info(1),      'uchar');
fwrite(fid, hdr.dime.dim(1:8),        'int16');
fwrite(fid, hdr.dime.intent_p1(1),  'float32');
fwrite(fid, hdr.dime.intent_p2(1),  'float32');
fwrite(fid, hdr.dime.intent_p3(1),  'float32');
fwrite(fid, hdr.dime.intent_code(1),  'int16');
fwrite(fid, hdr.dime.datatype(1),     'int16');
fwrite(fid, hdr.dime.bitpix(1),       'int16');
fwrite(fid, hdr.dime.slice_start(1),  'int16');
fwrite(fid, hdr.dime.pixdim(1:8),   'float32');
fwrite(fid, hdr.dime.vox_offset(1), 'float32');
fwrite(fid, hdr.dime.scl_slope(1),  'float32');
fwrite(fid, hdr.dime.scl_inter(1),  'float32');
fwrite(fid, hdr.dime.slice_end(1),    'int16');
fwrite(fid, hdr.dime.slice_code(1),   'uchar');
fwrite(fid, hdr.dime.xyzt_units(1),   'uchar');
fwrite(fid, hdr.dime.cal_max(1),    'float32');
fwrite(fid, hdr.dime.cal_min(1),    'float32');
fwrite(fid, hdr.dime.slice_duration(1), 'float32');
fwrite(fid, hdr.dime.toffset(1),    'float32');
fwrite(fid, hdr.dime.glmax(1),        'int32');
fwrite(fid, hdr.dime.glmin(1),        'int32');
fwrite(fid, padChar(hdr.hist.descrip, 80), 'uchar');
fwrite(fid, padChar(hdr.hist.aux_file, 24), 'uchar');
fwrite(fid, hdr.hist.qform_code,    'int16');
fwrite(fid, hdr.hist.sform_code,    'int16');
fwrite(fid, hdr.hist.quatern_b,   'float32');
fwrite(fid, hdr.hist.quatern_c,   'float32');
fwrite(fid, hdr.hist.quatern_d,   'float32');
fwrite(fid, hdr.hist.qoffset_x,   'float32');
fwrite(fid, hdr.hist.qoffset_y,   'float32');
fwrite(fid, hdr.hist.qoffset_z,   'float32');
fwrite(fid, hdr.hist.srow_x(1:4), 'float32');
fwrite(fid, hdr.hist.srow_y(1:4), 'float32');
fwrite(fid, hdr.hist.srow_z(1:4), 'float32');
fwrite(fid, padChar(hdr.hist.intent_name, 16), 'uchar');
fwrite(fid, padChar(hdr.hist.magic, 4), 'uchar');
if ~isequal(ftell(fid), 348), error('Header size is not 348 bytes.'); end
%end save_nii_hdrSub()

% Subfunction: pad or chop char to correct length. Called by save_nii_hdr
function buf = padChar(ch, len)
len1 = length(ch);
if len1 >= len,  
    buf = ch(1:len);
else
    buf = [ch zeros(1, len-len1, 'uint8')];
end
%end padChar()
