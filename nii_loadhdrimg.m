function [hdr, img] = nii_loadhdrimg(filename)
%load NIfTI (.nii, .nii.gz, .hdr/.img) image and header 
% filename: image to open
%To do:
%  endian: rare, currently detected and reported but not handled
%Examples
% hdr = nii_loadhdrimg('myimg.nii');
% [hdr, img] = nii_loadhdrimg('myimg.nii');
%Similar to following SPM functions (though handles gz)
% hdr = spm_vol(filename);
% img = spm_read_vols(hdr);
    
if ~exist('filename','var')  %fnmFA not specified
   [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select image');
   filename = [Apth, A];
end;
[fpth, fnam,fext] = fileparts(filename);
if strcmpi(fext,'.img') %hdr/img pair
    filename = fullfile(fpth, [fnam, '.hdr']);
end
if ~exist(filename, 'file')
    error('Unable to find file %s', filename);
end
%load data
if strcmpi(fext,'.gz') %unzip compressed data
	%http://undocumentedmatlab.com/blog/savezip-utility
    %http://www.mathworks.com/matlabcentral/fileexchange/39526-byte-encoding-utilities/content/encoder/gzipdecode.m
    streamCopier = com.mathworks.mlwidgets.io.InterruptibleStreamCopier.getInterruptibleStreamCopier;
    baos = java.io.ByteArrayOutputStream;
    fis  = java.io.FileInputStream(filename);
    zis  = java.util.zip.GZIPInputStream(fis);
    streamCopier.copyStream(zis,baos);
    fis.close;
    data = baos.toByteArray;
else
    fileID = fopen(filename);
    data = fread(fileID);
    data = uint8(data);
    fclose(fileID);
end
%read header
hdr = spm_vol_Sub(filename, data);
if nargout < 2, return; end; %only read image if requested
if strcmpi(fext,'.hdr') || strcmpi(fext,'.img') %analyze style .hdr and .img pairs
    if ~exist(hdr.fname, 'file')
        error('Unable to find image %s', hdr.fname);
    end
    fileID = fopen(hdr.fname);
    data = fread(fileID);
    data = uint8(data);
    fclose(fileID);
end
img = spm_read_vols_Sub(hdr, data);
%end nii_loadhdrimg()

function img = spm_read_vols_Sub(hdr, data)
% --- load NIfTI voxel data: mimics spm_read_vol without requiring SPM
nSamplesPerVoxel = 1;
switch hdr.dt(1)
   case   2,
      bitpix = 8;  myprecision = 'uint8';
   case   4,
      bitpix = 16; myprecision = 'int16';
   case   8,
      bitpix = 32; myprecision = 'int32';
   case  16,
      bitpix = 32; myprecision = 'single';%'float32';
   case  64,
      bitpix = 64; myprecision = 'double';%'float64';
   case   128,
      bitpix = 8;  myprecision = 'uint8';
      nSamplesPerVoxel = 3;
   case 512 
      bitpix = 16; myprecision = 'uint16';
   case 768 
      bitpix = 32; myprecision = 'uint32';
   otherwise
      error('This datatype is not supported %d', hdr.dt(1)); 
end
if numel(hdr.dim) > 3
    nVol = prod(hdr.dim(4:end));
else
    nVol = 1; %3D data has only a single volume
end
myvox = hdr.dim(1)*hdr.dim(2)*hdr.dim(3)*nVol*nSamplesPerVoxel;
%ensure file is large enough
imgbytes = myvox * (bitpix/8); %image bytes plus offset
offset = double(hdr.pinfo(3));
if (imgbytes+offset) > numel(data)
    fprintf('Error: expected %d but file has %d bytes %s',imgbytes, file_stats.bytes,hdr.fname);
    return;
end;
%read data

img = typecast(data(offset+1:offset+imgbytes),myprecision);%fread(fid, myvox, myprecision, 0, myformat);
img = double(img);
if nSamplesPerVoxel > 1
    %
else
    img = img(:).*hdr.pinfo(1)+hdr.pinfo(2); %apply scale slope and intercept
    img = reshape(img, hdr.dim(1), hdr.dim(2), hdr.dim(3), nVol);
end;
%end spm_read_vols_Sub()

function [Hdr] = spm_vol_Sub(filename, data)
[h, machine] = readHdrSub (data);
nDim = find(h.dime.dim > 1,1,'last') -1; %-1 since dim[2]=x, dim[3]=y, etc
if nDim < 3, nDim = 3; end;
Hdr.dim = ones(1,nDim);
for i = 1: nDim
    if (h.dime.dim(i+1) > 0), Hdr.dim(i) = h.dime.dim(i+1); end;
end
%Hdr.dim
%Hdr.dim = double([h.dime.dim(2) h.dime.dim(3) h.dime.dim(4)]);
%Hdr.dim
if (h.hist.sform_code == 0) && (h.hist.qform_code == 0)
    fprintf('Warning: no spatial transform detected. Perhaps Analyze rather than NIfTI format');
    Hdr.mat = fileUtils.nifti.hdr.hdr2m(h.dime.dim,h.dime.pixdim );
elseif (h.hist.sform_code == 0) && (h.hist.qform_code > 0) %use qform Quaternion only if no sform
    Hdr.mat = fileUtils.nifti.hdr.quarternion.hdrQ2m(h.hist,h.dime.dim,h.dime.pixdim );
else %precedence: get spatial transform from matrix (sform)
    Hdr.mat = [h.hist.srow_x; h.hist.srow_y; h.hist.srow_z; 0 0 0 1];
    Hdr.mat = Hdr.mat*[eye(4,3) [-1 -1 -1 1]']; % mimics SPM: Matlab arrays indexed from 1 not 0 so translate one voxel
end;
if strcmpi(machine, 'ieee-le')
	Hdr.dt = [h.dime.datatype 0];
else
	Hdr.dt = [h.dime.datatype 1];
end;
Hdr.pinfo = [h.dime.scl_slope; h.dime.scl_inter; h.dime.vox_offset];
if fileUtils.isExt('.hdr',filename)
	[pth, nam] = fileparts(filename);
    Hdr.fname =  fullfile(pth, [nam '.img']); %if file.hdr then set to file.img
else
	Hdr.fname =  filename;
end
Hdr.descrip = h.hist.descrip;
Hdr.n = [h.dime.dim(5) 1];
Hdr.private.hk = h.hk;
Hdr.private.dime = h.dime;
Hdr.private.hist = h.hist;
%end spm_vol_Sub()

function [h, machine] = readHdrSub (data)
machine = 'ieee-le';
%read header key
hk.sizeof_hdr = typecast(data(1:4),'int32');
if swapbytes(hk.sizeof_hdr) == 348
   error('%s error: NIfTI image has foreign endian (solution: convert with dcm2nii)',mfilename); 
end
if hk.sizeof_hdr ~= 348
    error('%s error: first byte of NIfTI image should be 348',mfilename);
end
hk.data_type =char(data(5:14));
hk.db_name =char(data(15:32));
hk.extents  = typecast(data(33:36),'int32');
hk.session_error = typecast(data(37:38),'int16');
hk.regular       = char(data(39));
hk.dim_info      = typecast(data(40),'uint8');
%next read dimensions
dime.dim        = typecast(data(41:56),'int16')';
dime.intent_p1  = typecast(data(57:60),'single')';
dime.intent_p2  = typecast(data(61:64),'single')';
dime.intent_p3  = typecast(data(65:68),'single')';
dime.intent_code= typecast(data(69:70),'int16')';
dime.datatype   = typecast(data(71:72),'int16')';
dime.bitpix     = typecast(data(73:74),'int16')';
dime.slice_start= typecast(data(75:76),'int16')';
dime.pixdim     = typecast(data(77:108),'single')';
dime.vox_offset = typecast(data(109:112),'single')';
dime.scl_slope  = typecast(data(113:116),'single')';
dime.scl_inter  = typecast(data(117:120),'single')';
dime.slice_end  = typecast(data(121:122),'int16')';
dime.slice_code = typecast(data(123),'uint8');
dime.xyzt_units = typecast(data(124),'uint8');
dime.cal_max    = typecast(data(125:128),'single')';
dime.cal_min    = typecast(data(129:132),'single')';
dime.slice_duration= typecast(data(133:136),'single')';
dime.toffset    = typecast(data(137:140),'single')';
dime.glmax      = typecast(data(141:144),'int32')';
dime.glmin      = typecast(data(145:148),'int32')';
%read history
hist.descrip     = char(data(149:228));
hist.aux_file    = char(data(229:252));
hist.qform_code  = typecast(data(253:254),'int16')';
hist.sform_code  = typecast(data(255:256),'int16')';
hist.quatern_b   = typecast(data(257:260),'single')';
hist.quatern_c   = typecast(data(261:264),'single')';
hist.quatern_d   = typecast(data(265:268),'single')';
hist.qoffset_x   = typecast(data(269:272),'single')';
hist.qoffset_y   = typecast(data(273:276),'single')';
hist.qoffset_z   = typecast(data(277:280),'single')';
hist.srow_x      = typecast(data(281:296),'single')';
hist.srow_y      = typecast(data(297:312),'single')';
hist.srow_z      = typecast(data(313:328),'single')';
hist.intent_name = char(data(329:344));
hist.magic       = char(data(345:347))';
if ~strcmp(hist.magic, 'n+1') && ~strcmp(hist.magic, 'ni1')
    hist.qform_code = 0;
    hist.sform_code = 0;
end %old analyze format image
hist.originator  = typecast(data(253:262),'int16')'; %used by SPM2 and earlier
h.hk = hk; h.dime = dime; h.hist = hist;
%end readHdrSub()