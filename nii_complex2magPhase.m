function nii_complex2magPhase = nii_smooth (realImg, imagImg)
%Given MRI Real and Imaginary images creates phase and magnitude
%  realImg : Real image in NIfTI format 
%  imageImg : Real image in NIfTI format
% Examples
%  nii_complex2magPhase ('real.nii','image.nii');
%  nii_complex2magPhase; %use GUI

isMagUint16 = true; %if true, data saved as 16-bit unsigned integer, else 32-bit real
if ~exist('realImg', 'var')
 realImg = spm_select(inf,'image','Select real image');
end;
if ~exist('imagImg', 'var')
 imagImg = spm_select(inf,'image','Select imaginary image');
end;
[pth,nm,ext] = spm_fileparts(realImg);

rhdr = spm_vol(realImg);
rimg = spm_read_vols(rhdr);
ihdr = spm_vol(imagImg);
iimg = spm_read_vols(ihdr);
%create magniude image
mimg = sqrt(rimg.^2 + iimg.^2);
mhdr = rhdr;
mhdr.pinfo = [1;0;0]; %slope=1, intercept=0 
if isMagUint16
    mhdr.dt(1) = 512; %16-bit uint
    mimg = mimg - min(mimg(:)); %set minimum to zero
    scalef = spm_type(mhdr.dt(1),'maxval')/max(mimg(:));
    mimg = mimg * scalef;
    fprintf('Magnitude image saved as UINT16, range %g..%g\n', min(mimg(:)), max(mimg(:)) );
else
    mhdr.dt(1) = 16; %32-bit float
    fprintf('Magnitude image range %g..%g\n', min(mimg(:)), max(mimg(:)) );
end
mhdr.fname = fullfile(pth, [nm '_mag' ext]);
spm_write_vol(mhdr,mimg);

%create phse image
pimg = atan2(iimg,rimg); %phase angle
phdr = rhdr;
phdr.dt(1) = 16; %32-bit float
phdr.pinfo = [1;0;0]; %slope=1, intercept=0 
phdr.fname = fullfile(pth, [nm '_ph' ext]);  
spm_write_vol(phdr,pimg);
fprintf('Phase image range (-pi..pi) %g..%g\n', min(pimg(:)), max(pimg(:)) );


