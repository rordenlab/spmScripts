function nii_make3d_new
%Make a NIfTI image using SPM

dim = [128, 128, 128, 1]; %image resolution in columns, rows, slices, volumes
dtype = 32; %precision of data
ofile    = 'test.nii';%spm_file(parfile,'path',opts.outdir,'ext',opts.ext);
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
mat(1,4) = -dim(1)/2;
mat(2,4) = -dim(2)/2;
mat(3,4) = -dim(3)/2;
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
for i=1:dim(4)
   N.dat(:,:,:,i) = makeBorgSub(dim(1), dim(2), dim(3));
end
%end nii_mk()

function img = makeBorgSub (nX, nY, nZ)
img = zeros(nX,nY,nZ); 
freq =7;
for z = 2 : (nZ-1) 
	z1 = freq*z/nZ;
	for y = 2 : (nY-1) 
		y1 = freq*y/nY;
		for x = 2 : (nX-1) 
			x1 = freq*x/nX;			
			img(x,y,z)=pi+sin(x1*y1)+sin(y1*z1)+sin(z1*x1);
		end;
	end;
end;
%end makeBorgSub()