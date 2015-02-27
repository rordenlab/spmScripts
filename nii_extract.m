function nii_extract (V,ClusterMM, V2)
%Extract object: remove noise from air.
%Runs Chris Rorden's Object Extraction, 2011 BSD license
%runs object extraction with default parameters
%  run executable from command line to learn about adjustable parameters
%  http://www.mccauslandcenter.sc.edu/mricrogl/?page_id=79
if nargin<1, V = spm_select(inf,'image','Select image[s] to extract'); end
if nargin<2, ClusterMM = 0; end
for i=1:size(V,1)
  r = deblank(V(i,:));
  [pth,nam,ext, vol] = spm_fileparts(r);
  ref = fullfile(pth,[nam ext]);
  [pth,nam,ext] = fileparts(mfilename('fullpath'));
	%OUTPUT OF COMPUTER	VERSION OF MATLAB (32 bit / 64 bit)
	%PCWIN			32 bit MATLAB on Windows
	%PCWIN64			64 bit MATLAB on Windows
	%GLNX86			32 bit MATLAB on Linux
	%GLNXA64			64 bit MATLAB on Linux
	%MACI64			64 bit MATLAB on Apple Mac OS X
  t = computer;
  switch t
	case 'PCWIN'
		exename = 'extract.exe';
	case 'PCWIN64'
		exename = 'extract.exe';
	case 'GLNX86'
		exename = 'extractlx32';
	case 'GLNXA64'
		exename = 'extractlx64';
	case 'MACI32'
		exename = 'extractosx';
	case 'MACI64'
		exename = 'extractosx';
	otherwise
		warning('Unexpected (probably unsupported) computer type.');
  end;
  if nargin>2, 
    r = deblank(V2(i,:));
      ref = [pth,filesep,exename,' -c ',num2str(ClusterMM) ,' ', ref, ' ',r];
  else
    ref = [pth,filesep,exename,' -c ',num2str(ClusterMM),' ', ref];
  end
  system(ref);
end;