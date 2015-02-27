function fnm = convert_ctpSub(fnm, modalityIsTime)
%Convert Siemens RGB CT-perfusion images to linear grayscale images
% V : name[s] of image[s] to convert (optional)
% modalityIsTime : Siemens uses different RGB schemes 
%              TRUE for MTT or TTP images (transit time, time to peak)
%              FALSE for CBF or CBV images (cerebral blow flow, volume)
%Examples
% convert_ctp('MTT.nii',true);
% convert_ctp('CBF.nii',false);
% convert_ctp(strvcat('MTT.nii','TTP.nii'),true);
% convert_ctp(strvcat('CBF.nii','CBV.nii'),false);
if ~exist('fnm','var') || isempty(fnm), return; end;
%if ~exist('fnm','var'), fnm = spm_select(1,'image','Select image to convert from RGB to scalar'); end
if ~exist('modalityIsTime','var'),modalityIsTime = true; end;
%fprintf('rgb->scalar %s\n',fnm);
fnm = deblank(fnm);
[pth nam ext] = spm_fileparts(fnm);
hdr = spm_vol(fnm);
if hdr.dt(1) ~= 128
fprintf('%s error: input file must be RGB format %s\n',mfilename,fnm);
return;
end
hdrM = hdr;
hdrM.dt(1)    =2; %2=8-bit char, 4= 16-bit integer; 16 =32-bit real datatype
hdrM.pinfo(3) = 352;
hdrM.dim(3) = hdrM.dim(3) * 3; %red,green,blue each saved on separate slice
imgRGB = spm_read_vols(hdrM);
imgR = imgRGB(:,:,(1 : 3: hdrM.dim(3))); %red plane
imgG = imgRGB(:,:,(2 : 3: hdrM.dim(3))); %green plane
imgB = imgRGB(:,:,(3 : 3: hdrM.dim(3))); %blue plane
if modalityIsTime == 1
img = ttp_rgb2scalarSub(imgR(:)', imgG(:)', imgB(:)');
else
img = cbf_rgb2scalarSub(imgR(:)', imgG(:)', imgB(:)');
img = img * 10; %we are saving as integers, so preserve precision
end
img = reshape(img,hdr.dim);
hdrOut = hdr;
hdrOut.fname = fullfile(pth, ['q' nam ext]);
hdrOut.dt(1)    =4; %2=8-bit char, 4= 16-bit integer; 16 =32-bit real datatype
hdrOut.pinfo(3) = 352;
spm_write_vol(hdrOut,img);
fnm = hdrOut.fname;
%end convert_ctpSub()

function s = ttp_rgb2scalarSub (R,G,B)
%converts Siemens RGB color scheme for MTT/TTP to scalar intensity
s = zeros(numel(R),1);
%next: voxels where Red is max 765...1020
v = 1020 - G;
idx =  (R == 255  ); 
s(idx) = v(idx);
%next voxels where red > blue and red < 255, 511..764
v = R + 510;
idx = intersect (find (R > B), find(R < 255 ) ); %only items 22:39
s(idx) = v(idx);
%next voxels where blue > red
v = (G+1)*2; %range 2..510
idx =intersect (find (G > 0), find (B >= R)); %only items 22:39
s(idx) = v(idx);
%end ttp_rgb2scalarSub()

function [s] = cbf_rgb2scalarSub (R,G,B)
%converts Siemens RGB color scheme for CBF/CBV to scalar intensity
%Siemens Red,Green,Blue colorscale has 127 values, but images have interpolated intermediate values
% Rv =[0 32 64 67 70 73 76 79 82 85 88 92 96 92 88 85 82 79 76 73 70 67 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 65 66 68 69 71 73 74 76 77 79 81 82 84 85 87 89 90 92 93 95 97 98 100 101 103 105 106 108 109 111 113 114 116 117 119 121 122 124 125 127 131 134 138 141 145 149 152 156 159 163 167 170 174 177 181 185 188 192 195 199 203 206 210 213 217 221 224 228 231 235 239 242 246 249 253 253 251 248 246 243 240 238 235 233 230 227 225];
% Gv =[0 0 0 0 0 0 0 0 0 0 0 0 0 6 12 19 25 32 38 45 51 58 63 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 68 73 77 82 87 92 97 101 106 111 118 121 125 130 135 140 145 149 154 159 164 169 173 178 183 188 193 197 202 207 212 217 221 226 231 236 241 245 250 254 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 255 245 226 207 188 169 150 131 112 93 74 55 41];
% Bv =[0 32 64 67 70 73 76 79 82 85 88 92 96 99 102 105 108 111 114 117 120 124 128 134 141 148 155 162 169 176 183 190 197 204 211 218 225 232 239 246 254 250 243 237 230 224 218 211 205 198 192 186 179 173 166 160 154 147 141 134 128 122 115 109 102 96 90 83 77 70 64 58 51 45 38 32 26 19 13 6 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 4 6 9 11 14 17 19 22 24 27 30];
s = zeros(numel(R),1);
%Indices 1:21 
v = (B-60)/3.2;
idx = intersect (find (B >= R), find(R > G ) ); %only items 1:21
s(idx) = v(idx);
%indices 22:39
v = (B-128)/7+21.14285;
idx = intersect (find (B > R), find(R == G ) ); %only items 22:39
s(idx) = v(idx);
%indices 40:79
v = (G-64)/4.75+39;
idx = intersect (find (G < 255), find(G > R ) ); %only items 40:79
s(idx) = v(idx);
%indices 80:114
v = (R-127)/3.6+79;
idx = find (G >= 255); %only items 80:114
s(idx) = v(idx);
%indices 115:126
v = (B-1)/2.636363636+115;
idx = intersect (find (R > G), find(R > B ) ); %only items 115:126
s(idx) = v(idx);
s(s<0) = 0;
%end cbf_rgb2scalarSub
