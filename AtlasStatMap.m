function AtlasStatMap(template, outname, indices, intensities);
%Create a MZ3 file where each region has specified intensity
% template : name of indexed mz3 file
% outname : name of created image
% intensities : brightness for each region
% indices : (optional) matches intensity to region number
%Examples
% AtlasStatMap('jhu.mz3','out.mz3',[], [3 2.5 5]); %create image where region 1 has intensity 3, region 2 has intensity 2.5 - jhu.mz3 MUST have precisely 3 regions
% AtlasStatMap('jhu.mz3','out.mz3',[17 22 32], [3.1 4 5]); %create image where region 17 has intensity 3.1 - jhu.mz3 MUST have at leat 32 regions

if isempty(which('fileUtils.mz3.readMz3')), error('Please install MRIcroS https://github.com/bonilhamusclab/MRIcroS'); end;
[f, v, c] = fileUtils.mz3.readMz3(template, false);
if size(c,2) ~= 5
   error('%s is not a template (it should have RGB and Scalar values\n', template);
end
cIndex = c(:,5);
if max(cIndex(:)) < max(indices(:))
    error('%s only has %d regions\n', max(cIndex(:)) );
end
cOut = zeros(size(v,1),1);%output color
cOut(:) = NaN;
for i = 1 : numel(indices)
    cOut(cIndex == indices(i)) = intensities(i);
end
%check that there is something to write
nSurvive = sum(isfinite(cOut));
if nSurvive == 0 %we have some undefined vertices - remove these
    error('Unable to find any of the regions in this template\n');
end
%write output
if nSurvive ~= numel(cOut) %not all vertices survive - remove unused vertices
    %collapse mesh to remove unused vertices
    %collapse faces
    vIdx = zeros(numel(cOut),1);
    vIdx(isfinite(cOut)) = [1:nSurvive];
    fOut = vIdx(f(:));
    fOut = reshape(fOut,size(f));
    fOut = fOut(min(fOut') ~= 0,:); %#ok<UDIM>
    %collapse vertices
    vOut = v(isfinite(cOut),:);
    %collapse colors
    cOut = cOut(isfinite(cOut));
else
    %all vertices survive - no need to collapse
    fOut = f;
    vOut = v;
end
fileUtils.mz3.writeMz3(outname, fOut, vOut, cOut);
%end AtlasStatMap()
