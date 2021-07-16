function [xmin,xmax,ymin,ymax,zmin,zmax,nn] = XR_subVolumeCoordinatesExtraction(volSize, varargin)
% function generates sub volumes for given volume, with approximate the
% same volumes
% 
% The idea is to try to split y, x, and z size to the sizeRang. If there are 
% multiple splitting methods, first keep the partition number as small as possible. 
% Then check if the total volume is within the maxSubVolume. If not, increase the
% partition number from the maximum sub range. 
%
% Xiongtao Ruan, 03/09/2020
% xruan (03/10/2020): add option to ensure the size of chunks nearly a
% number with only factors of 2, 3, 5, 7 (to avoid to much cropping in the deconvolution).

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('volSize', @(x) isvector(x) && numel(x) == 3);
ip.addParameter('ChunkSize', [1024,1024,1024], @(x) isvector(x) && numel(x) == 3 && all(x >= 1)); % in y, x, z
ip.addParameter('overlapSize', 50, @(x) isnumeric(x) && numel(x) <= 3);
ip.addParameter('maxSubVolume', 1.2e9, @isnumeric); % about 1GB
ip.addParameter('sizeRange', [0.75, 1.5], @isnumeric); % [minimum, maximum] size. 
ip.addParameter('useGoodNum', true, @isnumeric); % [minimum, maximum] size. 

ip.parse(volSize, varargin{:});

if any(volSize < 2)
    error('The volume size [%s] is not valid for a 3d volume.', num2str(volSize));
end

p = ip.Results;
ChunkSize = p.ChunkSize;
overlapSize = p.overlapSize;
maxSubVolume = p.maxSubVolume;
sizeRange = p.sizeRange;
useGoodNum = p.useGoodNum;

if numel(overlapSize) == 1
    overlapSize = ones(1, 3) * overlapSize;
elseif numel(overlapSize) == 2
    overlapSize = [overlapSize(1), overlapSize(1), overlapSize(2)];
end

volBlockRatio = volSize ./ ChunkSize;
minPartitionNum = ceil(volBlockRatio / sizeRange(2));
maxPartitionNum = ceil(volBlockRatio / sizeRange(1));

chunkSize_cell = cell(3, 1);
chunkRatio_cell = cell(3, 1);
partititionNum_cell = cell(3, 1);
for i = 1 : 3
    partitionNumRange = minPartitionNum(i) : maxPartitionNum(i);
    chunkSize_cell{i} = ceil(volSize(i) ./ partitionNumRange + overlapSize(i));
    chunkRatio_cell{i} = chunkSize_cell{i} ./ ChunkSize(i);
    partititionNum_cell{i} = partitionNumRange;
end

% start from the largest size to smallest size and check volumes
partitionNums = cellfun(@numel, chunkSize_cell);
inds = [1, 1, 1];
chosenChunkSize = [chunkSize_cell{1}(inds(1)), chunkSize_cell{2}(inds(2)), chunkSize_cell{3}(inds(3))];    

while prod(chosenChunkSize) > maxSubVolume
    partitionMaxRatios = arrayfun(@(x) chunkRatio_cell{x}(min(inds(x) + 1, partitionNums(x))), 1 : 3);
    [~, inds_sorted] = sort(partitionMaxRatios, 'descend');
    ct = 1;
    while inds(inds_sorted(ct)) >= numel(chunkSize_cell{inds_sorted(ct)})
        ct = ct + 1;
    end
    increaseInd = inds_sorted(ct);
    inds(increaseInd) = inds(increaseInd) + 1;
    chosenChunkSize = [chunkSize_cell{1}(inds(1)), chunkSize_cell{2}(inds(2)), chunkSize_cell{3}(inds(3))];    
end

chosenPartitionNum = [partititionNum_cell{1}(inds(1)), partititionNum_cell{2}(inds(2)), partititionNum_cell{3}(inds(3))];

% actual proposed chunk size.
csz = ceil((volSize + overlapSize .* (chosenPartitionNum - 1)) ./ chosenPartitionNum);

if useGoodNum
    for i = 1 : 3
        csz(i) = findGoodFactorNumber(csz(i), 1);
    end
end

% actual overlap in the left side.
aol = ceil((csz .* chosenPartitionNum - volSize) ./ (chosenPartitionNum - 1));
aol(chosenPartitionNum == 1) = 0;
% in some cases, the actual size is not multiple for csz, in this case we
% may need to shift some end chunks by some pixels. 
res = volSize - ((csz - aol) .* chosenPartitionNum + aol);
res(chosenPartitionNum == 1) = 0;

yl = 1 : csz(1) - aol(1) : volSize(1);
yl = yl(1 : chosenPartitionNum(1));
yl(end - res(1) + 1 : end) = yl(end - res(1) + 1 : end) + (1 : res(1));
yr = min(volSize(1), yl + csz(1) - 1);

xl = 1 : csz(2) - aol(2) : volSize(2);
xl = xl(1 : chosenPartitionNum(2));
xl(end - res(2) + 1 : end) = xl(end - res(2) + 1 : end) + (1 : res(2));
xr = min(volSize(2), xl + csz(2) - 1);

zl = 1 : csz(3) - aol(3) : volSize(3);
zl = zl(1 : chosenPartitionNum(3));
zl(end - res(3) + 1 : end) = zl(end - res(3) + 1 : end) + (1 : res(3));
zr = min(volSize(3), zl + csz(3) - 1);

[Ymin, Xmin, Zmin] = meshgrid(yl, xl, zl);
[Ymax, Xmax, Zmax] = meshgrid(yr, xr, zr);

ymin = Ymin(:);
ymax = Ymax(:);
xmin = Xmin(:);
xmax = Xmax(:);
zmin = Zmin(:);
zmax = Zmax(:);

nn = numel(ymin);

end

