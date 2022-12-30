function [batchBBoxes, regionBBoxes, localBBoxes] = XR_zarrChunkCoordinatesExtraction(volSize, varargin)
% generic framework for chunking zarr files in batchs (factored by batch
% size) with borderSize
% The outputs are batch bounding boxes and the effective output region
% bounding box
%
% Author: Xiongtao Ruan (11/11/2021)
% 
% xruan (02/18/2022): add support for bbox for input, input batch can go
% beyond bbox if bordersize is provided. 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('volSize', @(x) isvector(x) && numel(x) == 3);
ip.addParameter('BatchSize', [1024,1024,1024], @(x) isvector(x) && numel(x) == 3 && all(x >= 1)); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @(x) isnumeric(x) && numel(x) <= 3); 
ip.addParameter('bbox', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 6)); % bounding box for which region to process
ip.addParameter('SameBatchSize', true, @(x) islogical(x)); % same batch size before adding borderSize
ip.addParameter('BorderSize', 50, @(x) isnumeric(x) && numel(x) <= 3);

ip.parse(volSize, varargin{:});

if any(volSize < 2)
    error('The volume size [%s] is not valid for a 3d volume.', num2str(volSize));
end

pr = ip.Results;

BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
SameBatchSize = pr.SameBatchSize;
BorderSize = pr.BorderSize;
bbox = pr.bbox;

% use world start and end to unify the cases with or without bounding box
if isempty(bbox)
    wdStart = [1, 1, 1];
    wdEnd = volSize;    
else
    wdStart = bbox(1 : 3);
    wdEnd = bbox(4 : 6);
end

wdSize = wdEnd - wdStart + 1;

BatchSize = min(wdSize, BatchSize);
BlockSize = min(wdSize, BlockSize);

% normalize batch size so it contains complete blocksizes 
BatchSize = ceil(BatchSize ./ BlockSize) .* BlockSize;

bSubSz = ceil(wdSize ./ BatchSize);
numBatch = prod(bSubSz);

regionBBoxes = zeros(numBatch, 6);

[Y, X, Z] = ndgrid(1 : bSubSz(1), 1 : bSubSz(2), 1 : bSubSz(3));
bSubs = [Y(:), X(:), Z(:)];
clear Y X Z

regionBBoxes(:, 1 : 3) = (bSubs - 1) .* BatchSize + wdStart; 
regionBBoxes(:, 4 : 6) = min(regionBBoxes(:, 1 : 3) + BatchSize - 1, wdEnd);

% get batch bbox from region bbox
batchBBoxes = regionBBoxes;

if SameBatchSize
    batchBBoxes(:, 1 : 3) = batchBBoxes(:, 4 : 6) - BatchSize + 1;
end

batchBBoxes(:, 1 : 3) = max(1, batchBBoxes(:, 1 : 3) - BorderSize);
batchBBoxes(:, 4 : 6) = min(volSize, batchBBoxes(:, 4 : 6) + BorderSize);

% bbox crop from input to ouput batches
localBBoxes = zeros(numBatch, 6);
localBBoxes(:, 1 : 3) = regionBBoxes(:, 1 : 3) - batchBBoxes(:, 1 : 3) + 1;
localBBoxes(:, 4 : 6) = localBBoxes(:, 1 : 3) + (regionBBoxes(:, 4 : 6) - regionBBoxes(:, 1 : 3));

% reset regionBBoxes to start from bbox
regionBBoxes = regionBBoxes - [wdStart, wdStart] + 1;

end



