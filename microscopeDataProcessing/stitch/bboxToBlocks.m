function [blockInfo] = bboxToBlocks(bbox, bboxStart, blockSize, imageSize, borderSize)
% split bounding box to blocks with defined block size and image size
% 
% bbox: bouding box [ymin, xmin, zmin, ymax, xmax, zmax]
% bboxStart: the star index in the original image/tile [y, x, z]
% blockSize: block size in stitched image
% imageSize: image size of stitched image
% 
% Author: Xiongtao Ruan (10/04/2020)
% 
% xruan (10/17/2020): add border size for extra info for a block
% xruan (10/25/2020): add support for pad/crop (out of border)


if nargin < 5
    borderSize = [0, 0, 0];
end

bSubSz = ceil(imageSize ./ blockSize);
% numBlocks = prod(bSubSz);

bSubs_cell = cell(3, 1);
bCoords_cell = cell(3, 1);
for i = 1 : 3
    s = bbox(i);
    t = bbox(3 + i);
    bsz = blockSize(i);
    nb = bSubSz(i);
    bdsz = borderSize(i);
    
    [bSubs_i, bCoords_i] = singleAxisSplit(s, t, bsz, nb, bdsz);
    bSubs_cell{i} = bSubs_i;
    bCoords_cell{i} = bCoords_i;
end

% block subscripts
[Y, X, Z] = ndgrid(bSubs_cell{1}, bSubs_cell{2}, bSubs_cell{3});
blockSub = [Y(:), X(:), Z(:)];

% block coordinates
[Y, X, Z] = ndgrid(bCoords_cell{1}(1, :), bCoords_cell{2}(1, :), bCoords_cell{3}(1, :));
bCoords_s = [Y(:), X(:), Z(:)];
[Y, X, Z] = ndgrid(bCoords_cell{1}(2, :), bCoords_cell{2}(2, :), bCoords_cell{3}(2, :));
bCoords_t = [Y(:), X(:), Z(:)];
bCoords = [bCoords_s, bCoords_t];

% world coordinates and bonding box coordinates
% wCoords = repmat((blockSub - 1) .* blockSize, 1, 2) + bCoords;
wCoords = repmat((blockSub - 1) .* blockSize - (blockSub ~= 1) .* borderSize, 1, 2)  + bCoords;

bboxCoords = wCoords - repmat(bbox(1 : 3) - bboxStart, 1, 2);

% collect all info
blockInfo.blockSub = blockSub;
blockInfo.blockInd = sub2ind(bSubSz, blockSub(:, 1), blockSub(:, 2), blockSub(:, 3));
blockInfo.bCoords = bCoords;
blockInfo.wCoords = wCoords;
blockInfo.bboxCoords = bboxCoords;
blockInfo.bbox = bbox;
blockInfo.bboxStart = bboxStart;
blockInfo.blockSize = blockSize;
blockInfo.imageSize = imageSize;
blockInfo.bSubSize = bSubSz;
blockInfo.borderSize = borderSize;


end


function [bSubs, bCoords] = singleAxisSplit(s, t, bsz, nb, bdsz)
% split bbox for blocks for a single axis

bSubs = max(1, ceil((s - bdsz) / bsz)) : min(nb, ceil((t + bdsz) / bsz));
if bdsz > 0
    block_bdsz = (bSubs ~= 1) * bdsz;
    bcoords_s = [max(1, s - (bSubs(1) - 1) * bsz + block_bdsz(1)), ones(1, numel(bSubs) - 1)];
    bcoords_t = [ones(1, numel(bSubs) - 1) .* (bsz + bdsz + block_bdsz(1 : end - 1)), min(bsz + 2 * bdsz, t - (bSubs(end) - 1) * bsz + block_bdsz(end))];
else
    bcoords_s = [max(1, s - (bSubs(1) - 1) * bsz), ones(1, numel(bSubs) - 1)];
    bcoords_t = [ones(1, numel(bSubs) - 1) * bsz, min(bsz, t - (bSubs(end) - 1) * bsz)];
end
bCoords = [bcoords_s; bcoords_t];

% wCoords = (bSubs - 1) * bsz + bCoords;
% 
% bboxCoords = wCoords - s + 1;
% 
end
