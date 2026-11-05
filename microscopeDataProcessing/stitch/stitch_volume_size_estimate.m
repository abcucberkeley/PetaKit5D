function [nvSize, numBatches, nvSizeRange, numBatchesRange] = stitch_volume_size_estimate(coordinates, imSize, varargin)
% the size of the volume XR_stitching_frame_zarr_dev_v1 stitches from tiles
% of size imSize (as read; all tiles alike) at the stage coordinates (n x 3,
% um), without the shifts of the registration, and its number of batches;
% nvSizeRange (smallest; largest, 2 x 3) and numBatchesRange (fewest, most):
% with the registration (xcorrShift, more than one tile), whose shifts of up
% to its max offsets (voxels) can make the volume smaller or larger by up to
% them at both ends; for sizing the stitching job of the volume
% (XR_matlab_stitching_wrapper). The geometry is the stitching function's:
% axis and data orders, pixel sizes (resample type or factor, skew angle,
% scan type, DS/DSR), the crops and resampling of the tiles,
% stitch_process_coordinates, the padding of primary channels with
% xcorrShift, outBbox and the batch size (keep it in step with it; a test
% compares the two).


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('coordinates', @isnumeric);
ip.addRequired('imSize', @isnumeric);
ip.addParameter('tileIdx', [], @isnumeric);
ip.addParameter('flippedTile', [], @(x) isempty(x) || all(islogical(x) | isnumeric(x)));
ip.addParameter('axisOrder', 'x,y,z', @ischar);
ip.addParameter('dataOrder', 'y,x,z', @ischar);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('IOScan', false, @islogical);
ip.addParameter('DS', false, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('resampleType', 'xy_isotropic', @ischar);
ip.addParameter('resampleFactor', [], @isnumeric);
ip.addParameter('InputBbox', [], @isnumeric);
ip.addParameter('tileOutBbox', [], @isnumeric);
ip.addParameter('xcorrShift', true, @islogical);
ip.addParameter('xyMaxOffset', 300, @isnumeric);
ip.addParameter('zMaxOffset', 50, @isnumeric);
ip.addParameter('stitchPadSize', [2, 2, 2], @isnumeric);
ip.addParameter('outBbox', [], @isnumeric);
ip.addParameter('stitchMIP', [], @(x) isempty(x) || islogical(x));
ip.addParameter('batchSize', [500, 500, 500], @isnumeric);
ip.addParameter('blockSize', [500, 500, 500], @isnumeric);
ip.parse(coordinates, imSize, varargin{:});

pr = ip.Results;
nF = size(coordinates, 1);
DSR = pr.DSR && ~pr.IOScan;
imSize = [imSize(:)', ones(1, 3 - numel(imSize))];

% axis and data orders
axisOrder = strrep(pr.axisOrder, ' ', '');
order_sign_mat = zeros(2, 3);
axisOrder_pure = strrep(strrep(axisOrder, ',', ''), '-', '');
xyz_str = 'xyz';
for i = 1 : 3
    [~, order_sign_mat(1, :)] = sort(axisOrder_pure);
    order_sign_mat(2, i) = 1 - 2 * contains(axisOrder, ['-', xyz_str(i)]);
end
dataOrder_pure = strrep(strrep(pr.dataOrder, ',', ''), '-', '');
[~, data_order_mat] = sort(dataOrder_pure);
[~, data_order_reverse_mat] = sort(data_order_mat);

xyz = coordinates(:, order_sign_mat(1, :)) .* order_sign_mat(2, :);
xyz = xyz - min(xyz, [], 1);

tileIdx = pr.tileIdx;
if isempty(tileIdx)
    tileIdx = zeros(nF, 4);
end
tileNum = [numel(unique(tileIdx(:, 1))), numel(unique(tileIdx(:, 2))), numel(unique(tileIdx(:, 3)))];
tileNum = tileNum(order_sign_mat(1, :));

% pixel sizes of the stitching space
if pr.objectiveScan || pr.IOScan
    zAniso = pr.dz / pr.xyPixelSize;
else
    zAniso = sind(pr.SkewAngle) * pr.dz / pr.xyPixelSize;
end
theta = pr.SkewAngle * pi / 180;
resampleType = pr.resampleType;
resampleFactor = pr.resampleFactor;
if ~isempty(resampleFactor)
    resampleType = 'given';
end
switch resampleType
    case 'given'
        xf = resampleFactor(1);
        yf = resampleFactor(2);
        zf = resampleFactor(3);
    case 'method_1'
        zf = cot(abs(theta));
        yf = 1;
        xf = 1 / cos(abs(theta));
    case 'isotropic'
        zf = 1;
        yf = 1;
        xf = 1;
    case 'xy_isotropic'
        zf = sqrt((sin(theta) ^ 2 + zAniso ^ 2 * cos(theta) ^ 2) / (cos(theta) ^ 2 + zAniso ^ 2 * sin(theta) ^ 2));
        yf = 1;
        xf = 1;
end
if ~DSR
    zf = zAniso * zf;
end
xp = pr.xyPixelSize * xf;
yp = pr.xyPixelSize * yf;
zp = pr.xyPixelSize * zf;

% the tiles as stitched: cropped (InputBbox, then tileOutBbox) and, for
% skewed space, resampled
if isempty(resampleFactor)
    resampleFactor = [1, 1, 1];
end
if numel(pr.InputBbox) == 6
    imSize = pr.InputBbox(4 : 6) - pr.InputBbox(1 : 3) + 1;
end
if numel(pr.tileOutBbox) == 6
    imSize = pr.tileOutBbox(4 : 6) - pr.tileOutBbox(1 : 3) + 1;
end
if ~DSR && ~all(resampleFactor == 1)
    rs = resampleFactor(:)';
    rs = [ones(1, 4 - numel(rs)) * rs(1), rs(2 : end)];
    imSize = round(imSize ./ rs);
end
imSizes = repmat(imSize, nF, 1);

xyz = stitch_process_coordinates(xyz, imSizes, pr.flippedTile, order_sign_mat, data_order_mat, ...
    xp, yp, zp, theta, pr.IOScan, pr.objectiveScan, pr.DS, DSR);

% the stitched volume, padded for the volumes of all time points and
% channels to be of the same size
if nF > 1 && pr.xcorrShift
    pad_size = [xp, yp, zp] .* pr.stitchPadSize(data_order_mat) .* tileNum;
else
    pad_size = [0, 0, 0];
end
origin = min(xyz, [], 1) - pad_size;
dxyz = max(xyz + imSizes(:, data_order_mat) .* [xp, yp, zp], [], 1) + pad_size - origin;
vxyz = round(dxyz ./ [xp, yp, zp] + 2);
nvSize = vxyz(data_order_reverse_mat);
if any(pr.stitchMIP) || imSize(data_order_mat(3)) == 1
    nvSize(data_order_mat(3)) = 1;
end
outBbox = pr.outBbox;
if numel(outBbox) == 6
    outBbox = outBbox(:)';
    bboxEnd = outBbox(4 : 6);
    bboxEnd(isinf(bboxEnd)) = nvSize(isinf(bboxEnd));
    nvSize = bboxEnd - outBbox(1 : 3) + 1;
end

% (with registration: the tiles can move by up to the max offsets, x, y
% and z)
nvSizeRange = [nvSize; nvSize];
if nF > 1 && pr.xcorrShift && numel(outBbox) ~= 6
    maxOffset = [pr.xyMaxOffset(1), pr.xyMaxOffset(end), pr.zMaxOffset(1)];
    margin = 2 * maxOffset(data_order_reverse_mat);
    nvSizeRange = [max(nvSize - margin, 1); nvSize + margin];
    if any(pr.stitchMIP) || imSize(data_order_mat(3)) == 1
        nvSizeRange(:, data_order_mat(3)) = 1;
    end
end
numBatches = batch_number(nvSize, imSize, pr.batchSize, pr.blockSize);
numBatchesRange = [batch_number(nvSizeRange(1, :), imSize, pr.batchSize, pr.blockSize), ...
    batch_number(nvSizeRange(2, :), imSize, pr.batchSize, pr.blockSize)];

end


function [numBatches] = batch_number(nvSize, imSize, batchSize, blockSize)
% the batches of a volume, as the stitching function sets them

batchSize = min(nvSize, batchSize);
batchSize = min(imSize * 2, batchSize);
blockSize = min(batchSize, blockSize);
batchSize = ceil(batchSize ./ blockSize) .* blockSize;
numBatches = prod(ceil(nvSize ./ batchSize));

end
