function [slab] = MIP_zarr_merge_block_mips(flagFullpaths, axisInd, imSize, batchSize, dtype, maxMem)
% merge the MIPs of the blocks that the block tasks writing a zarr saved in
% their flags (MIP_zarr_save_block_mips) into the MIP slab of XR_MIP_zarr
% along axis axisInd, so that the zarr isn't read again: element k of the
% slab along the axis is the max over batch k of the zarr. The parts of
% blocks in the same batch are merged in their order along the axis, as max
% keeps the first of equal values (0 and -0) and the first NaN. Returns [] if
% a flag has no MIPs (e.g., from a task without them) or MIPs for other
% batches, if the blocks don't cover the zarr once, or if the merge would
% hold more than maxMem; then the caller reads the zarr.
%
% flagFullpaths: flag files of all block tasks that wrote the zarr
% axisInd: axis of the slab, 1, 2, 3 for y, x, z
% imSize: size of the zarr
% batchSize: batch size of the MIP slabs
% dtype: data type of the zarr
% maxMem: memory the merge may hold in GB (optional, default and NaN: no limit)


if nargin < 6
    maxMem = Inf;
end

axis_strs = {'y', 'x', 'z'};
varName = [axis_strs{axisInd}, 'MIPs'];
slabSize = imSize;
slabSize(axisInd) = ceil(imSize(axisInd) / batchSize(axisInd));

% flags without MIPs are expected (the zarr is read then)
warnState = warning('off', 'MATLAB:load:variableNotFound');
restoreWarnings = onCleanup(@() warning(warnState));

% the blocks of all flags, which must be saved for these batches
slab = [];
nF = numel(flagFullpaths);
bboxCell = cell(nF, 1);
for f = 1 : nF
    try
        a = load(flagFullpaths{f}, 'blockBboxes', 'mipBatchSize');
    catch
        return;
    end
    if ~all(isfield(a, {'blockBboxes', 'mipBatchSize'})) || ~isequal(a.mipBatchSize(:)', batchSize(:)')
        return;
    end
    bboxCell{f} = a.blockBboxes;
end
blockBboxes = cat(1, bboxCell{:});

% the blocks lie in the zarr and add up to its size (they don't overlap)
if isempty(blockBboxes) || any(blockBboxes(:, 1 : 3) < 1, 'all') || any(blockBboxes(:, 4 : 6) > imSize, 'all') ...
        || sum(prod(blockBboxes(:, 4 : 6) - blockBboxes(:, 1 : 3) + 1, 2)) ~= prod(imSize)
    return;
end

% the batches along the axis that each block overlaps, and the memory of
% the slab, the coverage count and the MIPs of the blocks
k1 = floor((blockBboxes(:, axisInd) - 1) / batchSize(axisInd)) + 1;
k2 = floor((blockBboxes(:, axisInd + 3) - 1) / batchSize(axisInd)) + 1;
mipSizes = blockBboxes(:, 4 : 6) - blockBboxes(:, 1 : 3) + 1;
mipSizes(:, axisInd) = k2 - k1 + 1;
byteNum = dataTypeToByteNumber(dtype);
if (prod(slabSize) * (byteNum + 1) + sum(prod(mipSizes, 2)) * byteNum) / 2^30 > maxMem
    return;
end

mipCell = cell(nF, 1);
for f = 1 : nF
    try
        a = load(flagFullpaths{f}, varName);
    catch
        return;
    end
    if ~isfield(a, varName) || numel(a.(varName)) ~= size(bboxCell{f}, 1)
        return;
    end
    mipCell{f} = a.(varName);
end
mips = cat(1, mipCell{:});
clear mipCell bboxCell a;

% the MIPs of the blocks by batch along the axis, with their region in the
% slab and the start of their block along the axis
nPart = sum(k2 - k1 + 1);
parts = cell(nPart, 1);
partBboxes = zeros(nPart, 6);
partStarts = zeros(nPart, 1);
p = 0;
idx = {':', ':', ':'};
for b = 1 : size(blockBboxes, 1)
    if ~isa(mips{b}, dtype) || ~isequal(size(mips{b}, 1 : 3), mipSizes(b, :))
        return;
    end
    for k = k1(b) : k2(b)
        p = p + 1;
        if k1(b) == k2(b)
            parts{p} = mips{b};
        else
            idx{axisInd} = k - k1(b) + 1;
            parts{p} = mips{b}(idx{:});
        end
        partBboxes(p, :) = blockBboxes(b, :);
        partBboxes(p, [axisInd, axisInd + 3]) = k;
        partStarts(p) = blockBboxes(b, axisInd);
    end
    mips{b} = [];
end

% parts of blocks in the same batch have the same region; the regions cover
% the slab once
[groupBboxes, ~, ic] = unique(partBboxes, 'rows');
covered = zeros(slabSize, 'uint8');
for g = 1 : size(groupBboxes, 1)
    gb = groupBboxes(g, :);
    covered(gb(1) : gb(4), gb(2) : gb(5), gb(3) : gb(6)) = covered(gb(1) : gb(4), gb(2) : gb(5), gb(3) : gb(6)) + 1;
end
if ~isequal(size(covered, 1 : 3), slabSize) || any(covered(:) ~= 1)
    return;
end
clear covered;

% merge the parts of each region in their order along the axis
slab = zeros(slabSize, dtype);
[~, order] = sortrows([ic, partStarts]);
for i = 1 : nPart
    p = order(i);
    if i == 1 || ic(p) ~= ic(order(i - 1))
        mip = parts{p};
    else
        mip = max(mip, parts{p});
    end
    parts{p} = [];
    if i == nPart || ic(order(i + 1)) ~= ic(p)
        gb = groupBboxes(ic(p), :);
        slab(gb(1) : gb(4), gb(2) : gb(5), gb(3) : gb(6)) = mip;
    end
end

end
