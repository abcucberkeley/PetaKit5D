function [batchSize] = MIP_zarr_batch_size(inSize, blockSize, dtype, varargin)
% batch size of the large-file MIPs (XR_MIP_zarr) of a zarr region of size
% inSize with blocks (zarr chunks or pooling blocks) of blockSize: the
% requested batch size rounded up to whole blocks, at most inSize, and for
% the final MIPs (not mipSlab) increased by a block along the axis of the
% largest intermediate MIP slab while a slab is over 100 GB (and the batch
% under 100 GB). The block functions of large-file processing save the MIPs
% of their blocks for these batches.
%
% Parameters:
%   batchSize : 1x3 vector (default: [2048, 2048, 2048]). Requested batch size in y, x, z.
%     mipSlab : true|false (default: false). Only the MIP slabs (no final MIPs), no increase for large slabs.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inSize', @isnumeric);
ip.addRequired('blockSize', @isnumeric);
ip.addRequired('dtype', @ischar);
ip.addParameter('batchSize', [2048, 2048, 2048], @isnumeric); % in y, x, z
ip.addParameter('mipSlab', false, @islogical);

ip.parse(inSize, blockSize, dtype, varargin{:});

pr = ip.Results;
batchSize = pr.batchSize;
mipSlab = pr.mipSlab;

byteNum = dataTypeToByteNumber(dtype);

batchSize = min(inSize, max(batchSize, ceil(batchSize ./ blockSize) .* blockSize));
if mipSlab
    return;
end

% in case of out size too large for very large data that causes oom for the main
% job, increase batch size if necessary.
bSubSz = ceil(inSize ./ batchSize);
outVolSizes = zeros(3, 1);
for i = 1 : 3
    outVolSizes(i) = prod([inSize(setdiff(1 : 3, i)), bSubSz(i)]) * byteNum / 1024^3;
end

% if the max intermediate MIP files is greater than 100 GB, increase the batchSize
% by the blockSize in the axis with largest outVolSizes
while any(outVolSizes > 100) && prod(batchSize) * byteNum / 1024^3 < 100
    [~, ind] = max(outVolSizes);
    batchSize = min(inSize, batchSize + blockSize .* ((1 : 3) == ind));
    bSubSz = ceil(inSize ./ batchSize);

    for i = 1 : 3
        outVolSizes(i) = prod([inSize(setdiff(1 : 3, i)), bSubSz(i)]) * byteNum / 1024^3;
    end
end

end
