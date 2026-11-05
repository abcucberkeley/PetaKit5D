function [mem] = MIP_zarr_merge_memory(imSize, regionSize, blockSize, dtype)
% memory (GB) that XR_MIP_zarr holds to merge the MIP slabs of a zarr of
% size imSize from the MIPs of its blocks (MIP_zarr_merge_block_mips, one
% axis at a time): the slab, its coverage count and the MIPs of the blocks
% along the axis, for the axis that needs the most. The blocks are the
% regions of the block tasks, a grid of regionSize; blockSize and dtype are
% the zarr's chunk size and data type (they give the MIP batch size). For a
% region grid along an axis, the MIPs of a column of blocks have at most one
% element per block and per boundary of MIP batches inside it, as the merge
% counts them exactly from the flags.

imSize = imSize(:)';
regionSize = min(imSize, regionSize(:)');
byteNum = dataTypeToByteNumber(dtype);
mipBatchSize = MIP_zarr_batch_size(imSize, blockSize, dtype);

mem = 0;
for i = 1 : 3
    nMIP = ceil(imSize(i) / mipBatchSize(i));
    nRegion = ceil(imSize(i) / regionSize(i));
    plane = prod(imSize) / imSize(i);
    slabVoxels = plane * nMIP;
    mipVoxels = plane * (nRegion + nMIP - 1);
    mem = max(mem, (slabVoxels * (byteNum + 1) + mipVoxels * byteNum) / 2^30);
end

end
