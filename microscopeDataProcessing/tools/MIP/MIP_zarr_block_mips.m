function [mips] = MIP_zarr_block_mips(block, bbox, mipBatchSize, dtype)
% y, x and z MIPs of a block written to region bbox of a zarr, for the MIP
% slabs of XR_MIP_zarr with batches of mipBatchSize (MIP_zarr_batch_size):
% along each axis j, the max over each part of the block within one batch.
% mips{j} has the block's size with one element per batch the block overlaps
% along j. They are cast to the zarr's data type dtype as writezarr casts
% the block (rounded and saturated, which keeps the order of values), so
% they are the MIPs of the block as written. MIP_zarr_merge_block_mips
% merges them into the slabs.
%
% block: the data written to the zarr
% bbox: the region of the block in the zarr, [ymin, xmin, zmin, ymax, xmax, zmax]
% mipBatchSize: batch size of the MIP slabs
% dtype: data type of the zarr


% range_max_3d_mex is the same as max, 2-3 times faster for uint16 (max if
% it isn't built for the platform)
useMex = true;
mips = cell(1, 3);
for j = 1 : 3
    % local starts of the parts in each batch along j
    bStart = floor((bbox(j) - 1) / mipBatchSize(j)) * mipBatchSize(j) + 1;
    starts = [bbox(j), bStart + mipBatchSize(j) : mipBatchSize(j) : bbox(j + 3)] - bbox(j) + 1;

    if useMex
        try
            mips{j} = range_max_3d_mex(block, j, starts);
        catch ME
            disp(ME);
            useMex = false;
        end
    end
    if ~useMex
        if isscalar(starts)
            mips{j} = max(block, [], j);
        else
            % one max per part (copies each part)
            sz = size(block, 1 : 3);
            ends = [starts(2 : end) - 1, sz(j)];
            idx = {':', ':', ':'};
            mipCell = cell(numel(starts), 1);
            for k = 1 : numel(starts)
                idx{j} = starts(k) : ends(k);
                mipCell{k} = max(block(idx{:}), [], j);
            end
            mips{j} = cat(j, mipCell{:});
        end
    end
    if ~isa(mips{j}, dtype)
        mips{j} = cast(mips{j}, dtype);
    end
end

end
