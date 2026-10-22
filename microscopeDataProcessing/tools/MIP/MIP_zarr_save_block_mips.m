function [] = MIP_zarr_save_block_mips(flagFullname, blockBboxes, mips, mipBatchSize, t)
% save the flag of a block task of large-file processing with the MIPs of
% its blocks (MIP_zarr_block_mips), which MIP_zarr_merge_block_mips merges
% into the MIP slabs. The flag marks the task as done, so it is written under
% a temporary name first; uncompressed, which is much faster (-v6 holds up
% to 2 GB).
%
% flagFullname: the task's flag file
% blockBboxes: the regions of the blocks in the zarr, one row per block
% mips: the y, x and z MIPs of the blocks, one row per block
% mipBatchSize: batch size of the MIP slabs
% t: elapsed time of the task (optional)


if nargin < 5
    t = [];
end

yMIPs = mips(:, 1);
xMIPs = mips(:, 2);
zMIPs = mips(:, 3);

flagTmpname = sprintf('%s_%s.mat', flagFullname(1 : end - 4), get_uuid());
mipInfo = whos('mips');
if mipInfo.bytes < 2^31
    save(flagTmpname, 't', 'blockBboxes', 'mipBatchSize', 'yMIPs', 'xMIPs', 'zMIPs', '-v6');
else
    save(flagTmpname, 't', 'blockBboxes', 'mipBatchSize', 'yMIPs', 'xMIPs', 'zMIPs', '-v7.3', '-nocompression');
end
movefile(flagTmpname, flagFullname);

end
