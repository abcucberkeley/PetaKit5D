function [tf] = stitch_in_process(numBatches, batchSize, dtype, blendMethod, imSizes, parseCluster, configFile)
% whether XR_stitching_frame_zarr_dev_v1 stitches a volume of numBatches
% batches in its own process instead of with jobs: a small volume (fewer
% than 30 batches) is, as before, if its batches (20 B per voxel of a batch
% of dtype) and, with feather blending, its distance maps (40 B per voxel of
% a tile; imSizes: n x 3 tile sizes) fit this process (master_task_fits: in
% a Slurm job, its memory). A stitching job that waits for its blocks is
% sized for what it holds (see XR_matlab_stitching_wrapper), so that it
% doesn't hold them; without the cluster (parseCluster false) always.

tf = false;
if numBatches >= 30
    return;
end
mem = prod(batchSize) * dataTypeToByteNumber(dtype) / 2^30 * 20;
if strcmpi(blendMethod, 'feather')
    mem = max(mem, max(prod(imSizes, 2)) * 4 / 2^30 * 10);
end
tf = ~parseCluster || master_task_fits(mem, 1, cluster_config_value(configFile, 'memPerCPU', 20.9), ...
    cluster_config_value(configFile, 'maxCPUNum', 24), 1);

end
