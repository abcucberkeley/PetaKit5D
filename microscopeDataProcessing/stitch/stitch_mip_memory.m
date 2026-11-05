function [mipMem] = stitch_mip_memory(clusterMode, configFile)
% the memory (GB) that XR_stitching_frame_zarr_dev_v1 measures the stitched
% volume against for its MIPs (the z MIP from the blending tasks for volumes
% under half of it, else x, y and z MIPs from cluster jobs): this process's
% memory, and for a stitching job of the cluster (clusterMode in a Slurm
% job) at least the largest job of the config (maxCPUNum x memPerCPU, else
% the job framework's defaults). Stitching jobs are sized for what they
% hold (see XR_matlab_stitching_wrapper); their MIPs are those of a job
% with memory for all its tiles, which stitching jobs got before (at most
% the largest job).

mipMem = getSystemMemory();
if clusterMode && ~isempty(getenv('SLURM_JOB_ID'))
    mipMem = max(mipMem, cluster_config_value(configFile, 'maxCPUNum', 24) * ...
        cluster_config_value(configFile, 'memPerCPU', 20.9));
end

end
