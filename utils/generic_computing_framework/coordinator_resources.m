function [cpusPerTask, memAllocate] = coordinator_resources(mergeMem, computes, taskCPUs, taskMem)
% the request of a job whose task runs tasks of its own as jobs (e.g., a
% large-file task and its block jobs): what the job holds itself, MATLAB
% (2 GB) and its merges (mergeMem, GB; XR_MIP_zarr merges its MIP slabs in
% up to half of the job's memory), with one CPU; and if it computes some of
% its tasks itself (computes: masterCompute, for tasks that run on the same
% kind of node), the CPUs the caller asks for its jobs (taskCPUs, e.g.
% cpusPerTask, which the blocks get too) and at least the memory one of
% these tasks gets as a job (taskMem GB; []: the memory of its CPUs), so
% that the job is one of the workers instead of an idle job next to them.
% The master runs a task only if it fits its job (master_task_fits).

memAllocate = max(2 + mergeMem, 2 * mergeMem);
cpusPerTask = 1;
if computes
    cpusPerTask = taskCPUs;
    if ~isempty(taskMem)
        memAllocate = max(memAllocate, taskMem);
    end
end

end
