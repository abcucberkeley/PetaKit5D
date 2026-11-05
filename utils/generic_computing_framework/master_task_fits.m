function [tf, taskMem, masterMem] = master_task_fits(memAllocate, cpusPerTask, memPerCPU, maxCPUNum, minCPUNum, masterMem)
% whether the master of a job framework call in cluster mode may run a task
% itself (masterCompute): outside a Slurm job always (the user's MATLAB has
% all of the machine, as before); in a Slurm job (a job that coordinates
% other jobs, or a pipeline's job) only if the task's memory fits the job's
% memory (getSystemMemory), so that the job isn't killed for memory with a
% task that its own job would have run. The task's memory is memAllocate (at
% most what its job gets, maxCPUNum CPUs), or without it the memory of its
% job's CPUs. memPerCPU is often a rounded value of the jobs' memory per CPU,
% so up to 5% over the job's memory counts as fitting; an unknown memory
% doesn't fit.
%
% Inputs (as the job framework has them; memAllocate and memPerCPU in GB):
%   memAllocate : empty or a number. Memory estimate of the task.
%   cpusPerTask : a number. CPUs of the task's job.
%     memPerCPU : a number. Memory per CPU of the jobs.
%     maxCPUNum : a number (optional, default: inf). Most CPUs of a job.
%     minCPUNum : a number (optional, default: 1). Fewest CPUs of a job.
%     masterMem : empty or a number (optional, default: []). The job's memory if
%                 known (getSystemMemory, read once by callers that check many
%                 tasks; in a job started with --mem-per-cpu it asks scontrol).
%
% Outputs:
%            tf : true if the master may run the task.
%       taskMem : the task's memory in GB.
%     masterMem : the master's memory in GB (getSystemMemory; inf outside Slurm).


if nargin < 4 || isempty(maxCPUNum)
    maxCPUNum = inf;
end
if nargin < 5 || isempty(minCPUNum)
    minCPUNum = 1;
end

if isempty(memAllocate)
    taskMem = min(maxCPUNum, max(minCPUNum, cpusPerTask)) * memPerCPU;
else
    taskMem = min(memAllocate, maxCPUNum * memPerCPU);
end

tf = true;
if isempty(getenv('SLURM_JOB_ID'))
    masterMem = inf;
    return;
end
if nargin < 6 || isempty(masterMem)
    masterMem = getSystemMemory();
end
tf = taskMem <= masterMem * 1.05;

end
