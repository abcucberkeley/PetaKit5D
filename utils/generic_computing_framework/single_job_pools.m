function [tf] = single_job_pools(used)
% whether a single job of this MATLAB session joined a pool of task runners
% (see generic_single_job_submit_wrapper); single_job_pools(true) records
% that one did. Kept apart from task_runner_run, so that the single jobs and
% status checks of a session without pools don't load it

persistent poolsUsed
if nargin > 0 && used
    poolsUsed = true;
end
tf = ~isempty(poolsUsed);

end
