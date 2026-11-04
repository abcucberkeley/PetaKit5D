function [tf] = single_job_pools(used, mark)
% whether a single job of this MATLAB session joined a pool of task runners
% (see generic_single_job_submit_wrapper); single_job_pools(true) records
% that one did. Kept apart from task_runner_run, so that the single jobs and
% status checks of a session without pools don't load it.
%
% The callers of the single jobs end the pools of a call with it:
%   mark = single_job_pools('mark')
%       at the start of the call: the number of pools this session made
%       so far (0 without pools)
%   single_job_pools('end', mark)
%       at its end (its tasks done or given up): the pools made since take
%       no more tasks, and the last of the master and their jobs to leave
%       them removes them (see task_runner_run 'single_end')

persistent poolsUsed
if nargin > 0 && ischar(used)
    tf = 0;
    if ~isempty(poolsUsed)
        switch used
            case 'mark'
                tf = task_runner_run('single_mark');
            case 'end'
                task_runner_run('single_end', mark);
        end
    end
    return;
end
if nargin > 0 && used
    poolsUsed = true;
end
tf = ~isempty(poolsUsed);

end
