function [job_status, is_unknown] = check_slurm_job_status(job_id, array_id)
% check if the job running, pending or not available
% job_status: -1: not available, 0: pending or other non-running status, 1: running
% is_unknown: true if squeue failed other than with "Invalid job id" (e.g.,
%   a timeout of a busy controller); the job may still run, and job_status
%   is 1 (it was -1, so that callers submitted the task again)
%
%
% Author: Xiongtao Ruan (02/18/2020)


if nargin < 2
    array_id = [];
end

if nargin < 1
    error('a job id has to be provided!');
end

if isempty(array_id)
    job_id_str = num2str(job_id);
else
    job_id_str = sprintf('%d_%d', job_id, array_id);
end

cmd = sprintf('squeue -j %s', job_id_str);

[status, cmdout] = system(cmd);

job_status = -1;
is_unknown = status ~= 0 && ~contains(cmdout, 'Invalid job id');
if is_unknown
    job_status = 1;
    return;
end
if contains(cmdout, job_id_str) && ~contains(cmdout, 'Invalid job id')
    if contains(cmdout, ' R ')
        job_status = 1;
    else
        job_status = 0;
    end
end

end
