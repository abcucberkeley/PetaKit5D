function [job_status] = check_slurm_job_status(job_id, array_id)
% check if the job running, pending or not available
% job_status: -1: not available, 0: pending or other non-running status, 1: running
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

[~, cmdout] = system(cmd);

job_status = -1;
if contains(cmdout, job_id_str) && ~contains(cmdout, 'Invalid job id')
    if contains(cmdout, ' R ')
        job_status = 1;
    else
        job_status = 0;
    end
end

end