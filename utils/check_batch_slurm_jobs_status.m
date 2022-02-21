function [job_status] = check_batch_slurm_jobs_status(job_ids, array_ids)
% check if a batch of jobs running, pending or not available
% job_status: -1: not available, 0: pending or other non-running status, 1: running
% 
%
% Author: Xiongtao Ruan (02/07/2022)


if nargin < 2
    array_ids = [];
end

if nargin < 1
    error('job ids have to be provided!');
end

% if job ids are 0, which is the default empty value, directly return -1
job_status = -ones(numel(job_ids), 1);
valid_inds = job_ids > 0;
if ~any(valid_inds)
    return;
end

job_ids = job_ids(valid_inds);
array_ids = array_ids(valid_inds);

if isempty(array_ids)
    job_ids_s = arrayfun(@num2str, job_ids(:)', 'unif', 0);
else
    job_ids_s = arrayfun(@(x) sprintf('%d_%d', job_ids(x), array_ids(x)), 1 : numel(job_ids), 'unif', 0);
end

job_ids_str = strjoin(job_ids_s, ',');

cmd = sprintf('squeue -j %s', job_ids_str);

[~, cmdout] = system(cmd);

if contains(cmdout, 'Invalid job id')
    return;
end

cmdout_cell = strsplit(strip(cmdout), '\n');
cmdout_cell(1) = [];
tmp = regexp(cmdout_cell, '^[ ]+(\d+_?\d+)[ ]+.* (R|PD|CG|CD|F|PR|S|ST) ', 'tokens');

cnum = cellfun(@(x) numel(x), tmp);
tmp(cnum ~= 1) = [];
job_id_s_m = cellfun(@(x) x{1}{1}, tmp, 'unif', 0);
status_s_m = cellfun(@(x) x{1}{2}, tmp, 'unif', 0);

for i = 1 : numel(job_id_s_m)
    ind = strcmp(job_ids_s, job_id_s_m{i});
    if any(ind)
        switch status_s_m{i}
            case 'R'
                job_status(ind) = 1;
            case 'PD'
                job_status(ind) = 0;
        end
    end
end

end

