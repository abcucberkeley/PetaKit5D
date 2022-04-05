function [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str, jobLogDir] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes)
% check if a slurm computing cluster exist, if so, setup properties, if not
% set parseCluster as false
% 
% Author: Xiongtao Ruan (10/04/2020)


% initial values
parseCluster = true; 
job_log_fname = '';
job_log_error_fname = '';
slurm_constraint_str = '';

% check whether the cluster exist
[status, ~] = system('sinfo');
if status ~= 0
    warning('A slurm-based computing cluster is not exist. Set parseCluster as false.')
    parseCluster = false;
    return;
end
if parseCluster && ~exist(jobLogDir, 'dir')
    warning('The job log directory does not exist, use %s/job_logs as job log directory', dataPath)
    jobLogDir = sprintf('%s/job_logs', dataPath);
    if ~exist(jobLogDir, 'dir')
        mkdir(jobLogDir);
        fileattrib(jobLogDir, '+w', 'a');
    end
end
dir_info = dir(jobLogDir);
jobLogDir = dir_info.folder;
job_log_fname = [jobLogDir, '/job_%A_%a.out'];
job_log_error_fname = [jobLogDir, '/job_%A_%a.err'];

if cpuOnlyNodes
    slurm_constraint_str = ' --constraint=c24 ';
else
    slurm_constraint_str = '';
end


end