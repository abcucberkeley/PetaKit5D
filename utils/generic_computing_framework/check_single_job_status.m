function [job_status] = check_single_job_status(job_id, task_id)
% the status of a single job that generic_single_job_submit_wrapper submitted
% for task task_id, as check_slurm_job_status gives it (-1: not available, 0:
% pending, 1: running), where the task counts as running while a task runner
% of another job has it (see task_runner_run), and as ended (-1) once a
% runner failed it, also while its job is pending (the job is cancelled).
% When it gives -1, no job or task runner has the task, and none will take
% it: the caller may run or submit it again. A pending job submitted to one
% of several places (SlurmAlternatives) may move to another one where it can
% start now (see slurm_alternatives).

array_id = rem(task_id, 5000);
job_status = check_slurm_job_status(job_id, array_id);
if job_status == 0
    slurm_alternatives('move', job_id, array_id);
elseif job_status < 0
    slurm_alternatives('forget', job_id, array_id);
end
% (a session without pools doesn't load the runner code)
if job_status < 1 && single_job_pools()
    job_status = task_runner_run('single_status', job_id, array_id, job_status);
end

end
