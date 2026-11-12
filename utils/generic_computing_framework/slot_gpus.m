function [gpus] = slot_gpus(n)
% the GPUs of up to n slots (or workers) of this process, one each, as
% CUDA_VISIBLE_DEVICES values: the first n of its own CUDA_VISIBLE_DEVICES
% (a Slurm job's GPUs; GPU k - 1 for slot k was another job's where Slurm
% doesn't renumber the GPUs of a job), else GPUs 0 to n - 1 of those MATLAB
% can use (usable_gpu_num); fewer if there are fewer, none if none is found
% (then the slots get no GPU of their own, as before). The same rule as the
% GPU slots of compiled jobs (task_runner.sh, the GNU parallel jobs of
% mcc_slurm_cluster_generic_computing_wrapper).

gpus = cell(1, 0);
if n < 1
    return;
end
list = strtrim(strsplit(getenv('CUDA_VISIBLE_DEVICES'), ','));
list = list(~cellfun(@isempty, list));
if ~isempty(list)
    gpus = list(1 : min(n, end));
else
    gpus = arrayfun(@(x) sprintf('%d', x), 0 : min(n, usable_gpu_num()) - 1, 'unif', 0);
end

end
