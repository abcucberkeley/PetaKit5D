function [nice] = coordinator_nice()
% the nice value (sbatch --nice) for a job that submits jobs of its own and
% waits for them (e.g., a stitching job and its blending jobs), so that the
% jobs it waits for go first: with Slurm's multifactor priority, a job's
% priority grows with its age (up to PriorityWeightAge) and its size (up to
% PriorityWeightJobSize), while fair share, partition and QOS are the same
% for a user's jobs of a partition and QOS. Older pending coordinating jobs
% outranked the jobs of the running ones, which then waited behind them
% while holding their resources. Lowered by the sum of the two weights, a
% coordinating job stays below the jobs it submits, whatever their ages.
% 0 without Slurm, or with another priority type.


nice = 0;
[status, out] = system('scontrol show config');
if status ~= 0 || ~contains(out, 'priority/multifactor')
    return;
end
weights = regexp(out, '^\s*PriorityWeight(Age|JobSize)\s*=\s*(\d+)', 'tokens', 'lineanchors');
for i = 1 : numel(weights)
    nice = nice + str2double(weights{i}{2});
end

end
