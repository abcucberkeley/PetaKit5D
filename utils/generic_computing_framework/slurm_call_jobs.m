function [varargout] = slurm_call_jobs(action, varargin)
% the Slurm jobs of a job framework call (slurm_cluster_ and
% mcc_slurm_cluster_generic_computing_wrapper). Each job has the call's tag
% and its batch as its comment (<tag>_<batch>), so that the job of a
% submission without an answer (a timeout, or a lost connection to the
% controller, after which Slurm may have made the job) is looked for with
% squeue before the batch is submitted again; such a batch was submitted
% again at once, and the first job stayed in the queue, unknown to the call.
% A call that ends with an error or is stopped (Ctrl+C) cancels its jobs,
% which stayed in the queue before.
%
% Usage:
%   [reg, cleanupObj] = slurm_call_jobs('start', uuid)
%       the record of a call with uuid (a handle: its tag, pk5d_<the first
%       12 letters and digits of uuid>, its jobs and its batches submitted
%       without an answer) and an onCleanup that cancels the jobs (scancel of
%       their ids) unless the call ended normally ('done')
%   cmd = slurm_call_jobs('tag', reg, cmd, b)
%       sbatch command cmd of batch b with the comment (after its --array)
%   [answer, jobId] = slurm_call_jobs('answer', reg, cmdout, b)
%       sbatch's answer cmdout for batch b: 'submitted' (the job is
%       recorded), 'unknown' (no answer from Slurm; the batch is recorded) or
%       'refused'
%   [job_ids, job_status_mat, trial_counter, unknown_b] = slurm_call_jobs('resolve', reg, unknown_b, job_ids, job_status_mat, trial_counter, taskBatchNum)
%       the batches submitted without an answer (unknown_b) whose jobs squeue
%       lists get them (as at a submission), and those it doesn't list 30 s
%       after the submission can be submitted again; none if squeue fails
%   slurm_call_jobs('done', reg)
%       at the call's normal end (it cancels its jobs itself): cancels the
%       jobs squeue lists of the batches still without an answer


switch action
    case 'start'
        uuid = regexprep(varargin{1}, '[^0-9A-Za-z]', '');
        if isempty(uuid)
            uuid = regexprep(get_uuid(), '[^0-9A-Za-z]', '');
        end
        reg = containers.Map('KeyType', 'char', 'ValueType', 'any');
        reg('tag') = ['pk5d_', uuid(1 : min(12, end))];
        reg('ids') = zeros(1, 0);
        % batches without an answer: [batch, time of the submission]
        reg('unknown') = zeros(0, 2);
        reg('done') = false;
        varargout = {reg, onCleanup(@() cancel_jobs(reg))};
    case 'tag'
        [reg, cmd, b] = varargin{1 : 3};
        varargout{1} = regexprep(cmd, '^(sbatch --array=\d+) ', sprintf('$1 --comment=%s_%d ', reg('tag'), b), 'once');
    case 'answer'
        [reg, cmdout, b] = varargin{1 : 3};
        jobId = [];
        % sbatch's errors after which Slurm may have made the job
        lost = {'Socket timed out', 'Unable to contact slurm controller', 'Zero Bytes were transmitted or received', ...
            'Connection timed out', 'Connection refused', 'Transport endpoint is not connected', ...
            'Communication connection failure'};
        tok = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens', 'once');
        if ~isempty(tok)
            answer = 'submitted';
            jobId = str2double(tok{1});
            reg('ids') = [reg('ids'), jobId];
        elseif contains(cmdout, lost)
            answer = 'unknown';
            reg('unknown') = [reg('unknown'); b, now_s()];
            fprintf('Task % 4d:    no answer from Slurm, its job is looked for before it is submitted again\n', b);
        else
            answer = 'refused';
        end
        varargout = {answer, jobId};
    case 'resolve'
        [reg, unknown_b, job_ids, job_status_mat, trial_counter, taskBatchNum] = varargin{1 : 6};
        u = reg('unknown');
        [ids, states, ok] = find_jobs(reg, u(:, 1));
        if ok
            nF = numel(job_ids);
            keep = false(size(u, 1), 1);
            for i = 1 : size(u, 1)
                b = u(i, 1);
                fs = (b - 1) * taskBatchNum + 1 : min(b * taskBatchNum, nF);
                if ids(i) > 0
                    fprintf('Task % 4d:    its job %d was submitted\n', b, ids(i));
                    job_ids(fs) = ids(i);
                    job_status_mat(fs, 1) = states(i);
                    trial_counter(fs) = trial_counter(fs) + 1;
                    reg('ids') = [reg('ids'), ids(i)];
                else
                    % (Slurm may make a job after sbatch gave up on its answer)
                    keep(i) = now_s() - u(i, 2) < 30;
                end
            end
            unknown_b(u(~keep, 1)) = false;
            reg('unknown') = u(keep, :);
        end
        varargout = {job_ids, job_status_mat, trial_counter, unknown_b};
    case 'done'
        reg = varargin{1};
        u = reg('unknown');
        if ~isempty(u)
            [ids, ~, ok] = find_jobs(reg, u(:, 1));
            cancel_ids(ids);
            if ~ok
                warning('slurm_call_jobs:unknownJobs', ['Jobs of %d task batches whose submissions got no answer ', ...
                    'may be in the queue (comments %s_<batch>).'], size(u, 1), reg('tag'));
            end
        end
        reg('done') = true;
    otherwise
        error('slurm_call_jobs:action', 'Unknown action %s.', action);
end

end


function [jobIds, states, ok] = find_jobs(reg, bs)
% the jobs (-1: none) and their states (1 running, 0 pending, -1 other) that
% squeue lists for batches bs of the call (by their comments); ok is false
% if squeue failed

jobIds = -ones(numel(bs), 1);
states = -ones(numel(bs), 1);
[status, cmdout] = system('squeue -u "$(id -un)" -h -r -o "%i %t %k"');
ok = status == 0;
if ~ok
    return;
end
tok = regexp(cmdout, sprintf('^(\\d+)\\S* (\\S+) %s_(\\d+)[ \\t]*$', reg('tag')), 'tokens', 'lineanchors');
for i = 1 : numel(tok)
    k = find(bs == str2double(tok{i}{3}), 1);
    if isempty(k)
        continue;
    end
    jobIds(k) = str2double(tok{i}{1});
    states(k) = double(strcmp(tok{i}{2}, 'R')) - double(~any(strcmp(tok{i}{2}, {'R', 'PD'})));
end

end


function [] = cancel_jobs(reg)
% (the onCleanup of 'start') a call that ended with an error or was stopped
% cancels the jobs it submitted and those squeue lists of its batches
% submitted without an answer

if reg('done')
    return;
end
try
    ids = reg('ids');
    u = reg('unknown');
    if ~isempty(u)
        found = find_jobs(reg, u(:, 1));
        ids = [ids, found(:)'];
    end
    ids = unique(ids(ids > 0));
    if ~isempty(ids)
        fprintf('The call ended before its tasks were done: cancel its jobs\n');
        cancel_ids(ids);
    end
catch ME
    warning('slurm_call_jobs:cancel', 'The jobs of the call could not be cancelled: %s', ME.message);
end

end


function [] = cancel_ids(ids)

ids = unique(ids(ids > 0));
if ~isempty(ids)
    system(sprintf('scancel %s', num2str(ids(:)')), '-echo');
end

end


function [t] = now_s()
% the time in seconds

t = posixtime(datetime('now'));

end
