function [records] = isolation_run_tasks(tasks, varargin)
% run tasks (task function strings, see isolation_prepare) one after another
% in this MATLAB process, with a reset after each, and record what each task
% changed in the process (isolation_state before and after it).
%
% Inputs:
%   tasks : struct array with fields id, func, funcStr and threads (as
%           isolation_prepare makes them), or the .mat file that holds it.
%
% Parameters (as 'specifier'-value pairs):
%      order : vector (default: all, in order). Indices of the tasks to run, in this order.
%      reset : char or function handle (default: 'none'). Reset after each task (see isolation_reset).
%   recordFn : char (default: ''). .mat file the records are saved to after each task.
%     noSlurm : true|false (default: true). Hide Slurm's commands while the tasks run, so a task string that asks for a cluster runs locally.
%
% Output: struct array, one per task run: id, func, ok, message, time (s),
%   resetTime (s), resetActions, changes (what the task changed, see
%   isolation_state_diff), leftover (what differs from the state before the
%   first task after the reset), and memory (RSS, peak, GPU) after the task
%   and after the reset.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tasks', @(x) isstruct(x) || ischar(x));
ip.addParameter('order', [], @isnumeric);
ip.addParameter('reset', 'none', @(x) ischar(x) || isa(x, 'function_handle'));
ip.addParameter('recordFn', '', @ischar);
ip.addParameter('noSlurm', true, @islogical);
ip.parse(tasks, varargin{:});

pr = ip.Results;
if ischar(tasks)
    a = load(tasks, 'tasks');
    tasks = a.tasks;
end
order = pr.order;
if isempty(order)
    order = 1 : numel(tasks);
end
resetLevel = pr.reset;

origPATH = getenv('PATH');
if pr.noSlurm && ~ispc
    stubDir = hide_slurm();
    cleanupObj = onCleanup(@() restore_slurm(stubDir, origPATH));
end

% a new process has no last warning
lastwarn('', '');
st0 = isolation_state();
records = struct('id', {}, 'func', {}, 'ok', {}, 'message', {}, 'time', {}, 'resetTime', {}, ...
    'resetActions', {}, 'changes', {}, 'leftover', {}, 'rssGB', {}, 'peakGB', {}, 'gpuMemMB', {}, ...
    'rssAfterResetGB', {}, 'gpuMemAfterResetMB', {});
for k = 1 : numel(order)
    task = tasks(order(k));
    fprintf('\n=== isolation task %d / %d: %s (%s)\n', k, numel(order), task.id, task.func);
    before = isolation_state();
    if isfield(task, 'threads') && ~isempty(task.threads)
        maxNumCompThreads(task.threads);
    end
    ok = true;
    message = '';
    t0 = tic;
    try
        feval(str2func(['@()', task.funcStr]));
    catch ME
        ok = false;
        message = getReport(ME, 'extended', 'hyperlinks', 'off');
        fprintf(2, 'Task %s failed:\n%s\n', task.id, message);
    end
    t = toc(t0);
    after = isolation_state();
    resetInfo = isolation_reset(resetLevel, st0);
    afterReset = isolation_state();

    r.id = task.id;
    r.func = task.func;
    r.ok = ok;
    r.message = message;
    r.time = t;
    r.resetTime = resetInfo.time;
    r.resetActions = resetInfo.actions;
    r.changes = isolation_state_diff(before, after);
    r.leftover = isolation_state_diff(st0, afterReset);
    r.rssGB = after.rssGB;
    r.peakGB = after.peakGB;
    r.gpuMemMB = after.gpuMemMB;
    r.rssAfterResetGB = afterReset.rssGB;
    r.gpuMemAfterResetMB = afterReset.gpuMemMB;
    records(end + 1) = r; %#ok<AGROW>
    fprintf('=== isolation task %s: %s in %.2f s, reset %.3f s; changed: %s\n', task.id, ...
        string(ok), t, resetInfo.time, strjoin(r.changes, '; '));
    if ~isempty(pr.recordFn)
        save(pr.recordFn, 'records', 'order');
    end
end

end


function [stubDir] = hide_slurm()
% failing sinfo, sbatch, squeue and scancel first on the PATH

stubDir = strrep(tempname, '\', '/');
mkdir(stubDir);
for cmd = {'sinfo', 'sbatch', 'squeue', 'scancel'}
    fn = [stubDir, '/', cmd{1}];
    fid = fopen(fn, 'w');
    fprintf(fid, '#!/bin/sh\nexit 1\n');
    fclose(fid);
    fileattrib(fn, '+x');
end
setenv('PATH', [stubDir, pathsep, getenv('PATH')]);

end


function [] = restore_slurm(stubDir, origPATH)

setenv('PATH', origPATH);
if exist(stubDir, 'dir')
    rmdir(stubDir, 's');
end

end
