function [out] = task_runner(action, varargin)
% runs several job framework tasks in one process (a task runner): after a
% task, the process claims the next task (first of its own batch, then
% unstarted tasks of the framework call's other batches), resets itself to
% the state it had before its first task and runs it. The job script
% (task_runner.sh, written for each call by task_runner_run) claims the
% tasks and starts the processes. A process ends when no task is left,
% after a task failed (the job script runs that task again in a process of
% its own), when it holds more memory, open files or threads than it had at
% its start, or after 10 minutes of tasks (the job script starts a new
% process for the next task).
%
% The reset gives back the current folder, MATLAB path, warning states (and
% an empty last warning), the random number generator's default state (as
% at the start of MATLAB) and the environment variables; it removes the
% variables tasks made in the base workspace and as globals, closes the
% files they opened and the figures they left, shuts down a parallel pool,
% clears functions and their persistent variables (clear functions; mex
% files that locked themselves stay loaded) and deselects the GPU (freeing
% its memory) if a task used it. A task that changes the number of
% computational threads ends the process: results depend on the thread
% count, and MATLAB fixes at a process's first FFT whether its FFTs split
% lines among threads, so all tasks of a process run with the thread count
% of its start (which is that of the processes the tasks would have had
% without runners: the job script takes only tasks of jobs with the same
% CPUs and GPU, see task_runner_run).
%
% Usage:
%   r = task_runner('start')
%       in mccMaster before its task, if the job script set PETAKIT5D_RUNNER:
%       the runner's settings and the state of the new process ([] if the
%       process is no runner)
%   task_runner('continue', r, taskSec)
%       in mccMaster after its task (which took taskSec seconds): runs the
%       next tasks
%   task_runner('matlab', D, own, steal, st, L, f)
%       a MATLAB job's runner (the job script starts it after setup([])):
%       runs task f and then the next tasks, with the call folder D, its
%       batch or task list own, stealing (0, 1, 2), its status file st and
%       local folder L
%   v = task_runner('version')
%       the version of the runner protocol (compile_mccMaster writes it with
%       a build; the framework uses runners only with builds of this version)
%
% The status file of a process (node-local) has lines 'c <task>' (claimed),
% 'd <task>' (done), 's <task>' (skipped, its output exists, or given back),
% 'e <task>' (failed), 't <seconds>' (a task's time) and 'x <reason>' at the
% end (none: no task left; recycle; failed). PETAKIT5D_RUNNER_MEM_GB sets the
% memory margin (default 0.5 GB), PETAKIT5D_RUNNER_MINUTES how long a
% process runs tasks (default 10; see next_task).


out = [];
switch action
    case 'start'
        out = start_runner();
    case 'continue'
        run_mcc(varargin{:});
    case 'matlab'
        run_matlab(varargin{:});
    case 'version'
        out = 1;
    otherwise
        error('task_runner:unknownAction', 'Unknown action %s.', action);
end

end


function [r] = start_runner()
% the runner's settings (PETAKIT5D_RUNNER: <call folder>|<own>|<steal>|<slot>|
% <status file>|<local folder>|<task>) and the state of this new process

r = [];
cfg = getenv('PETAKIT5D_RUNNER');
if isempty(cfg)
    return;
end
% the tasks see the environment a process of their own would have
try
    unsetenv('PETAKIT5D_RUNNER');
catch
    setenv('PETAKIT5D_RUNNER', '');
end
c = strsplit(cfg, '|');
if numel(c) < 7
    return;
end
r.D = c{1};
r.own = c{2};
r.steal = str2double(c{3});
r.slot = str2double(c{4});
r.st = c{5};
r.L = c{6};
r.task = str2double(c{7});
r.t0 = tic;
try
    r.state = capture_state(false);
catch ME
    % (a MATLAB version without what the reset needs) the process runs
    % only its task, and the job script a new process for the next one
    fprintf('Task runner %s: this process runs only task %d (%s)\n', runner_name(r), r.task, ME.message);
    r.state = [];
end

end


function [s] = capture_state(keepBase)
% the state of this process that tasks can change; keepBase: the base
% workspace has variables to keep (a compiled master's is empty)

s.cwd = cd;
s.path = matlabpath;
s.warning = warning;
s.env = getenv();
s.nThreads = maxNumCompThreads;
s.fids = open_files();
s.globals = who('global');
s.base = {};
if keepBase
    s.base = evalin('base', 'who');
end
s.figs = figures();
s.rssKB = rss_kb();

end


function [] = run_mcc(r, taskSec)
% compiled runner: the job script started this process with the task line
% of task r.task, which ran in mccMaster; then the next tasks

append_status(r.st, sprintf('d %d\nt %.3f', r.task, taskSec));
if isempty(r.state)
    append_status(r.st, 'x recycle');
    return;
end
fprintf('Task runner %s: task %d done in %.3f s, memory %.2f GB (%.2f GB before it)\n', runner_name(r), r.task, ...
    taskSec, rss_kb() / 1048576, r.state.rssKB / 1048576);
est = taskSec;
dumpFn = sprintf('%s/dump_%d', r.L, r.slot);
while true
    if threads_changed(r)
        return;
    end
    t0 = tic;
    [f, r] = next_task(r, est);
    if isempty(f)
        return;
    end
    [envs, args] = read_dump(dumpFn);
    if isempty(args)
        % its arguments couldn't be read: the claim goes back, and the job
        % script runs the task with its line in a new process
        fprintf('Task runner %s: no arguments for task %d, a new process runs it\n', runner_name(r), f);
        delete(sprintf('%s/claims/t%d', r.D, f));
        append_status(r.st, sprintf('s %d\nx recycle', f));
        return;
    end
    % the environment its task line gives the task
    names = {'PETAKIT5D_TASK_METRICS', 'PETAKIT5D_TASK', 'PETAKIT5D_TASK_GPU'};
    for i = 1 : numel(names)
        if ~isempty(envs{i}) || ~isempty(getenv(names{i}))
            setenv(names{i}, envs{i});
        end
    end
    fprintf('Task runner %s: task %d (claimed, process reset in %.3f s)\n', runner_name(r), f, toc(t0));
    t1 = tic;
    try
        mccMaster(args{:});
        ok = true;
    catch ME
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        ok = false;
    end
    t = toc(t1);
    est = max(est, t);
    if ~ok
        fprintf('Task runner %s: task %d failed after %.3f s, a process of its own runs it again\n', runner_name(r), f, t);
        append_status(r.st, sprintf('e %d\nx failed', f));
        return;
    end
    append_status(r.st, sprintf('d %d\nt %.3f', f, t));
    fprintf('Task runner %s: task %d done in %.3f s, memory %.2f GB\n', runner_name(r), f, t, rss_kb() / 1048576);
end

end


function [] = run_matlab(D, own, steal, st, L, f)
% MATLAB job runner: runs task f (claimed by the job script) and the next
% tasks of the call (task strings in D/tasks.mat)

r.D = D;
r.own = own;
r.steal = steal;
r.slot = 1;
r.st = st;
r.L = L;
r.t0 = tic;
r.state = [];
try
    r.state = capture_state(true);
catch ME
    % (see start_runner) only this task
    fprintf('Task runner %s: this process runs only task %d (%s)\n', runner_name(r), f, ME.message);
end
a = load([D, '/tasks.mat'], 'funcStrs', 'runDir', 'GPUJob');
est = 0;
while true
    funcStr = a.funcStrs{f};
    fprintf('Task runner %s: task %d: %s\n', runner_name(r), f, funcStr);
    tm = task_metrics('start', funcStr, a.runDir, f, 'slurm_matlab', a.GPUJob);
    t1 = tic;
    try
        run_task_string(funcStr);
        task_metrics('stop', tm);
        ok = true;
    catch ME
        task_metrics('stop', tm, ME);
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        ok = false;
    end
    t = toc(t1);
    est = max(est, t);
    if ~ok
        fprintf('Task runner %s: task %d failed after %.3f s\n', runner_name(r), f, t);
        append_status(st, sprintf('e %d\nx failed', f));
        return;
    end
    append_status(st, sprintf('d %d\nt %.3f', f, t));
    fprintf('Task runner %s: task %d done in %.3f s, memory %.2f GB\n', runner_name(r), f, t, rss_kb() / 1048576);
    if isempty(r.state)
        append_status(st, 'x recycle');
        return;
    end
    if threads_changed(r)
        return;
    end
    t0 = tic;
    [f, r] = next_task(r, est);
    if isempty(f)
        return;
    end
    fprintf('Task runner %s: claimed task %d, process reset in %.3f s\n', runner_name(r), f, toc(t0));
end

end


function [f, r] = next_task(r, est)
% claims the next task and resets the process for it; [] if there is none,
% or if the process holds more memory than a process that ran a task of
% the call holds (it then ends: recycle). That memory (warm) is what the
% process held after its first task and the reset (code, libraries such as
% the GPU's), at most twice the margin more than at its start; later tasks
% may add up to the margin to it. A process also ends after it ran tasks for
% PETAKIT5D_RUNNER_MINUTES (default 10): whatever builds up in a process
% (tasks got up to 8% slower over 100 tasks of a demo tile) stays bounded,
% for one MATLAB Runtime start per 10 minutes

f = [];
margin = str2double(getenv('PETAKIT5D_RUNNER_MEM_GB'));
if isnan(margin)
    margin = 0.5;
end
maxMinutes = str2double(getenv('PETAKIT5D_RUNNER_MINUTES'));
if isnan(maxMinutes)
    maxMinutes = 10;
end
if toc(r.t0) >= maxMinutes * 60
    fprintf('Task runner %s: tasks for %.1f min, a new process runs the next task\n', runner_name(r), toc(r.t0) / 60);
    append_status(r.st, 'x recycle');
    return;
end
nextFn = sprintf('%s/next_%d', r.L, r.slot);
cmd = sprintf('PK_PID=%d bash "%s/run_tasks.sh" claim "%s" %d %d %d "%s" "%s" "%s"', feature('getpid'), ...
    r.D, r.own, r.steal, r.slot, ceil(est), r.st, r.L, nextFn);
status = system(cmd);
txt = '';
if status == 0 && isfile(nextFn)
    txt = strtrim(fileread(nextFn));
end
if startsWith(txt, 'recycle')
    fprintf('Task runner %s: %s, a new process runs the next task\n', runner_name(r), txt(9 : end));
    append_status(r.st, 'x recycle');
    return;
end
f = sscanf(txt, '%d', 1);
if isempty(f)
    append_status(r.st, 'x none');
    return;
end
append_status(r.st, sprintf('c %d', f));
reset_state(r.state);
% the memory the process holds without the task's variables and caches
rss = rss_kb();
if ~isfield(r, 'warmKB')
    r.warmKB = rss;
    limit = r.state.rssKB + 2 * margin * 1048576;
else
    limit = r.warmKB + margin * 1048576;
end
if rss > limit
    fprintf('Task runner %s: memory %.2f GB (%.2f GB at its start), a new process runs task %d\n', ...
        runner_name(r), rss / 1048576, r.state.rssKB / 1048576, f);
    % the claim goes back for a new process (of this slot or another)
    delete(sprintf('%s/claims/t%d', r.D, f));
    append_status(r.st, sprintf('s %d\nx recycle', f));
    f = [];
end

end


function [tf] = threads_changed(r)
% whether the task changed the number of computational threads (then the
% process ends, see the header)

tf = maxNumCompThreads ~= r.state.nThreads;
if tf
    fprintf('Task runner %s: the task changed the number of threads (%d, %d at the start), a new process runs the next task\n', ...
        runner_name(r), maxNumCompThreads, r.state.nThreads);
    append_status(r.st, 'x recycle');
end

end


function [] = run_task_string(funcStr__)
% a task string in a workspace of its own (no variables of the runner)

try
    task__ = str2func(['@()', funcStr__]);
catch
    eval(funcStr__);
    return;
end
feval(task__);

end


function [envs, args] = read_dump(fn)
% the environment (PETAKIT5D_TASK_METRICS, PETAKIT5D_TASK,
% PETAKIT5D_TASK_GPU) and the arguments the task line gives the compiled
% master, as the stand-in master wrote them (NUL separated)

envs = {};
args = {};
fid = fopen(fn, 'r');
if fid < 0
    return;
end
raw = fread(fid, inf, 'uint8=>uint8')';
fclose(fid);
if isempty(raw) || raw(end) ~= 0
    return;
end
parts = strsplit(native2unicode(raw(1 : end - 1), 'UTF-8'), char(0), 'CollapseDelimiters', false);
if numel(parts) < 4
    return;
end
envs = parts(1 : 3);
args = parts(4 : end);

end


function [] = reset_state(s)
% back to the state of the process before its first task (see the header)

try
    cd(s.cwd);
catch
end
if ~strcmp(matlabpath, s.path)
    path(s.path);
end
warning(s.warning);
lastwarn('', '');
rng('default');
env = getenv();
if ~isequal(env, s.env)
    restore_env(env, s.env);
end
fids = open_files();
for i = 1 : numel(fids)
    if ~any(fids(i) == s.fids)
        fclose(fids(i));
    end
end
if isempty(s.base)
    evalin('base', 'clear');
else
    names = evalin('base', 'who');
    names = names(~ismember(names, s.base));
    if ~isempty(names)
        evalin('base', ['clear ', strjoin(names, ' ')]);
    end
end
names = who('global');
if ~isempty(names)
    names = names(~ismember(names, s.globals));
    if ~isempty(names)
        clear('global', names{:});
    end
end
% a parallel pool a task started (if one could have: its code is loaded)
if any(strcmp(inmem, 'gcp'))
    try
        delete(gcp('nocreate'));
    catch
    end
end
clear functions;
figs = figures();
for i = 1 : numel(figs)
    if ~any(figs(i) == s.figs)
        close(figs(i), 'force');
    end
end
reset_gpu();

end


function [] = restore_env(env, env0)
% the environment variables of env0

names = keys(env);
for i = 1 : numel(names)
    if ~isKey(env0, names(i))
        unsetenv(names(i));
    end
end
names = keys(env0);
for i = 1 : numel(names)
    if ~isKey(env, names(i)) || ~strcmp(env(names(i)), env0(names(i)))
        setenv(names(i), env0(names(i)));
    end
end

end


function [fids] = open_files()
% identifiers of the files opened with fopen

try
    fids = openedFiles();
catch
    fids = fopen('all');
end

end


function [figs] = figures()
% the open figures (none without Java)

figs = [];
if usejava('jvm')
    try
        figs = groot().Children;
    catch
    end
end

end


function [] = reset_gpu()
% frees the GPU memory and device state of the tasks (no device selected,
% as in a new process), if a task used the GPU

try
    if parallel.internal.gpu.isAnyDeviceSelected
        gpuDevice([]);
    end
catch
end

end


function [kb] = rss_kb()
% resident memory of this process in kB (Linux; NaN elsewhere)

kb = NaN;
fid = fopen('/proc/self/status', 'r');
if fid < 0
    return;
end
txt = fread(fid, [1, inf], '*char');
fclose(fid);
k = strfind(txt, 'VmRSS:');
if ~isempty(k)
    kb = sscanf(txt(k(1) + 6 : end), '%f', 1);
end

end


function [name] = runner_name(r)

name = sprintf('%s_%s_%d', getenv('SLURM_JOB_ID'), getenv('SLURM_ARRAY_TASK_ID'), r.slot);

end


function [] = append_status(fn, s)
% a line of the process's status file (node-local, this process only)

fid = fopen(fn, 'a');
if fid < 0
    return;
end
fprintf(fid, '%s\n', s);
fclose(fid);

end
