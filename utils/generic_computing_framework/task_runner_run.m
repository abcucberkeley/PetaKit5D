function [varargout] = task_runner_run(action, varargin)
% sets up and manages the task runners of a job framework call (see
% task_runner): the files of the call folder for the job script
% (run_tasks.sh, from task_runner.sh), the job commands, and the claims the
% master reads before it submits or resubmits a batch.
%
% Usage:
%   tf = task_runner_run('supported', MCCMasterStr)
%       whether the compiled master has task runners of this code's version
%       (task_runner.txt next to run_mccMaster.sh, written by compile_mccMaster)
%   rc = task_runner_run('setup', D, tlineStrs, MCCMasterStr, taskBatchNum, tmpInfo, excluded, keys)
%       compiled mode: the runner files in the call folder D for the task
%       lines tlineStrs (one per task, '' for tasks done at the start);
%       tmpInfo (local runs, else []) has the tasks' .tmp files (tmpFns)
%       and outputs (outFns) and the wait time in minutes (unitWaitTime);
%       excluded tasks (missing inputs) run only in their own batch; keys
%       (one per batch, see 'keys') say which batches' tasks a runner may
%       take: those with the key of its own batch. rc.on is false if
%       runners can't be used (the framework then runs as before)
%   rc = task_runner_run('setup_matlab', D, funcStrs, outputFullpaths, finalOutFullpath, isDone, taskBatchNum, MatlabLaunchStr, runDir, GPUJob, keys)
%       MATLAB mode: the same for MATLAB jobs
%   keys = task_runner_run('keys', nB, cpus, GPUJob, masterBatch)
%       the batches' keys: the CPUs (cpus, the threads of the processes)
%       and GPU of the jobs that run the batches without runners, and
%       'master' for the batch the master runs itself (masterBatch, 0 if
%       none). A process's results depend on its thread count, and on the
%       thread count of its first FFT (MATLAB fixes then whether its FFTs
%       split lines among threads), so a runner only takes tasks whose own
%       process would have had the same
%   cmd = task_runner_run('job', rc, own, slots, gpu, steal, retries)
%       the job script command of a compiled job; own is a batch number or
%       a list of task indices in a cell (written to a file)
%   cmd = task_runner_run('matlab', rc, own, steal)
%       the job script command of a MATLAB job
%   [claimed, live, cache] = task_runner_run('claims', rc, fs, cache)
%       for tasks fs, whether they are claimed, and whether by a process
%       that may still run them: the owner's Slurm job is in the queue, or
%       (no job, or the master's job) its process runs on this node; a
%       claim marked failed is not live. cache (optional) keeps the claims
%       and the owners' jobs read before, for the other batches of a poll
%   task_runner_run('release', rc, fs, b)
%       removes the claims of tasks fs and their failed marks; with b, marks
%       batch b resubmitted (other jobs leave its tasks to its own job)
%   ok = task_runner_run('claim', rc, fs)
%       claims tasks fs for this process (the master runs them itself)
%   rc = task_runner_run('waiting', rc, waiting, t, taskBatchNum)
%       lists for the job script (waiting.txt) the batches whose tasks
%       (waiting: logical per task) have waited in the queue or for their
%       submission since 30 s or more at time t (s); runners start where
%       their jobs can take such tasks
%
% Single jobs (generic_single_job_submit_wrapper): the single jobs this
% MATLAB session submits with the same function, resources and compiled
% master (or MATLAB launch) share a call folder, a pool, in the job log
% folder (from the second of them on; the first runs as before): each job is
% a batch of one task (its own), and tasks are added as their jobs are
% submitted. A job runs its own task (compiled jobs: in a process of its
% own, as before; MATLAB jobs: in a runner), then unstarted tasks of the
% pool's other jobs in the order of their submission (in a task runner)
% while Slurm keeps those jobs in the queue (pending, not just submitted),
% and cancels them; a job whose task another job's runner has ends at once.
% The master takes a task back from the pool before it runs or submits it
% again (as before, without runners).
%   t = task_runner_run('single_task', funcStr)
%       the pool task of funcStr ([] if it has none, or the master took it
%       back): pool key, index i and whether it has a job (submitted)
%   t = task_runner_run('single_add', funcStr, key, parentDir, line, MCCMasterStr, MatlabLaunchStr)
%       adds the task to the pool of key, which its first task makes in
%       parentDir for the compiled master MCCMasterStr, or ('' and the MATLAB
%       launch string) for MATLAB jobs; line is the task line (compiled) or
%       task string (MATLAB); [] if the task can't be added (it then runs
%       as before)
%   cmd = task_runner_run('single_cmd', funcStr)
%       the job script command of the task's job
%   task_runner_run('single_job', funcStr, jobId, arrayId)
%       records the job of the task; jobId [] if it couldn't be submitted
%       (then the master takes the task back)
%   [live, claimed] = task_runner_run('single_live', funcStr)
%       whether a process other than the task's own job has it and may
%       still run it, and whether it is claimed
%   ok = task_runner_run('single_take', funcStr)
%       the master takes the task back (no runner takes it then); false if
%       another process has it
%   status = task_runner_run('single_status', jobId, arrayId, status)
%       the status of a single job whose Slurm job has status (see
%       check_single_job_status): after the job ended, 1 while another
%       process has its task, -1 (after the master took the task back) if
%       none has; while it is pending, -1 if a runner of another job failed
%       the task (the master took it back and cancelled the job), and a
%       runner that ended without finishing it leaves it to the job


persistent S
if isempty(S) && startsWith(action, 'single')
    S = struct('session', '', 'pools', containers.Map('KeyType', 'char', 'ValueType', 'any'), ...
        'tasks', containers.Map('KeyType', 'char', 'ValueType', 'any'), ...
        'jobs', containers.Map('KeyType', 'char', 'ValueType', 'any'));
end

varargout = cell(1, max(1, nargout));
switch action
    case 'waiting'
        varargout{1} = waiting_batches(varargin{:});
    case 'supported'
        varargout{1} = is_supported(varargin{:});
    case 'setup'
        varargout{1} = setup_mcc(varargin{:});
    case 'setup_matlab'
        varargout{1} = setup_matlab(varargin{:});
    case 'keys'
        varargout{1} = batch_keys(varargin{:});
    case 'job'
        varargout{1} = job_cmd(varargin{:});
    case 'matlab'
        varargout{1} = matlab_cmd(varargin{:});
    case 'claims'
        [varargout{1 : max(1, nargout)}] = claim_states(varargin{:});
    case 'release'
        release_claims(varargin{:});
    case 'claim'
        varargout{1} = claim_tasks(varargin{:});
    case 'single_task'
        varargout{1} = single_task(S, varargin{:});
    case 'single_add'
        [varargout{1}, S] = single_add(S, varargin{:});
    case 'single_cmd'
        varargout{1} = single_cmd(S, varargin{:});
    case 'single_job'
        S = single_job(S, varargin{:});
    case 'single_live'
        [varargout{1 : max(1, nargout)}] = single_live(S, varargin{:});
    case 'single_take'
        [varargout{1}, S] = single_take(S, varargin{:});
    case 'single_status'
        [varargout{1}, S] = single_status(S, varargin{:});
    otherwise
        error('task_runner_run:unknownAction', 'Unknown action %s.', action);
end

end


function [tf] = is_supported(MCCMasterStr)
% task runners need bash (Linux) and a compiled master built with them

tf = false;
if ~isunix || ismac
    return;
end
fn = [fileparts(MCCMasterStr), '/task_runner.txt'];
if ~isfile(fn)
    return;
end
v = sscanf(fileread(fn), '%d', 1);
tf = ~isempty(v) && v == task_runner_version();

end


function [rc] = new_rc(D)
% the call's runner settings; D without the trailing slash. The job script
% and its commands take paths without quoting

% a leading ~ is the home folder, as for MATLAB's file functions (the job
% script has D in quotes, where bash doesn't expand it)
if strcmp(D, '~') || startsWith(D, '~/')
    D = [getenv('HOME'), D(2 : end)];
end
D = regexprep(strrep(D, '\', '/'), '/+', '/');
D = regexprep(D, '/$', '');
rc = struct('on', false, 'D', D, 'claimsDir', [D, '/claims'], 'script', [D, '/run_tasks.sh']);
if any(ismember(D, ['''"$`\| ', newline, char(9)]))
    error('task_runner_run:path', 'The call folder %s has characters the job script does not take.', D);
end

end


function [rc] = setup_mcc(D, tlineStrs, MCCMasterStr, taskBatchNum, tmpInfo, excluded, keys)
% compiled mode: task lines, the master's run_mccMaster.sh with a stand-in
% master that writes its arguments, and the job script

if nargin < 7
    keys = {};
end
rc = struct('on', false);
try
    rc = new_rc(D);
    D = rc.D;
    if any(contains(MCCMasterStr, {'''', '|', newline}))
        error('task_runner_run:path', 'The compiled master %s has characters the job script does not take.', MCCMasterStr);
    end
    tmpl = script_template();
    mkdir_new([D, '/claims']);
    mkdir_new([D, '/argdump']);
    copyfile(MCCMasterStr, [D, '/argdump/run_mccMaster.sh']);
    fileattrib([D, '/argdump/run_mccMaster.sh'], '+x');
    write_lines([D, '/argdump/mccMaster'], arg_writer());
    fileattrib([D, '/argdump/mccMaster'], '+x');
    nF = numel(tlineStrs);
    write_task_lines(D, tlineStrs(:));
    write_numbers([D, '/runnable.txt'], find(~cellfun('isempty', tlineStrs(:)) & ~excluded(:)));
    if ~isempty(keys)
        write_lines([D, '/keys.txt'], keys(:));
    end
    localRun = ~isempty(tmpInfo);
    uwt = 0;
    if localRun
        write_lines([D, '/tmps.txt'], cellfun(@(x, y) sprintf('%s\t%s', x, y), tmpInfo.tmpFns(:), tmpInfo.outFns(:), 'unif', 0));
        uwt = tmpInfo.unitWaitTime;
    end
    write_script(rc, tmpl, nF, taskBatchNum, MCCMasterStr, '', localRun, uwt);
    rc.on = true;
catch ME
    warning('task_runner_run:setup', 'The tasks run without task runners: %s', ME.message);
end

end


function [rc] = setup_matlab(D, funcStrs, outputFullpaths, finalOutFullpath, isDone, taskBatchNum, MatlabLaunchStr, runDir, GPUJob, keys)
% MATLAB mode: the task strings (tasks.mat), lines that check the tasks'
% outputs as compiled task lines do, and the job script

if nargin < 10
    keys = {};
end
rc = struct('on', false);
try
    rc = new_rc(D);
    D = rc.D;
    if any(contains([outputFullpaths(:); {finalOutFullpath}], {'"', '$', '`', '\', newline}))
        error('task_runner_run:path', 'Output paths with ", $, `, \ or line breaks are not supported.');
    end
    tmpl = script_template();
    mkdir_new(D);
    mkdir_new([D, '/claims']);
    mkdir_new([D, '/argdump']);
    write_lines([D, '/argdump/run_mccMaster.sh'], {'#!/bin/sh'; 'shift 1'; 'exec "$(dirname "$0")/mccMaster" "$@"'});
    fileattrib([D, '/argdump/run_mccMaster.sh'], '+x');
    write_lines([D, '/argdump/mccMaster'], arg_writer());
    fileattrib([D, '/argdump/mccMaster'], '+x');
    if isempty(runDir)
        runDir = '';
    end
    save([D, '/tasks.mat'], 'funcStrs', 'runDir', 'GPUJob');
    nF = numel(funcStrs);
    lines = repmat({''}, nF, 1);
    for f = reshape(find(~isDone), 1, [])
        if isempty(finalOutFullpath)
            lines{f} = sprintf(['OUT="%s"; if [ ! -f "$OUT" ] && [ ! -d "$OUT" ]; then "%s/argdump/run_mccMaster.sh" x %d ; ', ...
                'else echo output "$OUT" already exists; fi'], outputFullpaths{f}, D, f);
        else
            lines{f} = sprintf(['OUT="%s"; FINALOUT="%s"; if [ ! -f "$OUT" ] && [ ! -d "$OUT" ] && [ ! -f "$FINALOUT" ] ', ...
                '&& [ ! -d "$FINALOUT" ]; then "%s/argdump/run_mccMaster.sh" x %d ; else echo output "$OUT" already exists; fi'], ...
                outputFullpaths{f}, finalOutFullpath, D, f);
        end
    end
    write_task_lines(D, lines);
    write_numbers([D, '/runnable.txt'], find(~isDone(:)));
    if ~isempty(keys)
        write_lines([D, '/keys.txt'], keys(:));
    end
    write_script(rc, tmpl, nF, taskBatchNum, '', MatlabLaunchStr, false, 0);
    rc.on = true;
catch ME
    warning('task_runner_run:setup', 'The tasks run without task runners: %s', ME.message);
end

end


function [keys] = batch_keys(nB, cpus, GPUJob, masterBatch)
% per batch: '<CPUs>' or '<CPUs>gpu' (the jobs' processes), 'master' for
% the master's batch

key = sprintf('%d', cpus);
if GPUJob
    key = [key, 'gpu'];
end
keys = repmat({key}, nB, 1);
if masterBatch >= 1 && masterBatch <= nB
    keys{masterBatch} = 'master';
end

end


function [tmpl] = script_template()
% the job script's text (task_runner.sh, next to task_runner.m; compiled
% masters have it in their archive)

fn = [fileparts(which('task_runner')), '/task_runner.sh'];
if ~isfile(fn) && isdeployed
    dir_info = dir(fullfile(ctfroot, '**', 'task_runner.sh'));
    if ~isempty(dir_info)
        fn = fullfile(dir_info(1).folder, dir_info(1).name);
    end
end
tmpl = fileread(fn);

end


function [lines] = arg_writer()
% the stand-in master: writes the metrics environment of the task line and
% the arguments it gets, NUL separated, to $PK_DUMP

lines = {'#!/bin/bash'
    ['{ printf ''%s\0'' "${PETAKIT5D_TASK_METRICS-}" "${PETAKIT5D_TASK-}" "${PETAKIT5D_TASK_GPU-}"; ', ...
    'printf ''%s\0'' "$@"; } > "$PK_DUMP"']};

end


function [] = write_script(rc, tmpl, nF, taskBatchNum, MCCMasterStr, MatlabLaunchStr, localRun, uwt, waitAll)
% the job script: the call's settings, then the template

if nargin < 9
    waitAll = false;
end
q = @(x) strrep(x, '''', '''\''''');
header = {'#!/bin/bash'
    sprintf('D=''%s''', q(rc.D))
    sprintf('NT=%d', nF)
    sprintf('BS=%d', taskBatchNum)
    sprintf('MASTER=''%s''', q(MCCMasterStr))
    sprintf('MLAUNCH=''%s''', q(MatlabLaunchStr))
    sprintf('LOCALRUN=%d', localRun)
    sprintf('UWTSEC=%d', round(uwt * 60))};
if waitAll
    header{end + 1} = 'WAITALL=1';
end
fid = fopen(rc.script, 'w');
fprintf(fid, '%s\n', header{:});
fwrite(fid, tmpl);
fclose(fid);

end


function [cmd] = job_cmd(rc, own, P, gpu, steal, retries)
% the job script command of a compiled job

cmd = sprintf('bash %s job %s %d %d %d %d', rc.script, own_str(rc, own), P, gpu, steal, retries);

end


function [cmd] = matlab_cmd(rc, own, steal)
% the job script command of a MATLAB job

cmd = sprintf('bash %s matlab %s %d', rc.script, own_str(rc, own), steal);

end


function [s] = own_str(rc, own)
% a batch number, or ({indices}) the file of a list of task indices

if ~iscell(own)
    s = sprintf('%d', own);
    return;
end
[~, name] = fileparts(tempname);
fn = sprintf('%s/own_%s.txt', rc.D, name);
write_numbers(fn, own{1}(:));
s = fn;

end


function [claimed, live, cache] = claim_states(rc, fs, cache)
% the claims of tasks fs and whether their owners may still run them; cache
% has the claims and the owners' jobs read before (one call per poll)

if nargin < 3 || isempty(cache)
    cache = read_claims(rc);
end
fs = fs(:);
claimed = false(size(fs));
live = false(size(fs));
if isempty(fs) || isempty(cache.names)
    return;
end
[isClaim, k] = ismember(arrayfun(@(x) sprintf('t%d', x), fs, 'unif', 0), cache.names);
isFailed = ismember(arrayfun(@(x) sprintf('f%d', x), fs, 'unif', 0), cache.names);
claimed = isClaim;
owners = repmat({''}, numel(fs), 1);
owners(isClaim) = cache.targets(k(isClaim));

% owners: <job>_<array task>:<node>:<pid>:<slot>
tok = regexp(owners, '^(\d+)_(\d+):([^:]*):(\d+):', 'tokens', 'once');
[myJob, myArray] = own_job();
myNode = node_name();
jobKeys = repmat({''}, numel(fs), 1);
for i = reshape(find(isClaim & ~isFailed), 1, [])
    if isempty(tok{i})
        % an owner the job script didn't write: taken as live
        live(i) = true;
        continue;
    end
    job = str2double(tok{i}{1});
    arr = str2double(tok{i}{2});
    if job > 0 && ~(job == myJob && arr == myArray)
        jobKeys{i} = sprintf('%d_%d', job, arr);
    else
        % a process of no job, or of the master's job: one the master
        % started has ended (the master waits for them), others run if
        % their process is alive on this node
        pid = str2double(tok{i}{4});
        live(i) = strcmp(tok{i}{3}, myNode) && pid ~= feature('getpid') && isfolder(['/proc/', tok{i}{4}]);
    end
end
inJob = ~cellfun(@isempty, jobKeys);
if ~any(inJob)
    return;
end
% the owners' jobs not queried yet, with one squeue call
newKeys = unique(jobKeys(inJob));
newKeys = newKeys(~isKey(cache.jobStatus, newKeys));
if ~isempty(newKeys)
    ids = sscanf(strjoin(strrep(newKeys, '_', ' '), ' '), '%d');
    status = check_batch_slurm_jobs_status(ids(1 : 2 : end), ids(2 : 2 : end));
    for i = 1 : numel(newKeys)
        cache.jobStatus(newKeys{i}) = status(i);
    end
end
live(inJob) = cellfun(@(x) cache.jobStatus(x) >= 0, jobKeys(inJob));

end


function [cache] = read_claims(rc)
% the claims folder's entries and their link targets (owners)

cache = struct('names', {{}}, 'targets', {{}}, 'jobStatus', containers.Map('KeyType', 'char', 'ValueType', 'double'));
[status, txt] = system(sprintf('find "%s" -maxdepth 1 -printf ''%%f\\t%%l\\n''', rc.claimsDir));
if status ~= 0
    return;
end
c = regexp(txt, '^(?<name>[tf]\d+)\t(?<target>[^\n]*)$', 'names', 'lineanchors');
if ~isempty(c)
    cache.names = {c.name};
    cache.targets = {c.target};
end

end


function [job, arr] = own_job()
% this process's Slurm job and array task (0 without)

job = str2double(getenv('SLURM_ARRAY_JOB_ID'));
if isnan(job)
    job = str2double(getenv('SLURM_JOB_ID'));
end
arr = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isnan(job)
    job = 0;
end
if isnan(arr)
    arr = 0;
end

end


function [node] = node_name()
% the node name the job script records

node = getenv('SLURMD_NODENAME');
if isempty(node)
    [~, node] = system('uname -n');
    node = strtrim(node);
end

end


function [] = release_claims(rc, fs, b)
% removes the claims and failed marks of tasks fs; marks batch b resubmitted

if ~isempty(fs)
    names = [arrayfun(@(x) sprintf('t%d', x), fs(:)', 'unif', 0), arrayfun(@(x) sprintf('f%d', x), fs(:)', 'unif', 0)];
    for i = 1 : 500 : numel(names)
        system(sprintf('cd "%s" && rm -f %s', rc.claimsDir, strjoin(names(i : min(i + 499, end)), ' ')));
    end
end
if nargin > 2 && ~isempty(b)
    fid = fopen(sprintf('%s/r%d', rc.claimsDir, b), 'w');
    if fid >= 0
        fclose(fid);
    end
end

end


function [rc] = waiting_batches(rc, waiting, t, taskBatchNum)
% writes the batches waiting since 30 s or more (when they change)

if ~isfield(rc, 'waitT0')
    rc.waitT0 = nan(numel(waiting), 1);
    rc.waitBatches = [];
end
waiting = waiting(:);
rc.waitT0(~waiting) = NaN;
rc.waitT0(waiting & isnan(rc.waitT0)) = t;
b = unique(ceil(find(waiting & t - rc.waitT0 >= 30) / taskBatchNum));
if isequal(b(:), rc.waitBatches(:))
    return;
end
rc.waitBatches = b(:);
try
    fn = [rc.D, '/waiting.txt'];
    write_numbers([fn, '.tmp'], b(:));
    movefile([fn, '.tmp'], fn);
catch ME
    warning('task_runner_run:waiting', 'The list of waiting batches could not be written: %s', ME.message);
end

end


function [ok] = claim_tasks(rc, fs)
% claims tasks fs for this process; ok(i) is false where another process
% has task fs(i)

fs = fs(:);
ok = false(size(fs));
[job, arr] = own_job();
owner = sprintf('%d_%d:%s:%d:0', job, arr, node_name(), feature('getpid'));
[status, txt] = system(sprintf('cd "%s" && for f in %s; do ln -s "%s" "t$f" 2> /dev/null && echo "$f"; done', ...
    rc.claimsDir, sprintf('%d ', fs), owner));
if status ~= 0
    return;
end
ok = ismember(fs, sscanf(txt, '%d'));

end


function [t] = single_task(S, funcStr)
% the pool task of funcStr, [] if it has none or the master took it back

t = [];
if isKey(S.tasks, funcStr)
    t = S.tasks(funcStr);
    if strcmp(t.state, 'taken')
        t = [];
    end
end

end


function [t, S] = single_add(S, funcStr, key, parentDir, line, MCCMasterStr, MatlabLaunchStr)
% adds the task (its task line or task string) to the pool of key, made
% with its first task

t = [];
if isKey(S.tasks, funcStr)
    t = single_task(S, funcStr);
    return;
end
% (one line per task in tasks.txt; the job script reads them by offset)
if any(line == newline | line == char(13))
    return;
end
try
    if isKey(S.pools, key)
        p = S.pools(key);
        if ~p.on
            return;
        end
    else
        if isempty(S.session)
            % the time (also its microseconds) and this MATLAB's process id
            % (a first datetime, uuid or RandStream call took 70-250 ms,
            % before the first job's submission)
            S.session = sprintf('%04d%02d%02d_%02d%02d%02d_%d_%06d', fix(clock), feature('getpid'), mod(tic, 1e6));
        end
        p = single_pool(sprintf('%s/task_runner_single_%s_%d', parentDir, S.session, S.pools.Count + 1), ...
            MCCMasterStr, MatlabLaunchStr);
    end
    % the line, its byte offset, then the task's index: the job script takes
    % only indices in runnable.txt, so it never reads a line being written
    i = p.n + 1;
    bytes = unicode2native([line, newline], 'UTF-8');
    append_text([p.rc.D, '/tasks.txt'], bytes);
    append_text([p.rc.D, '/tasks.idx'], sprintf('%012d\n', p.bytes));
    append_text([p.rc.D, '/runnable.txt'], sprintf('%d\n', i));
    p.n = i;
    p.bytes = p.bytes + numel(bytes);
    S.pools(key) = p;
    t = struct('key', key, 'i', i, 'state', 'new', 'jobKey', '');
    S.tasks(funcStr) = t;
catch ME
    warning('task_runner_run:single', 'The single job runs without task runners: %s', ME.message);
    % (the pool's files may be incomplete: its later tasks run as before)
    if exist('p', 'var') && isstruct(p)
        p.on = false;
    else
        p = struct('on', false);
    end
    S.pools(key) = p;
    t = [];
end

end


function [p] = single_pool(D, MCCMasterStr, MatlabLaunchStr)
% the call folder of a pool: the job script for batches of one task, whose
% count grows (NT), with all of a pool's tasks counted as waiting (their
% jobs are submitted one by one, so tasks a job can take have jobs that
% have not started), the stand-in master and the (empty) task files

rc = new_rc(D);
D = rc.D;
if any(contains(MCCMasterStr, {'''', '|', newline}))
    error('task_runner_run:path', 'The compiled master %s has characters the job script does not take.', MCCMasterStr);
end
tmpl = script_template();
% the folders, a probe that the claims folder takes symbolic links (claims
% are symbolic links), the empty task files and the compiled master's
% run_mccMaster.sh, in one shell (MATLAB's mkdir and copyfile took 70 ms
% more at their first use, before the pool's first job was submitted)
probe = sprintf('%s/claims/.probe_%d', D, feature('getpid'));
cmd = sprintf(['mkdir -p "%s/claims" "%s/argdump" && ln -s probe "%s" && rm -f "%s" ', ...
    '&& : > "%s/tasks.txt" && : > "%s/tasks.idx" && : > "%s/runnable.txt"'], D, D, probe, probe, D, D, D);
if ~isempty(MCCMasterStr)
    cmd = sprintf('%s && cp ''%s'' "%s/argdump/run_mccMaster.sh" && chmod +x "%s/argdump/run_mccMaster.sh"', ...
        cmd, MCCMasterStr, D, D);
end
[status, out] = system(cmd);
if status ~= 0
    error('task_runner_run:claims', 'The pool folder %s can''t be made, or its claims folder takes no symbolic links: %s', ...
        D, strtrim(out));
end
if isempty(MCCMasterStr)
    write_lines([D, '/argdump/run_mccMaster.sh'], {'#!/bin/sh'; 'shift 1'; 'exec "$(dirname "$0")/mccMaster" "$@"'});
    fileattrib([D, '/argdump/run_mccMaster.sh'], '+x');
    % the MATLAB runner takes each task string from the arguments its line
    % gives the stand-in master
    funcStrs = {};
    runDir = '';
    GPUJob = false;
    fromDump = true;
    save([D, '/tasks.mat'], 'funcStrs', 'runDir', 'GPUJob', 'fromDump');
end
write_lines([D, '/argdump/mccMaster'], arg_writer());
fileattrib([D, '/argdump/mccMaster'], '+x');
write_script(rc, tmpl, 999999999, 1, MCCMasterStr, MatlabLaunchStr, false, 0, true);
rc.on = true;
p = struct('rc', rc, 'mcc', ~isempty(MCCMasterStr), 'n', 0, 'bytes', 0, 'on', true);

end


function [cmd] = single_cmd(S, funcStr)
% the job script command of a single job: its own task first (the job ends
% if another process has it), then other tasks of the pool whose jobs wait
% in the queue, in the order of their submission (as the queue would start
% their jobs, so that results come in about the order of before); one try
% per task (the caller submits a failed task again)

t = S.tasks(funcStr);
p = S.pools(t.key);
if p.mcc
    cmd = job_cmd(p.rc, t.i, 1, 0, 3, 1);
else
    cmd = matlab_cmd(p.rc, t.i, 3);
end

end


function [S] = single_job(S, funcStr, jobId, arrayId)
% records the job of the task, also for the runner that takes the task
% (jobs.txt: it cancels the job if it is pending); without one, the master
% takes the task back

if isempty(jobId)
    [~, S] = single_take(S, funcStr);
    return;
end
t = S.tasks(funcStr);
t.state = 'submitted';
t.jobKey = sprintf('%d_%d', jobId, arrayId);
S.tasks(funcStr) = t;
S.jobs(t.jobKey) = funcStr;
p = S.pools(t.key);
try
    append_text([p.rc.D, '/jobs.txt'], sprintf('%d %s\n', t.i, t.jobKey));
catch ME
    % (the job then starts and ends at once if a runner took its task)
    warning('task_runner_run:single', 'The job of a single job''s task is not recorded: %s', ME.message);
end

end


function [live, claimed] = single_live(S, funcStr)
% whether a process other than the task's own job has its claim and may
% still run it, and whether the task is claimed

t = S.tasks(funcStr);
p = S.pools(t.key);
[owner, failed] = read_claim(p.rc, t.i);
live = owner_live(owner, failed, t.jobKey);
claimed = ~isempty(owner);

end


function [ok, S] = single_take(S, funcStr)
% the master takes the task back: it removes a claim whose owner ended (and
% the failed mark), claims the task and marks it resubmitted, so that no
% runner takes it; false if another process has it

t = S.tasks(funcStr);
ok = strcmp(t.state, 'taken');
if ok
    return;
end
p = S.pools(t.key);
[owner, failed] = read_claim(p.rc, t.i);
if owner_live(owner, failed, t.jobKey)
    return;
end
% a claim whose owner ended is replaced at once (rename), so that no runner
% claims the task in between (it stays as it is until then: only its owner,
% which ended, and this master remove claims); without a claim, the master's
% fails if a runner just claimed the task
[job, arr] = own_job();
me = sprintf('%d_%d:%s:%d:0', job, arr, node_name(), feature('getpid'));
if isempty(owner)
    claimCmd = sprintf('ln -s "%s" t%d 2> /dev/null', me, t.i);
else
    tmpName = sprintf('.take_t%d_%d', t.i, feature('getpid'));
    claimCmd = sprintf('rm -f %s && ln -s "%s" %s && mv -fT %s t%d', tmpName, me, tmpName, tmpName, t.i);
end
[status, txt] = system(sprintf('cd "%s" && rm -f f%d && %s && : > r%d && echo taken', p.rc.claimsDir, t.i, claimCmd, t.i));
ok = status == 0 && contains(txt, 'taken');
if ~ok
    % another process just claimed it, or (if no live process has it now)
    % the claims folder failed us, and then no runner can claim it either
    [owner, failed] = read_claim(p.rc, t.i);
    ok = ~owner_live(owner, failed, t.jobKey);
end
if ok
    t.state = 'taken';
    S.tasks(funcStr) = t;
end

end


function [status, S] = single_status(S, jobId, arrayId, status)
% the status of a pool task's own job (status from Slurm): after the job
% ended, 1 while another process has the task, else the master takes it
% back (-1); while the job is pending, see single_pending (-1 after a
% runner failed the task: the master took it back and cancelled the job)

if status == 1
    return;
end
jobKey = sprintf('%d_%d', jobId, arrayId);
if ~isKey(S.jobs, jobKey)
    return;
end
funcStr = S.jobs(jobKey);
t = S.tasks(funcStr);
if strcmp(t.state, 'taken') || ~strcmp(t.jobKey, jobKey)
    return;
end
if status == 0
    [state, S] = single_pending(S, funcStr);
    if strcmp(state, 'failed')
        [ok, S] = single_take(S, funcStr);
        if ok
            system(sprintf('scancel %s', jobKey));
            status = -1;
        end
    end
    return;
end
[ok, S] = single_take(S, funcStr);
if ~ok
    status = 1;
end

end


function [state, S] = single_pending(S, funcStr)
% while the task's own job is pending: 'live' if a runner of another job
% has the task, 'failed' if one failed it, 'released' if one had it and
% ended without finishing it (killed): its claim is removed, so that the own
% job runs the task when it starts (or another runner takes it), as without
% runners; else 'none'. The pool's claims and their owners' jobs are read
% at most every 10 s (only this master removes claims of processes that
% ended, so such a claim is still there)

state = 'none';
t = S.tasks(funcStr);
p = S.pools(t.key);
if ~isfield(p, 'claimCache') || isempty(p.claimCache) || toc(p.claimTime) > 10
    p.claimCache = [];
    p.claimTime = tic;
end
[claimed, live, p.claimCache] = claim_states(p.rc, t.i, p.claimCache);
if ~claimed
    S.pools(t.key) = p;
    return;
end
if live
    state = 'live';
elseif ismember(sprintf('f%d', t.i), p.claimCache.names)
    state = 'failed';
else
    state = 'released';
    system(sprintf('rm -f "%s/t%d"', p.rc.claimsDir, t.i));
    p.claimCache.names(strcmp(p.claimCache.names, sprintf('t%d', t.i))) = {''};
end
S.pools(t.key) = p;

end


function [owner, failed] = read_claim(rc, i)
% the owner of task i's claim ('' if none) and whether it is marked failed;
% opening the claims folder first revalidates it on NFS, so that claims
% other nodes just made are seen

[~, txt] = system(sprintf('C="%s"; : < "$C"; printf ''t %%s\\n'' "$(readlink "$C/t%d" 2> /dev/null)"; [ -L "$C/f%d" ] && echo f', ...
    rc.claimsDir, i, i));
lines = splitlines(txt);
k = find(startsWith(lines, 't '), 1);
owner = '';
if ~isempty(k)
    owner = strtrim(lines{k}(3 : end));
end
failed = any(strcmp(strtrim(lines), 'f'));

end


function [live] = owner_live(owner, failed, ownJobKey)
% whether the owner of a claim may still run its task: its Slurm job is in
% the queue or (no job, or this process's job) its process runs on this
% node; not a claim marked failed, of the task's own job (the caller has
% its status) or of this process

live = false;
if isempty(owner) || failed
    return;
end
tok = regexp(owner, '^(\d+)_(\d+):([^:]*):(\d+):', 'tokens', 'once');
if isempty(tok)
    % an owner the job script didn't write: taken as live
    live = true;
    return;
end
if strcmp(sprintf('%s_%s', tok{1}, tok{2}), ownJobKey)
    return;
end
job = str2double(tok{1});
arr = str2double(tok{2});
[myJob, myArray] = own_job();
if job > 0 && ~(job == myJob && arr == myArray)
    live = check_batch_slurm_jobs_status(job, arr) >= 0;
else
    live = strcmp(tok{3}, node_name()) && str2double(tok{4}) ~= feature('getpid') && isfolder(['/proc/', tok{4}]);
end

end


function [] = append_text(fn, x)
% appends text or bytes to a file

fid = fopen(fn, 'a');
if fid < 0
    error('task_runner_run:write', 'Cannot write %s.', fn);
end
fwrite(fid, x);
if fclose(fid) ~= 0
    error('task_runner_run:write', 'Cannot write %s.', fn);
end

end


function [] = mkdir_new(d)

if ~isfolder(d)
    mkdir(d);
end

end


function [] = write_task_lines(D, lines)
% the task lines (tasks.txt) and the byte offset of each (tasks.idx: 12
% digits and a line break per task), so that the job script reads a task's
% line without reading the lines before it

bytes = unicode2native([strjoin(lines(:)', newline), newline], 'UTF-8');
offsets = [0, find(bytes(1 : end - 1) == 10)];
fid = fopen([D, '/tasks.txt'], 'w');
if fid < 0
    error('task_runner_run:write', 'Cannot write %s/tasks.txt.', D);
end
fwrite(fid, bytes, 'uint8');
fclose(fid);
fid = fopen([D, '/tasks.idx'], 'w');
if fid < 0
    error('task_runner_run:write', 'Cannot write %s/tasks.idx.', D);
end
fprintf(fid, '%012d\n', offsets);
fclose(fid);

end


function [] = write_numbers(fn, x)
% integers, one per line

fid = fopen(fn, 'w');
if fid < 0
    error('task_runner_run:write', 'Cannot write %s.', fn);
end
fprintf(fid, '%d\n', x);
fclose(fid);

end


function [] = write_lines(fn, lines)
% lines of text, each ending with a line break

fid = fopen(fn, 'w');
if fid < 0
    error('task_runner_run:write', 'Cannot write %s.', fn);
end
if ~isempty(lines)
    % (fprintf would skip empty lines)
    fwrite(fid, [strjoin(lines, newline), newline]);
end
fclose(fid);

end
