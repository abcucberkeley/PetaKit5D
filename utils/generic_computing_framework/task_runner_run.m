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
tf = ~isempty(v) && v == task_runner('version');

end


function [rc] = new_rc(D)
% the call's runner settings; D without the trailing slash. The job script
% and its commands take paths without quoting

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


function [] = write_script(rc, tmpl, nF, taskBatchNum, MCCMasterStr, MatlabLaunchStr, localRun, uwt)
% the job script: the call's settings, then the template

q = @(x) strrep(x, '''', '''\''''');
header = {'#!/bin/bash'
    sprintf('D=''%s''', q(rc.D))
    sprintf('NT=%d', nF)
    sprintf('BS=%d', taskBatchNum)
    sprintf('MASTER=''%s''', q(MCCMasterStr))
    sprintf('MLAUNCH=''%s''', q(MatlabLaunchStr))
    sprintf('LOCALRUN=%d', localRun)
    sprintf('UWTSEC=%d', round(uwt * 60))};
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
