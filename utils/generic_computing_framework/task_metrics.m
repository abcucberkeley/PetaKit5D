function [out] = task_metrics(action, varargin)
% records the resources of a job framework task when the framework's
% taskMetrics option is on (task_metrics_run sets up the metrics of a
% framework call, read_task_metrics reads them). Each process that runs
% tasks appends one JSON line per task to its own file in the metrics folder
% of the call. Errors are ignored (one warning per process), so the metrics
% never stop or fail a task.
%
% A record has the function, mode (local, master, pool, slurm_matlab,
% slurm_mcc, mcc), node, process id, Slurm job and array task, GNU parallel
% slot, start and end (seconds since 1970, UTC), wall and CPU time (all
% threads), exit status and error; on Linux also user and system CPU, the
% CPU of waited-for child processes, resident memory at start and end and its
% peak during the task (VmHWM, reset at the start through
% /proc/self/clear_refs), bytes read and written, page faults, the CPUs the
% process may use and the process start (for start-up costs); the GPU for
% GPU tasks.
%
% Compiled tasks run this code in a new MATLAB Runtime each, where every
% first call costs: the code avoids constructs that are slow the first time
% (^ loads a math library, ~22 ms; setenv, cell arrays, ~3-5 ms).
%
% Usage:
%   tm = task_metrics('start', funcName, runDir, taskInd, mode, GPUJob)
%       the state before a task ([] if runDir is empty). Compiled tasks give
%       only the function name: their settings are in the environment their
%       task line sets: PETAKIT5D_TASK_METRICS (the metrics folder),
%       PETAKIT5D_TASK (task index) and PETAKIT5D_TASK_GPU (1 for GPU tasks)
%   task_metrics('stop', tm, err)
%       appends the record of the task; err is its MException or exit status
%       (empty or 0 if it finished)
%   p = task_metrics('process')
%       this process's node and id, and the task running in it (current)


out = [];
try
    switch action
        case 'start'
            out = start_task(varargin{:});
        case 'stop'
            stop_task(varargin{:});
        case 'process'
            out = process_info();
            out.current = current_task();
    end
catch ME
    warn_once(ME);
end

end


function [tm] = start_task(funcName, runDir, taskInd, mode, GPUJob)
% the state of this process before a task

t0 = tic;
tm = [];
if nargin < 2
    % a compiled task: the settings are in its environment
    runDir = getenv('PETAKIT5D_TASK_METRICS');
    if isempty(runDir)
        return;
    end
    taskInd = sscanf(getenv('PETAKIT5D_TASK'), '%f');
    GPUJob = strcmp(getenv('PETAKIT5D_TASK_GPU'), '1');
    mode = 'mcc';
    if ~isempty(getenv('SLURM_JOB_ID'))
        mode = 'slurm_mcc';
    end
else
    if isempty(runDir)
        return;
    end
    if nargin < 3
        taskInd = [];
    end
    if nargin < 4
        mode = '';
    end
    if nargin < 5
        GPUJob = false;
    end
end

linux = isunix && ~ismac;
tm.proc0 = proc_snapshot(linux);
p = process_info(tm.proc0(11));
k = find(funcName == '(', 1);
if ~isempty(k)
    funcName = strtrim(funcName(1 : k - 1));
end
if isempty(taskInd)
    taskInd = NaN;
end
% the run is the name of the metrics folder
d = runDir(1 : find(runDir ~= '/' & runDir ~= '\', 1, 'last'));
tm.runId = d(find(d == '/' | d == '\', 1, 'last') + 1 : end);
tm.runDir = runDir;
tm.task = taskInd;
tm.func = funcName;
tm.mode = mode;
tm.gpu = GPUJob;
tm.threads = maxNumCompThreads;
tm.parent = current_task();
if p.linux
    % the peak resident memory from here on is the task's
    fid = fopen('/proc/self/clear_refs', 'w');
    if fid >= 0
        fwrite(fid, '5');
        fclose(fid);
    end
end
current_task(sprintf('%s/%d', tm.runId, taskInd(1)));
tm.metricsSec = toc(t0);
tm.cpu0 = cputime;
tm.wall0 = tic;

end


function [] = stop_task(tm, err)
% appends the record of the task to this process's file in the metrics folder

if isempty(tm)
    return;
end
wallSec = toc(tm.wall0);
cpuSec = cputime - tm.cpu0;
t0 = tic;
if nargin < 2
    err = [];
end
p = process_info();
current_task(tm.parent);
v0 = tm.proc0;
v1 = proc_snapshot(p.linux);

status = 'ok';
exitCode = 0;
errStr = '';
if isnumeric(err) && ~isempty(err) && err(1) ~= 0
    status = 'failed';
    exitCode = double(err(1));
elseif ~isnumeric(err)
    status = 'error';
    exitCode = 1;
    errStr = err.message;
    if ~isempty(err.identifier)
        errStr = [err.identifier, ': ', errStr];
    end
    errStr = errStr(1 : min(end, 1000));
end
gpuName = '';
gpuMem = [NaN, NaN];
cudaDevices = getenv('CUDA_VISIBLE_DEVICES');
if tm.gpu || ~isempty(cudaDevices)
    [gpuName, gpuMem] = gpu_info();
end

% memory in GB from kB, CPU in s from clock ticks of 10 ms (see proc_snapshot);
% the numbers in one sprintf, with null for NaN (no field name has NaN or Inf)
startTime = epoch_seconds(tm.wall0, p);
numFmt = ['"version":1,"task":%.16g,"batchSize":%.16g,"exitCode":%.16g,"pid":%.16g,"slot":%.16g,', ...
    '"seq":%.16g,"start":%.16g,"end":%.16g,"wallSec":%.16g,"cpuSec":%.16g,"userSec":%.16g,', ...
    '"sysSec":%.16g,"childCpuSec":%.16g,"rssStartGB":%.16g,"rssEndGB":%.16g,"rssPeakGB":%.16g,', ...
    '"rssPeakBeforeGB":%.16g,"readBytes":%.16g,"writeBytes":%.16g,"rchar":%.16g,"wchar":%.16g,', ...
    '"minorFaults":%.16g,"majorFaults":%.16g,"threads":%.16g,"cpus":%.16g,"slurmCpus":%.16g,', ...
    '"osThreads":%.16g,"gpuMemGB":%.16g,"gpuFreeGB":%.16g,"processStart":%.16g,"memTotalGB":%.16g,', ...
    '"metricsSec":%.16g}'];
nums = [tm.task(1), numel(tm.task), exitCode, p.pid, env_num('PARALLEL_JOBSLOT'), env_num('PARALLEL_SEQ'), ...
    startTime, startTime + wallSec, wallSec, cpuSec, (v1(7) - v0(7)) / 100, (v1(8) - v0(8)) / 100, ...
    (v1(9) + v1(10) - v0(9) - v0(10)) / 100, v0(1) / 1048576, v1(1) / 1048576, v1(2) / 1048576, ...
    v0(2) / 1048576, v1(14) - v0(14), v1(15) - v0(15), v1(12) - v0(12), v1(13) - v0(13), v1(5) - v0(5), ...
    v1(6) - v0(6), tm.threads, v0(4), env_num('SLURM_CPUS_PER_TASK'), v1(3), gpuMem(1), gpuMem(2), ...
    p.processStart, p.memTotalGB, tm.metricsSec + toc(t0)];
js = ['{"kind":"task","runId":"', json_str(tm.runId), '","func":"', json_str(tm.func), ...
    '","mode":"', tm.mode, '","status":"', status, '","error":"', json_str(errStr), ...
    '","node":"', json_str(p.host), '","jobId":"', getenv('SLURM_JOB_ID'), ...
    '","arrayJobId":"', getenv('SLURM_ARRAY_JOB_ID'), '","arrayTaskId":"', getenv('SLURM_ARRAY_TASK_ID'), ...
    '","gpu":"', json_str(gpuName), '","cudaDevices":"', json_str(cudaDevices), ...
    '","parent":"', json_str(tm.parent), '",', strrep(strrep(sprintf(numFmt, nums), 'Inf', 'null'), 'NaN', 'null')];
append_line([tm.runDir, '/', p.file], js);

end


function [v] = proc_snapshot(linux)
% resident memory and its peak (kB), threads, allowed CPUs, minor and major
% page faults, user, system and waited-for children's CPU (clock ticks of
% 10 ms), the process start (ticks after boot) and bytes read and written
% (rchar, wchar, read_bytes, write_bytes), from /proc/self; NaN elsewhere

v = NaN(1, 15);
if ~linux
    return;
end
s = read_text('/proc/self/status');
v(1) = status_field(s, 'VmRSS:');
v(2) = status_field(s, 'VmHWM:');
v(3) = status_field(s, 'Threads:');
k = strfind(s, 'Cpus_allowed:');
if ~isempty(k)
    % the CPUs set in the hexadecimal mask
    mask = lower(sscanf(s(k(1) + 13 : end), '%s', 1));
    d = double(mask(mask ~= ',')) - 48;
    d(d > 9) = d(d > 9) - 39;
    bits = [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4];
    v(4) = sum(bits(d(d >= 0 & d <= 15) + 1));
end
s = read_text('/proc/self/stat');
k = find(s == ')', 1, 'last');
if ~isempty(k) && numel(s) > k + 4
    % fields from the 4th on, after the name and the state
    x = sscanf(s(k + 4 : end), '%f');
    if numel(x) >= 19
        v(5 : 11) = x([7, 9, 11, 12, 13, 14, 19]);
    end
end
x = sscanf(read_text('/proc/self/io'), 'rchar: %f wchar: %f syscr: %f syscw: %f read_bytes: %f write_bytes: %f');
if numel(x) == 6
    v(12 : 15) = x([1, 2, 5, 6]);
end

end


function [x] = status_field(s, name)
% a number after name in /proc/self/status, NaN if missing

x = NaN;
k = strfind(s, name);
if ~isempty(k)
    y = sscanf(s(k(1) + numel(name) : end), '%f', 1);
    if ~isempty(y)
        x = y;
    end
end

end


function [p] = process_info(startTicks)
% this process's node, id and record file, found once (startTicks: the
% process start in clock ticks after boot, if known)

persistent info
if ~isempty(info)
    p = info;
    return;
end
p.linux = isunix && ~ismac;
if nargin < 1
    v = proc_snapshot(p.linux);
    startTicks = v(11);
end
p.pid = feature('getpid');
% MATLAB's tic counts microseconds since 1970 on Linux
t = double(tic) * 1e-6;
p.ticEpoch = t > 1e9 && t < 5e9;
p.host = getenv('SLURMD_NODENAME');
if isempty(p.host) && p.linux
    p.host = strtrim(read_text('/proc/sys/kernel/hostname'));
end
if isempty(p.host)
    p.host = getenv('COMPUTERNAME');
end
if isempty(p.host)
    p.host = getenv('HOSTNAME');
end
if isempty(p.host)
    [~, host] = system('hostname');
    p.host = strtrim(host);
end
p.processStart = NaN;
p.memTotalGB = NaN;
nowSec = epoch_seconds(tic, p);
token = sprintf('%d', round(nowSec * 1000));
if p.linux
    % the process start from its start in clock ticks after boot
    uptime = sscanf(read_text('/proc/uptime'), '%f', 1);
    if ~isnan(startTicks) && ~isempty(uptime)
        p.processStart = nowSec - (uptime - startTicks / 100);
        token = sprintf('%d', startTicks);
    end
    p.memTotalGB = status_field(read_text('/proc/meminfo', 256), 'MemTotal:') / 1048576;
end
% file names: letters, digits, '_', '-' and '.'
host = p.host;
host(~((host >= 'a' & host <= 'z') | (host >= 'A' & host <= 'Z') | (host >= '0' & host <= '9') ...
    | host == '_' | host == '-' | host == '.')) = '_';
p.file = sprintf('task_%s_%d_%s.ndjson', host, p.pid, token);
info = p;

end


function [t] = epoch_seconds(ticValue, p)
% seconds since 1970 (UTC) at the time of ticValue

if p.ticEpoch
    t = double(ticValue) * 1e-6;
else
    t = posixtime(datetime('now', 'TimeZone', 'local')) - toc(ticValue);
end

end


function [id] = current_task(newId)
% the task running in this process, the parent of framework calls it makes

persistent cur
if isempty(cur)
    cur = '';
end
id = cur;
if nargin > 0
    cur = newId;
end

end


function [name, mem] = gpu_info()
% the name, memory and free memory (GB) of the GPU the task used, without
% selecting one

name = '';
mem = [NaN, NaN];
try
    if parallel.internal.gpu.isAnyDeviceSelected
        d = gpuDevice;
        name = d.Name;
        mem = [d.TotalMemory, d.AvailableMemory] / 1073741824;
    end
catch
end

end


function [x] = env_num(name)

x = sscanf(getenv(name), '%f', 1);
if isempty(x)
    x = NaN;
end

end


function [s] = json_str(s)
% a string for JSON, escaped

if any(s == '\' | s == '"' | s < 32)
    s = strrep(strrep(s, '\', '\\'), '"', '\"');
    s = strrep(strrep(strrep(s, char(10), '\n'), char(13), '\r'), char(9), '\t');
    s(s < 32) = ' ';
end

end


function [s] = read_text(fn, n)
% the contents (at most n characters) of a text file, '' if it can't be read

if nargin < 2
    n = inf;
end
s = '';
fid = fopen(fn, 'r');
if fid >= 0
    s = fread(fid, [1, n], '*char');
    fclose(fid);
end

end


function [] = append_line(fn, s)
% appends a line to a file only this process writes

fid = fopen(fn, 'a');
if fid < 0
    % the folder is gone (or not yet seen on this node)
    [~, ~] = mkdir(fileparts(fn));
    fid = fopen(fn, 'a');
    if fid < 0
        error('task_metrics:write', 'Cannot write %s.', fn);
    end
end
fprintf(fid, '%s\n', s);
fclose(fid);

end


function [] = warn_once(ME)
% the first metrics error of this process, as a warning

persistent warned
if isempty(warned)
    warned = true;
    warning('task_metrics:failed', 'Task metrics are incomplete: %s', ME.message);
end

end
