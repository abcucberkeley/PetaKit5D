function [taskTable, runTable] = read_task_metrics(metricsPath, varargin)
% collects the task metrics the job framework records with its taskMetrics
% option (see task_metrics) into tables: taskTable has a row per task run (a
% task that ran twice has two), runTable a row per framework call (or per
% MATLAB session for single jobs). metricsPath is the metrics folder of a
% run, a task_metrics folder or a job log folder (all runs in it).
%
% Besides the recorded fields (times in seconds since 1970, memory in GB;
% compiled tasks the master ran itself have mode 'master', or 'local'
% without a cluster), task rows have
%         startupSec : process start to the start of its first task, for
%                      tasks in their own process (MATLAB or Runtime start;
%                      setup in MATLAB jobs); NaN otherwise
%   lineSec, exitSec : for compiled tasks GNU parallel ran, the run time of
%                      the task line and the time from the end of the task
%                      to the end of its process (Runtime exit), from the
%                      joblog. Task lines without a record get a row with
%                      status 'skipped' (exit status 0, e.g., the output
%                      existed), 'failed' or 'killed' (signal)
%     input, output, cpusPerTask, memAllocate : from the framework call
%   jobSubmit, jobStart, jobEnd, jobState, queueSec, jobHeadSec, jobTailSec :
%                      the task's Slurm job (sacct), with slurmAccounting:
%                      its times, its state, its wait in the queue, its start
%                      to its first task and its last task to its end (MATLAB
%                      or Runtime exit; jobs the framework cancels once all
%                      outputs exist end early). sacct has whole seconds, so
%                      these are good to about a second (and can be negative)
%
% Parameters (as 'specifier'-value pairs):
%   slurmAccounting : true|false (default: false). Add the Slurm job of each task from sacct.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('metricsPath', @ischar);
ip.addParameter('slurmAccounting', false, @islogical);
ip.parse(metricsPath, varargin{:});

pr = ip.Results;
slurmAccounting = pr.slurmAccounting;

runDirs = find_run_dirs(metricsPath);
tasks = cell(0, 1);
runs = cell(0, 1);
lines = cell(0, 1);
for i = 1 : numel(runDirs)
    % the folder name (function names in it can have dots)
    runId = runDirs{i}(1 : end - 1);
    runId = runId(find(runId == '/', 1, 'last') + 1 : end);
    runRecs = read_records([runDirs{i}, 'run.ndjson']);
    taskRecs = {};
    dir_info = dir([runDirs{i}, 'task_*.ndjson']);
    for j = 1 : numel(dir_info)
        taskRecs = [taskRecs; read_records([runDirs{i}, dir_info(j).name])]; %#ok<AGROW>
    end
    dir_info = dir([runDirs{i}, 'parallel_*.joblog']);
    for j = 1 : numel(dir_info)
        lines = [lines; read_joblog([runDirs{i}, dir_info(j).name], runId)]; %#ok<AGROW>
    end
    tasks = [tasks; taskRecs(:)]; %#ok<AGROW>
    runs{end + 1, 1} = struct('runId', runId, 'folder', runDirs{i}, 'records', {runRecs}, 'nRecords', numel(taskRecs)); %#ok<AGROW>
end

taskTable = task_table(tasks);
taskTable = add_startup(taskTable);
taskTable = add_joblog(taskTable, lines);
[taskTable, runTable] = add_run_settings(taskTable, runs);
if slurmAccounting
    taskTable = add_slurm_accounting(taskTable);
end

end


function [runDirs] = find_run_dirs(metricsPath)
% the run folders: metricsPath itself, or the runs in it or in its task_metrics folder

p = [strrep(metricsPath, '\', '/'), '/'];
p = strrep(p, '//', '/');
if ~isempty(dir([p, 'run.ndjson'])) || ~isempty(dir([p, 'task_*.ndjson']))
    runDirs = {p};
    return;
end
if exist([p, 'task_metrics'], 'dir')
    p = [p, 'task_metrics/'];
end
dir_info = dir(p);
dir_info = dir_info([dir_info.isdir] & ~startsWith({dir_info.name}, '.'));
runDirs = cellfun(@(x) [p, x, '/'], {dir_info.name}, 'unif', 0);

end


function [recs] = read_records(fn)
% the JSON records of a file, one per line; lines that don't decode (e.g., cut
% off by a killed process) are left out

recs = cell(0, 1);
if ~exist(fn, 'file')
    return;
end
strs = splitlines(fileread(fn));
for i = 1 : numel(strs)
    if isempty(strtrim(strs{i}))
        continue;
    end
    try
        recs{end + 1, 1} = jsondecode(strs{i}); %#ok<AGROW>
    catch
    end
end

end


function [rows] = read_joblog(fn, runId)
% the lines of a GNU parallel joblog with the task index from their metrics
% environment and the Slurm job from the file name (parallel_<job>_<array>)

rows = cell(0, 1);
[~, name] = fileparts(fn);
jobId = regexp(name, '^parallel_(\d+)_', 'tokens', 'once');
if isempty(jobId)
    jobId = {''};
end
jobId = string(jobId);
strs = splitlines(fileread(fn));
for i = 1 : numel(strs)
    c = strsplit(strs{i}, char(9));
    if numel(c) < 9 || strcmp(c{1}, 'Seq')
        continue;
    end
    task = regexp(c{9}, 'PETAKIT5D_TASK=(\d+)', 'tokens', 'once');
    if isempty(task)
        continue;
    end
    rows{end + 1, 1} = struct('runId', runId, 'task', str2double(task{1}), 'jobId', jobId{1}, ...
        'lineStart', str2double(c{3}), 'lineSec', str2double(c{4}), 'exitCode', str2double(c{7}), ...
        'signal', str2double(c{8})); %#ok<AGROW>
end

end


function [t] = task_table(recs)
% the task records as a table, missing fields as NaN or ""

strNames = {'runId', 'func', 'mode', 'status', 'error', 'node', 'jobId', 'arrayJobId', 'arrayTaskId', ...
    'gpu', 'cudaDevices', 'parent'};
numNames = {'task', 'batchSize', 'exitCode', 'pid', 'slot', 'seq', 'start', 'end', 'wallSec', 'cpuSec', ...
    'userSec', 'sysSec', 'childCpuSec', 'rssStartGB', 'rssEndGB', 'rssPeakGB', 'rssPeakBeforeGB', 'readBytes', ...
    'writeBytes', 'rchar', 'wchar', 'minorFaults', 'majorFaults', 'threads', 'cpus', 'slurmCpus', 'osThreads', ...
    'gpuMemGB', 'gpuFreeGB', 'processStart', 'memTotalGB', 'metricsSec'};
order = {'runId', 'task', 'func', 'mode', 'status', 'exitCode', 'error', 'node', 'pid', 'jobId', 'arrayJobId', ...
    'arrayTaskId', 'slot', 'seq', 'start', 'end', 'wallSec', 'cpuSec', 'userSec', 'sysSec', 'childCpuSec', ...
    'rssStartGB', 'rssEndGB', 'rssPeakGB', 'rssPeakBeforeGB', 'readBytes', 'writeBytes', 'rchar', 'wchar', ...
    'minorFaults', 'majorFaults', 'threads', 'cpus', 'slurmCpus', 'osThreads', 'gpu', 'gpuMemGB', 'gpuFreeGB', ...
    'cudaDevices', 'processStart', 'memTotalGB', 'batchSize', 'parent', 'metricsSec'};

isTask = cellfun(@(x) isstruct(x) && isfield(x, 'kind') && strcmp(x.kind, 'task'), recs);
recs = recs(isTask);
n = numel(recs);
s = struct();
for i = 1 : numel(strNames)
    s.(strNames{i}) = strings(n, 1);
end
for i = 1 : numel(numNames)
    s.(numNames{i}) = NaN(n, 1);
end
for k = 1 : n
    r = recs{k};
    for i = 1 : numel(strNames)
        if isfield(r, strNames{i}) && ischar(r.(strNames{i}))
            s.(strNames{i})(k) = r.(strNames{i});
        end
    end
    for i = 1 : numel(numNames)
        if isfield(r, numNames{i}) && isnumeric(r.(numNames{i})) && ~isempty(r.(numNames{i}))
            s.(numNames{i})(k) = r.(numNames{i})(1);
        end
    end
end
t = struct2table(s);
t = t(:, order);

end


function [t] = add_startup(t)
% process start to its first task, for tasks in a process of their own

t.startupSec = NaN(height(t), 1);
own = ismember(t.mode, ["slurm_matlab", "slurm_mcc", "mcc"]) & ~isnan(t.processStart);
keys = t.node + "_" + string(t.pid) + "_" + string(round(t.processStart));
for key = unique(keys(own))'
    inds = find(own & keys == key);
    [~, k] = min(t.start(inds));
    t.startupSec(inds(k)) = t.start(inds(k)) - t.processStart(inds(k));
end

end


function [t] = add_joblog(t, lines)
% GNU parallel's line run time and exit for compiled tasks; lines without a
% record get a row

t.lineSec = NaN(height(t), 1);
t.exitSec = NaN(height(t), 1);
extra = t([], :);
for i = 1 : numel(lines)
    l = lines{i};
    lineEnd = l.lineStart + l.lineSec;
    k = find(t.runId == l.runId & t.task == l.task & (t.jobId == l.jobId | l.jobId == "") ...
        & t.start >= l.lineStart - 1 & t.end <= lineEnd + 1 & isnan(t.lineSec), 1);
    if ~isempty(k)
        t.lineSec(k) = l.lineSec;
        t.exitSec(k) = lineEnd - t.end(k);
        continue;
    end
    row = task_table({struct('kind', 'task')});
    row.startupSec = NaN;
    row.lineSec = NaN;
    row.exitSec = NaN;
    row.task(1) = l.task;
    row.runId(1) = l.runId;
    row.jobId(1) = l.jobId;
    row.start(1) = l.lineStart;
    row.end(1) = lineEnd;
    row.lineSec(1) = l.lineSec;
    row.exitCode(1) = l.exitCode;
    if l.signal > 0
        row.status(1) = "killed";
    elseif l.exitCode ~= 0
        row.status(1) = "failed";
    else
        row.status(1) = "skipped";
    end
    extra = [extra; row]; %#ok<AGROW>
end
t = [t; extra];

end


function [t, runTable] = add_run_settings(t, runs)
% the task's input, output and requested resources from its framework call
% (or single-job submission), and a row per run

n = height(t);
t.input = strings(n, 1);
t.output = strings(n, 1);
t.cpusPerTask = NaN(n, 1);
t.memAllocate = NaN(n, 1);

nr = numel(runs);
rt = struct('runId', strings(nr, 1), 'kind', strings(nr, 1), 'func', strings(nr, 1), ...
    'backend', strings(nr, 1), 'node', strings(nr, 1), 'start', NaN(nr, 1), 'end', NaN(nr, 1), ...
    'elapsedSec', NaN(nr, 1), 'nTasks', NaN(nr, 1), 'nDone', NaN(nr, 1), 'nRecords', NaN(nr, 1), ...
    'parseCluster', NaN(nr, 1), 'mccMode', NaN(nr, 1), 'cpusPerTask', NaN(nr, 1), 'memAllocate', NaN(nr, 1), ...
    'taskBatchNum', NaN(nr, 1), 'paraJobNum', NaN(nr, 1), 'GPUJob', NaN(nr, 1), 'parent', strings(nr, 1), ...
    'folder', strings(nr, 1));
for i = 1 : nr
    rt.runId(i) = runs{i}.runId;
    rt.folder(i) = runs{i}.folder;
    rt.nRecords(i) = runs{i}.nRecords;
    rt.kind(i) = "framework";
    inRun = t.runId == runs{i}.runId;
    recs = runs{i}.records;
    for j = 1 : numel(recs)
        r = recs{j};
        if ~isstruct(r) || ~isfield(r, 'kind')
            continue;
        end
        switch r.kind
            case 'run_start'
                rt.node(i) = get_str(r, 'node');
                rt.start(i) = get_num(r, 'time') - max(0, get_num(r, 'elapsedSec'), 'omitnan');
                rt.backend(i) = get_str(r, 'backend');
                rt.parent(i) = get_str(r, 'parent');
                names = {'parseCluster', 'mccMode', 'cpusPerTask', 'memAllocate', 'taskBatchNum', 'paraJobNum', 'GPUJob'};
                for k = 1 : numel(names)
                    rt.(names{k})(i) = get_num(r, names{k});
                end
                inputs = get_cell(r, 'inputFullpaths');
                outputs = get_cell(r, 'outputFullpaths');
                funcs = get_cell(r, 'functionStrs');
                rt.nTasks(i) = numel(outputs);
                if ~isempty(funcs)
                    rt.func(i) = strtok(first_str(funcs{1}), '(');
                end
                % tasks of the run: their input and output and the request
                rows = find(inRun);
                for k = rows(:)'
                    f = t.task(k);
                    if f >= 1 && f <= numel(outputs)
                        t.output(k) = first_str(outputs{f});
                    end
                    if f >= 1 && f <= numel(inputs)
                        t.input(k) = first_str(inputs{f});
                    end
                end
                t.cpusPerTask(inRun) = rt.cpusPerTask(i);
                t.memAllocate(inRun) = rt.memAllocate(i);
                % compiled tasks in the master's Slurm job (or, like the
                % master, in none) are the ones the master ran itself
                byMaster = inRun & (t.mode == "slurm_mcc" | t.mode == "mcc") & t.jobId == get_str(r, 'jobId');
                if rt.parseCluster(i) == 0
                    t.mode(byMaster) = "local";
                else
                    t.mode(byMaster) = "master";
                end
            case 'run_end'
                rt.end(i) = get_num(r, 'time');
                rt.elapsedSec(i) = get_num(r, 'elapsedSec');
                rt.nTasks(i) = get_num(r, 'nTasks');
                rt.nDone(i) = get_num(r, 'nDone');
            case 'submit'
                % a single job: its task by the job and array task
                rt.kind(i) = "single";
                rt.node(i) = get_str(r, 'node');
                rt.start(i) = min(rt.start(i), get_num(r, 'time'), 'omitnan');
                rt.nTasks(i) = sum([rt.nTasks(i), 1], 'omitnan');
                rows = inRun & t.arrayJobId == string(sprintf('%d', get_num(r, 'jobId'))) ...
                    & t.task == get_num(r, 'task');
                t.cpusPerTask(rows) = get_num(r, 'cpusPerTask');
                t.memAllocate(rows) = get_num(r, 'memAllocate');
                if ~any(rows)
                    continue;
                end
                funcStr = get_str(r, 'functionStr');
                if strlength(funcStr) > 0
                    rt.func(i) = strtok(funcStr, '(');
                end
        end
    end
end
runTable = struct2table(rt);

end


function [t] = add_slurm_accounting(t)
% each task's Slurm job from sacct: submit, start and end (seconds since
% 1970), state, queue wait, and the job's time before its first and after
% its last task

n = height(t);
t.jobSubmit = NaN(n, 1);
t.jobStart = NaN(n, 1);
t.jobEnd = NaN(n, 1);
t.jobState = strings(n, 1);
t.queueSec = NaN(n, 1);
t.jobHeadSec = NaN(n, 1);
t.jobTailSec = NaN(n, 1);
ids = unique(t.jobId(t.jobId ~= ""));
fmt = 'yyyy-MM-dd''T''HH:mm:ss';
for b = 1 : 500 : numel(ids)
    bids = ids(b : min(b + 499, numel(ids)));
    [status, out] = system(sprintf('sacct -n -P -X --format=JobIDRaw,Submit,Start,End,State -j %s', ...
        strjoin(bids, ',')));
    if status ~= 0
        warning('read_task_metrics:sacct', 'sacct failed: %s', strtrim(out));
        break;
    end
    for l = splitlines(strtrim(out))'
        c = strsplit(l{1}, '|');
        if numel(c) < 5
            continue;
        end
        rows = t.jobId == c{1};
        times = NaN(1, 3);
        for k = 1 : 3
            try
                times(k) = posixtime(datetime(c{k + 1}, 'InputFormat', fmt, 'TimeZone', 'local'));
            catch
            end
        end
        t.jobSubmit(rows) = times(1);
        t.jobStart(rows) = times(2);
        t.jobEnd(rows) = times(3);
        t.jobState(rows) = strtok(c{5});
        t.queueSec(rows) = times(2) - times(1);
    end
end
% the job's start to its first task and its last task to its end (sacct
% times are whole seconds)
for id = ids(:)'
    rows = t.jobId == id;
    t.jobHeadSec(rows) = min(t.start(rows)) - t.jobStart(find(rows, 1));
    t.jobTailSec(rows) = t.jobEnd(find(rows, 1)) - max(t.end(rows));
end

end


function [s] = get_str(r, name)

s = "";
if isfield(r, name) && ischar(r.(name))
    s = string(r.(name));
end

end


function [x] = get_num(r, name)

x = NaN;
if isfield(r, name) && (isnumeric(r.(name)) || islogical(r.(name))) && ~isempty(r.(name))
    x = double(r.(name)(1));
end

end


function [c] = get_cell(r, name)
% a list of the run record (jsondecode gives a char for one element)

c = {};
if isfield(r, name)
    c = r.(name);
    if ~iscell(c)
        c = {c};
    end
end

end


function [s] = first_str(x)
% the first path of an input (inputs can be lists)

while iscell(x) && ~isempty(x)
    x = x{1};
end
s = "";
if ischar(x) || isstring(x)
    s = string(x);
end

end
