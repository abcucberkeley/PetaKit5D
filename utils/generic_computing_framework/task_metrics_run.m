function [out] = task_metrics_run(action, varargin)
% sets up the task metrics of a job framework call (taskMetrics option): the
% metrics folder <jobLogDir>/task_metrics/<runId>/, the master's records of
% the call in run.ndjson there (settings and tasks at the start, the result
% at the end, single job submissions), and the metrics parts of job commands
% and compiled task lines. The tasks record themselves with task_metrics;
% read_task_metrics collects both. Errors are ignored: no metrics rather
% than a failed run.
%
% Usage:
%   runDir = task_metrics_run('run_dir', jobLogDir, inputFullpaths, funcStrs)
%       a new metrics folder for a framework call: in jobLogDir if it is an
%       existing absolute path, otherwise in <first input folder>/job_logs;
%       '' if it can't be made
%   runDir = task_metrics_run('backend_dir', taskMetrics, jobLogDir, inputFullpaths, outputFullpaths, funcStrs, backend)
%       the metrics folder of a framework backend: taskMetrics itself if it
%       is a folder (from generic_computing_frameworks_wrapper), a new one
%       with a record of the tasks if it is true, else ''
%   runDir = task_metrics_run('session_dir', jobLogDir, taskMetrics)
%       the metrics folder of the single jobs this MATLAB submits to jobLogDir
%       ('' if taskMetrics is false)
%   task_metrics_run('record', runDir, kind, info)
%       appends the struct info as a record of the given kind to run.ndjson
%   str = task_metrics_run('matlab_task', runDir, taskInd, funcStr, GPUJob)
%       the task statement with its metrics calls, for a MATLAB job's command
%   str = task_metrics_run('env', runDir, taskInd, GPUJob)
%       the metrics environment of a compiled task line (sh, or cmd on Windows)
%
% Without metrics, a framework call run by a task with metrics clears the
% task's metrics environment (PETAKIT5D_TASK_METRICS, PETAKIT5D_TASK,
% PETAKIT5D_TASK_GPU), so that the tasks it starts don't record into the
% task's run.


out = [];
try
    switch action
        case 'run_dir'
            out = run_dir(varargin{:});
        case 'backend_dir'
            out = backend_dir(varargin{:});
        case 'session_dir'
            out = session_dir(varargin{:});
        case 'record'
            run_record(varargin{:});
        case 'matlab_task'
            out = matlab_task(varargin{:});
        case 'env'
            out = task_env(varargin{:});
        otherwise
            error('task_metrics_run:unknownAction', 'Unknown action %s.', action);
    end
catch ME
    % no metrics rather than a failed run: a MATLAB job runs the plain task
    % statement, a compiled task line gets no metrics environment
    switch action
        case 'matlab_task'
            out = varargin{3};
        case {'run_dir', 'backend_dir', 'session_dir', 'env'}
            out = '';
    end
    warning('task_metrics_run:failed', 'Task metrics (%s) are incomplete: %s', action, ME.message);
end

end


function [runDir] = run_dir(jobLogDir, inputFullpaths, funcStrs)
% a new folder <base>/task_metrics/<date>_<time>_<function>_<id>/

firstInput = inputFullpaths;
while iscell(firstInput) && ~isempty(firstInput)
    firstInput = firstInput{1};
end
dataPath = strrep(fileparts(firstInput), '\', '/');
if ~is_absolute(dataPath)
    dataPath = [strrep(pwd, '\', '/'), '/', dataPath];
end
if ischar(jobLogDir) && is_absolute(jobLogDir) && exist(jobLogDir, 'dir')
    base = strrep(jobLogDir, '\', '/');
else
    base = [dataPath, '/job_logs'];
end
funcStr = funcStrs;
while iscell(funcStr) && ~isempty(funcStr)
    funcStr = funcStr{1};
end
name = safe_name(function_name(funcStr));
runDir = new_dir(base, name(1 : min(end, 40)));

end


function [runDir] = backend_dir(taskMetrics, jobLogDir, inputFullpaths, outputFullpaths, funcStrs, backend)
% the metrics folder of a backend; when called directly, a new one with a
% record of the backend's tasks

runDir = '';
if ischar(taskMetrics) && ~isempty(taskMetrics)
    runDir = [taskMetrics(1 : find(taskMetrics ~= '/', 1, 'last')), '/'];
elseif ~ischar(taskMetrics) && ~isempty(taskMetrics) && taskMetrics(1)
    runDir = run_dir(jobLogDir, inputFullpaths, funcStrs);
    run_record(runDir, 'run_start', struct('backend', backend, 'jobLogDir', jobLogDir, ...
        'inputFullpaths', {inputFullpaths}, 'outputFullpaths', {outputFullpaths}, 'functionStrs', {funcStrs}));
else
    clear_task_env();
end

end


function [] = clear_task_env()
% the metrics environment a task got from its task line, if any, is not for
% the tasks of a framework call without metrics

if ~isempty(getenv('PETAKIT5D_TASK_METRICS'))
    setenv('PETAKIT5D_TASK_METRICS', '');
    setenv('PETAKIT5D_TASK', '');
    setenv('PETAKIT5D_TASK_GPU', '');
end

end


function [runDir] = session_dir(jobLogDir, taskMetrics)
% one folder per MATLAB session and job log folder for single jobs

persistent bases dirs
if nargin > 1 && ~taskMetrics
    runDir = '';
    clear_task_env();
    return;
end
if isempty(bases)
    bases = {};
    dirs = {};
end
k = find(strcmp(bases, jobLogDir), 1);
if ~isempty(k) && exist(dirs{k}, 'dir')
    runDir = dirs{k};
    return;
end
runDir = '';
base = strrep(jobLogDir, '\', '/');
if ~is_absolute(base)
    base = [strrep(pwd, '\', '/'), '/', base];
end
if ~exist(base, 'dir')
    return;
end
runDir = new_dir(base, 'single');
if isempty(runDir)
    return;
end
if isempty(k)
    k = numel(bases) + 1;
end
bases{k} = jobLogDir;
dirs{k} = runDir;

end


function [runDir] = new_dir(base, name)
% a new folder <base>/task_metrics/<date>_<time>_<name>_<id>/, '' if it can't be made

[~, token] = fileparts(tempname);
runId = sprintf('%s_%s_%s', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')), name, token(end - 7 : end));
runDir = sprintf('%s/task_metrics/%s/', base(1 : find(base ~= '/', 1, 'last')), runId);
if ~mkdir(runDir)
    runDir = '';
end

end


function [] = run_record(runDir, kind, info)
% appends a record of the run (master only) to run.ndjson

if isempty(runDir)
    return;
end
p = task_metrics('process');
d = runDir(1 : find(runDir ~= '/', 1, 'last'));
% the master's Slurm job tells the tasks it runs itself from those of jobs
r = struct('kind', kind, 'version', 1, 'runId', d(find(d == '/', 1, 'last') + 1 : end), ...
    'time', posixtime(datetime('now', 'TimeZone', 'local')), 'node', p.host, 'pid', p.pid, ...
    'jobId', getenv('SLURM_JOB_ID'), 'parent', p.current);
fns = fieldnames(info);
for i = 1 : numel(fns)
    r.(fns{i}) = info.(fns{i});
end
fid = fopen([d, '/run.ndjson'], 'a');
if fid < 0
    error('task_metrics_run:write', 'Cannot write %s/run.ndjson.', d);
end
fprintf(fid, '%s\n', jsonencode(r));
fclose(fid);

end


function [str] = matlab_task(runDir, taskInd, funcStr, GPUJob)
% the task statement between its metrics calls; the calls are in try blocks
% of their own, so that the task also runs without task_metrics on the path

q = @(x) strrep(x, '''', '''''');
str = sprintf(['try;pkm_=task_metrics(''start'',''%s'',''%s'',%d,''slurm_matlab'',%d);catch;pkm_=[];end;', ...
    'try;%s;catch pke_;try;task_metrics(''stop'',pkm_,pke_);catch;end;rethrow(pke_);end;', ...
    'try;task_metrics(''stop'',pkm_);catch;end'], q(function_name(funcStr)), q(runDir), ...
    taskInd, GPUJob, funcStr);

end


function [str] = task_env(runDir, taskInd, GPUJob)
% the environment of a compiled task line

if nargin < 3
    GPUJob = false;
end
if isunix
    str = sprintf('PETAKIT5D_TASK_METRICS="%s" PETAKIT5D_TASK=%d ', runDir, taskInd);
    if GPUJob
        str = [str, 'PETAKIT5D_TASK_GPU=1 '];
    end
else
    str = sprintf('set "PETAKIT5D_TASK_METRICS=%s" & set "PETAKIT5D_TASK=%d" & ', runDir, taskInd);
    if GPUJob
        str = [str, 'set "PETAKIT5D_TASK_GPU=1" & '];
    end
end

end


function [name] = function_name(funcStr)
% the function a task string calls

k = find(funcStr == '(', 1);
if isempty(k)
    name = strtrim(funcStr);
else
    name = strtrim(funcStr(1 : k - 1));
end

end


function [s] = safe_name(s)
% for file names: letters, digits, '_', '-' and '.'

s(~((s >= 'a' & s <= 'z') | (s >= 'A' & s <= 'Z') | (s >= '0' & s <= '9') | s == '_' | s == '-' | s == '.')) = '_';

end


function [tf] = is_absolute(p)

tf = ~isempty(p) && (p(1) == '/' || p(1) == '\' || (numel(p) > 1 && p(2) == ':'));

end
