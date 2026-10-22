function [is_done_flag] = matlab_pool_generic_computing_wrapper(inputFullpaths, outputFullpaths, funcStrs, varargin)
% run the tasks (matlab function strings) concurrently on the local MATLAB
% process pool (Parallel Computing Toolbox), at most workerNum at a time. Uses
% the current process pool, or starts one with workerNum workers that share
% the computational threads of this MATLAB. Each task gets the threads divided
% by the number of tasks running at once. Tasks whose output exists or whose
% input is missing are skipped, and no task runs once the final output
% exists. A task that fails is reported and left unfinished, for the caller to
% retry. Without a usable pool (no toolbox or license, a thread-based pool,
% compiled mode) no task runs. The pool is kept for later calls.
%
% Parameters (as 'specifier'-value pairs):
%   finalOutFullpath : char (default: ''). Stop once this path exists (another run finished the work).
%          workerNum : a number (default: 2). Maximum number of tasks running at once.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('outputFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('functionStrs', @(x) iscell(x) || ischar(x));
ip.addParameter('finalOutFullpath', '', @ischar);
ip.addParameter('workerNum', 2, @isnumeric);

ip.parse(inputFullpaths, outputFullpaths, funcStrs, varargin{:});

pr = ip.Results;
finalOutFullpath = pr.finalOutFullpath;
workerNum = max(1, round(pr.workerNum));

if ischar(outputFullpaths)
    inputFullpaths = {inputFullpaths};
    outputFullpaths = {outputFullpaths};
    funcStrs = {funcStrs};
end
nF = numel(outputFullpaths);
if ~isempty(finalOutFullpath) && batch_file_exist({finalOutFullpath})
    is_done_flag = true(nF, 1);
    return;
end
is_done_flag = batch_file_exist(outputFullpaths, [], true);
is_done_flag = is_done_flag(:);
if all(is_done_flag)
    return;
end
% tasks with missing inputs are skipped, as the framework does
if iscellstr(inputFullpaths)
    input_exist_mat = batch_file_exist(inputFullpaths, [], true);
else
    input_exist_mat = cellfun(@(x) all(cellfun(@(y) exist(y, 'file') > 0, cellstr(x))), inputFullpaths);
end
task_inds = find(~is_done_flag & input_exist_mat(:));
nT = numel(task_inds);
if nT == 0
    return;
end

p = get_local_process_pool(workerNum, max(1, floor(maxNumCompThreads / workerNum)));
if isempty(p)
    return;
end
workerNum = min([workerNum, p.NumWorkers, nT]);
threadNum = max(1, floor(maxNumCompThreads / workerNum));
fprintf('Run %d tasks on the local pool (%d at a time)...\n', nT, workerNum);
ts = tic;
fs = parallel.FevalFuture.empty;
run_inds = zeros(0, 1);
next = 1;
nF_done = sum(is_done_flag);
while next <= nT || ~isempty(fs)
    % keep up to workerNum tasks running, in task order
    while numel(fs) < workerNum && next <= nT
        f = task_inds(next);
        fs(end + 1) = parfeval(p, @run_task_str, 2, funcStrs{f}, threadNum);
        run_inds(end + 1, 1) = f;
        next = next + 1;
    end

    try
        [k, t1, msg] = fetchNext(fs, 1);
    catch ME
        % the worker failed (e.g., crashed); the task stays unfinished
        k = find(arrayfun(@(x) strcmp(x.State, 'finished') && ~isempty(x.Error), fs), 1, 'first');
        if isempty(k)
            % the pool can't run the tasks: leave them to the caller
            disp(ME);
            cancel(fs);
            break;
        end
        t1 = NaN;
        msg = getReport(ME, 'extended', 'hyperlinks', 'off');
    end
    if isempty(k)
        continue;
    end
    f = run_inds(k);
    fprintf('%s', fs(k).Diary);
    fs(k) = [];
    run_inds(k) = [];

    if ~isempty(msg)
        fprintf('Task %d failed: %s\n', f, msg);
    end
    is_done_flag(f) = batch_file_exist(outputFullpaths(f), [], true);
    fprintf('Done! Elapsed time is %f seconds.\n', t1);
    if nF_done < sum(is_done_flag)
        nF_done = sum(is_done_flag);
        fprintf('Time %0.2f s: %d / %d (%0.2f%%) are finished!\n', toc(ts), nF_done, nF, nF_done / nF * 100);
    end

    % another run made the final output: stop, as the framework does
    if ~isempty(finalOutFullpath) && batch_file_exist({finalOutFullpath})
        cancel(fs);
        is_done_flag = true(nF, 1);
        fprintf('Time %0.2f s: The final output file %s already exists!\n', toc(ts), finalOutFullpath);
        break;
    end
end

end


function [p] = get_local_process_pool(workerNum, threadNum)
% the current process pool, or a new one with workerNum workers of threadNum
% threads; empty if none can be used

p = [];
if isdeployed || ~license('test', 'Distrib_Computing_Toolbox') || isempty(which('gcp'))
    return;
end
try
    p = gcp('nocreate');
    if ~isempty(p)
        % thread-based workers cannot run mex functions
        if isa(p, 'parallel.ThreadPool') || isa(p, 'parallel.BackgroundPool')
            p = [];
        end
        return;
    end
    try
        c = parcluster('Processes');
    catch
        c = parcluster('local');
    end
    c.NumWorkers = max(c.NumWorkers, workerNum);
    c.NumThreads = threadNum;
    % job files of its own: MATLAB sessions that start pools at the same time
    % with the default folder in a shared (NFS) home can corrupt each other's
    jobStorageLocation = tempname;
    mkdir(jobStorageLocation);
    c.JobStorageLocation = jobStorageLocation;
    p = parpool(c, workerNum);
catch ME
    disp(ME);
    p = [];
end

end


function [t, msg] = run_task_str(func_str, threadNum)
% run a task string on a worker with threadNum computational threads; errors
% are returned as text

msg = '';
t0 = tic;
nt = maxNumCompThreads(threadNum);
try
    feval(str2func(['@()', func_str]));
catch ME
    msg = getReport(ME, 'extended', 'hyperlinks', 'off');
end
maxNumCompThreads(nt);
t = toc(t0);

end
