function [info] = isolation_reset(level, st0)
% reset the process between two tasks toward st0, an isolation_state snapshot
% taken before the first task, so that the next task finds what a new process
% would give it.
%
%   level : 'none' | 'prototype' | 'full' | 'full_gpu', or a function handle
%           called as fcn(st0).
%           'prototype': what the 2026-10-22 task runner prototype restored:
%             the default random numbers, the warning states, the current
%             folder, and closing the files the task left open.
%           'full': 'prototype' and lastwarn, persistent variables and
%             unlocked mex files (clear functions), global variables and
%             variables in the base workspace, the path, the environment,
%             the number of computational threads, figures, timers and a
%             parallel pool the task started.
%           'full_gpu': 'full' and no GPU selected (gpuDevice([]), which
%             frees its memory) when a task selected one.
%     st0 : struct. isolation_state snapshot of the process before its tasks.
%
% Output: struct with fields time (s) and actions (cell array of what was reset).


info.actions = {};
t0 = tic;
if isa(level, 'function_handle')
    level(st0);
    info.actions = {'custom'};
    info.time = toc(t0);
    return;
end
if strcmp(level, 'none')
    info.time = toc(t0);
    return;
end
if ~any(strcmp(level, {'prototype', 'full', 'full_gpu'}))
    error('Unknown reset level %s.', level);
end

% the prototype's reset
rng('default');
warning(st0.warnings);
if ~strcmp(pwd, st0.pwd)
    cd(st0.pwd);
    info.actions{end + 1} = 'cd';
end
fids = setdiff(isolation_open_fids(), [st0.files{:, 1}]);
for i = 1 : numel(fids)
    fclose(fids(i));
end
if ~isempty(fids)
    info.actions{end + 1} = sprintf('fclose %d', numel(fids));
end
if strcmp(level, 'prototype')
    info.time = toc(t0);
    return;
end

lastwarn('', '');
clear('functions');
clearvars('-global');
names = setdiff(evalin('base', 'who'), st0.baseVars);
if ~isempty(names)
    evalin('base', ['clear ', strjoin(names, ' ')]);
    info.actions{end + 1} = 'base workspace';
end
if ~strcmp(path, st0.path)
    path(st0.path);
    info.actions{end + 1} = 'path';
end
restore_env(st0.env);
if maxNumCompThreads ~= st0.threads
    maxNumCompThreads(st0.threads);
    info.actions{end + 1} = 'threads';
end
try
    figs = findall(groot, 'Type', 'figure');
    if numel(figs) > st0.figures
        close(figs, 'force');
        info.actions{end + 1} = 'figures';
    end
catch
end
try
    t = timerfindall;
    if numel(t) > st0.timers
        stop(t);
        delete(t);
        info.actions{end + 1} = 'timers';
    end
catch
end
try
    if ~st0.pool && ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'));
        info.actions{end + 1} = 'pool';
    end
catch
end
if strcmp(level, 'full_gpu')
    % no device selected, as in a new process (its memory is freed)
    try
        if parallel.internal.gpu.isAnyDeviceSelected
            gpuDevice([]);
            info.actions{end + 1} = 'gpu';
        end
    catch
    end
end
info.time = toc(t0);

end


function [] = restore_env(env0)
% set the environment back to env0 ({name, value} rows)

e = getenv();
added = setdiff(cellstr(e.keys), env0(:, 1));
for i = 1 : numel(added)
    if exist('unsetenv', 'builtin') || exist('unsetenv', 'file')
        unsetenv(added{i});
    else
        setenv(added{i}, '');
    end
end
for i = 1 : size(env0, 1)
    v = getenv(env0{i, 1});
    if ~strcmp(v, env0{i, 2})
        setenv(env0{i, 1}, env0{i, 2});
    end
end

end
