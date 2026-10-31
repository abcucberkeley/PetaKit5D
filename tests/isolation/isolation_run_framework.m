function [is_done_flag] = isolation_run_framework(tasks, varargin)
% run tasks (see isolation_prepare) through the job framework
% (generic_computing_frameworks_wrapper) in one call, in the given order, so
% that a framework that runs several tasks per process (compiled workers,
% MATLAB jobs, local runs) can be checked against fresh processes with
% isolation_compare. The framework options decide where and how the tasks
% run, e.g. {'parseCluster', false} (this MATLAB) or {'parseCluster', true,
% 'mccMode', true, 'configFile', '<cluster config>'} (compiled jobs).
%
% Inputs:
%   tasks : struct array (isolation_prepare) or the .mat file that holds it.
%
% Parameters (as 'specifier'-value pairs):
%        order : vector (default: all, in order). Indices of the tasks to run, in this order.
%   frameworkArgs : cell array (default: {'parseCluster', false}). Name-value arguments of generic_computing_frameworks_wrapper.
%      noSlurm : true|false (default: true). Hide Slurm's commands (false to run the tasks as Slurm jobs).
%
% Output: is_done_flag of the framework call, in the order of the tasks run.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tasks', @(x) isstruct(x) || ischar(x));
ip.addParameter('order', [], @isnumeric);
ip.addParameter('frameworkArgs', {'parseCluster', false}, @iscell);
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
tasks = tasks(order);

% the framework skips a task whose input is missing; a task without an
% output it can check gets a flag file, written by the task string after the
% task (MATLAB mode only: compiled tasks are one function call)
flagDir = [tempname, '/'];
mkdir(flagDir);
funcStrs = {tasks.funcStr}';
outputs = {tasks.output}';
inputs = repmat({flagDir}, numel(tasks), 1);
if isfield(tasks, 'input')
    hasInput = ~cellfun(@isempty, {tasks.input});
    inputs(hasInput) = {tasks(hasInput).input};
end
for i = 1 : numel(tasks)
    if isempty(outputs{i})
        outputs{i} = sprintf('%s%s.done', flagDir, tasks(i).id);
        funcStrs{i} = sprintf('%s;fclose(fopen(''%s'',''w''))', funcStrs{i}, outputs{i});
    end
end
if pr.noSlurm && ~ispc
    origPATH = getenv('PATH');
    stubDir = [flagDir, 'slurm_stubs'];
    mkdir(stubDir);
    for cmd = {'sinfo', 'sbatch', 'squeue', 'scancel'}
        fid = fopen([stubDir, '/', cmd{1}], 'w');
        fprintf(fid, '#!/bin/sh\nexit 1\n');
        fclose(fid);
        fileattrib([stubDir, '/', cmd{1}], '+x');
    end
    setenv('PATH', [stubDir, pathsep, origPATH]);
    cleanupObj = onCleanup(@() setenv('PATH', origPATH));
end
is_done_flag = generic_computing_frameworks_wrapper(inputs, outputs, funcStrs, pr.frameworkArgs{:});

end
