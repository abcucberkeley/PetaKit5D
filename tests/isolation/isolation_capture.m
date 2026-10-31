function [calls] = isolation_capture(fcn, varargin)
% run fcn (a data wrapper call) with the job framework replaced by a recorder,
% and return the tasks it hands to the framework: the state of the data
% folder when the call stops is the state just before those tasks run.
%
% The recorder replaces generic_computing_frameworks_wrapper: calls before
% stopAtCall run their tasks here, one after another (as the framework does
% without a cluster), and call stopAtCall is recorded and stops fcn. With
% singleJob, generic_single_job_submit_wrapper (one job per task, e.g. the
% decon wrapper in cluster mode) is recorded instead of run, checkSlurmCluster
% reports a cluster, and check_slurm_job_status reports the jobs as gone, so
% the wrapper hands out its tasks as it would on the cluster and returns.
% Slurm's commands fail while fcn runs, so nothing is ever submitted.
%
% Parameters (as 'specifier'-value pairs):
%   stopAtCall : a number (default: 1). The framework call that is recorded
%                and stops fcn (inf: run all calls).
%    singleJob : true|false (default: false). Record single job submissions
%                (see above).
%
% Output: cell array of structs with fields inputFullpaths, outputFullpaths,
%   funcStrs, options (name-value arguments) and type ('framework' or
%   'singleJob'), in call order.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fcn', @(x) isa(x, 'function_handle'));
ip.addParameter('stopAtCall', 1, @isnumeric);
ip.addParameter('singleJob', false, @islogical);
ip.parse(fcn, varargin{:});

pr = ip.Results;

shadowDir = strrep(tempname, '\', '/');
mkdir(shadowDir);
cfg.stopAtCall = pr.stopAtCall;
cfg.singleJob = pr.singleJob;
save([shadowDir, '/config.mat'], 'cfg');

write_lines([shadowDir, '/generic_computing_frameworks_wrapper.m'], {
    'function [is_done_flag] = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, funcStrs, varargin)'
    '% recorder written by isolation_capture'
    'd = fileparts(mfilename(''fullpath''));'
    'a = load([d, ''/config.mat'']);'
    'calls = {};'
    'if exist([d, ''/calls.mat''], ''file'')'
    '    b = load([d, ''/calls.mat'']);'
    '    calls = b.calls;'
    'end'
    'if ischar(funcStrs)'
    '    funcStrs = {funcStrs};'
    'end'
    'if ischar(inputFullpaths)'
    '    inputFullpaths = {inputFullpaths};'
    'end'
    'if ischar(outputFullpaths)'
    '    outputFullpaths = {outputFullpaths};'
    'end'
    'calls{end + 1} = struct(''type'', ''framework'', ''inputFullpaths'', {inputFullpaths}, ''outputFullpaths'', {outputFullpaths}, ...'
    '    ''funcStrs'', {funcStrs}, ''options'', {varargin});'
    'save([d, ''/calls.mat''], ''calls'');'
    'nFramework = sum(cellfun(@(x) strcmp(x.type, ''framework''), calls));'
    'if nFramework >= a.cfg.stopAtCall'
    '    error(''isolation:stop'', ''isolation_capture: stopped at framework call %d'', nFramework);'
    'end'
    'for i = 1 : numel(funcStrs)'
    '    if exist(outputFullpaths{i}, ''file'') || exist(outputFullpaths{i}, ''dir'')'
    '        continue;'
    '    end'
    '    feval(str2func([''@()'', funcStrs{i}]));'
    'end'
    'is_done_flag = true(numel(funcStrs), 1);'
    'end'});

if pr.singleJob
    write_lines([shadowDir, '/generic_single_job_submit_wrapper.m'], {
        'function [job_id, job_status, submit_status] = generic_single_job_submit_wrapper(funcStr, job_id, task_id, varargin)'
        '% recorder written by isolation_capture'
        'd = fileparts(mfilename(''fullpath''));'
        'calls = {};'
        'if exist([d, ''/calls.mat''], ''file'')'
        '    b = load([d, ''/calls.mat'']);'
        '    calls = b.calls;'
        'end'
        'calls{end + 1} = struct(''type'', ''singleJob'', ''inputFullpaths'', {{}}, ''outputFullpaths'', {{}}, ...'
        '    ''funcStrs'', {{funcStr}}, ''options'', {varargin});'
        'save([d, ''/calls.mat''], ''calls'');'
        'job_id = numel(calls);'
        'job_status = -1;'
        'submit_status = 1;'
        'end'});
    write_lines([shadowDir, '/checkSlurmCluster.m'], {
        'function [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str, jobLogDir] = checkSlurmCluster(dataPath, jobLogDir)'
        '% recorder written by isolation_capture'
        'parseCluster = true;'
        'if nargin < 2 || isempty(jobLogDir)'
        '    jobLogDir = tempdir;'
        'end'
        'job_log_fname = [jobLogDir, ''/job_%A_%a.out''];'
        'job_log_error_fname = [jobLogDir, ''/job_%A_%a.err''];'
        'slurm_constraint_str = '''';'
        'end'});
    write_lines([shadowDir, '/check_slurm_job_status.m'], {
        'function [job_status] = check_slurm_job_status(job_id, array_id)'
        '% recorder written by isolation_capture'
        'job_status = -1;'
        'end'});
end

% Slurm's commands fail, so that nothing is submitted
origPATH = getenv('PATH');
if ~ispc
    for cmd = {'sinfo', 'sbatch', 'squeue', 'scancel'}
        write_lines([shadowDir, '/', cmd{1}], {'#!/bin/sh', 'exit 1'});
        fileattrib([shadowDir, '/', cmd{1}], '+x');
    end
    setenv('PATH', [shadowDir, pathsep, origPATH]);
end
addpath(shadowDir, '-begin');
cleanupObj = onCleanup(@() remove_shadow(shadowDir, origPATH));

try
    fcn();
catch ME
    if ~strcmp(ME.identifier, 'isolation:stop')
        rethrow(ME);
    end
end

calls = {};
if exist([shadowDir, '/calls.mat'], 'file')
    a = load([shadowDir, '/calls.mat']);
    calls = a.calls;
end

end


function [] = remove_shadow(shadowDir, origPATH)

rmpath(shadowDir);
setenv('PATH', origPATH);
clear generic_computing_frameworks_wrapper generic_single_job_submit_wrapper checkSlurmCluster check_slurm_job_status
if exist(shadowDir, 'dir')
    rmdir(shadowDir, 's');
end

end


function [] = write_lines(fn, lines)

fid = fopen(fn, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);

end
