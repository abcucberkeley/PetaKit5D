function [is_done_flag] = slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, funcStrs, varargin)
% generic computing framework that use a function handle/string as input
% for the computing
% 
% 
% Author: Xiongtao Ruan (10/30/2020)
% xruan (12/27/2020): add support for a batch of input filenames for each instance
% xruan (01/20/2021): add support to wait several seconds for just finished
% jobs, in case that the results files are instantly recognized by the master job. 
% xruan (06/24/2021): add support to retry matlab job by saving func str to
% the disk and load it in matlab when func str in wrap fails.
% xruan (07/15/2021): improving the transition of a job finished while the
% result is not yet visible to the master job.
% xruan(08/26/2021): add support for limiting active jobs (pending/running). 
%   also add support for batch tasks running, i.e., run several tasks within one job
% xruan(02/06/2021): change to first check all results at once instead of once
%   per file to reduce the IO requirement. 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('outputFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('functionStrs', @(x) iscell(x) || ischar(x));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('tmpDir', '', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);
ip.addParameter('maxJobNum', inf, @isnumeric); % submit limited number of jobs (pending/running)
ip.addParameter('taskBatchNum', 1, @isnumeric); % aggragate several tasks together
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2021a; matlab -nodisplay -nosplash -nodesktop -nojvm -r', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);
ip.addParameter('language', 'matlab', @ischar); % support matlab, bash

ip.parse(inputFullpaths, outputFullpaths, funcStrs, varargin{:});

% move to the root path
paths = split(which('slurm_cluster_generic_computing_wrapper'), 'LLSM3DTools');
cd(paths{1});

pr = ip.Results;
 % Resolution = pr.Resolution;
jobLogDir = pr.jobLogDir;
tmpDir = pr.tmpDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes = pr.cpuOnlyNodes;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
maxJobNum = pr.maxJobNum;
taskBatchNum = pr.taskBatchNum;
uuid = pr.uuid;
SlurmParam = pr.SlurmParam;
MatlabLaunchStr = pr.MatlabLaunchStr;
language = pr.language;

if isempty(uuid)
    uuid = get_uuid();
end

[dataPath, ~] = fileparts(inputFullpaths{1});

% check if a slurm-based computing cluster exist
if parseCluster 
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
end

nF = numel(inputFullpaths);
is_done_flag = false(nF, 1);
input_exist_mat = batch_file_exist(inputFullpaths);
outputFullpaths = strip(outputFullpaths, 'right', filesep);
output_exist_mat = batch_file_exist(outputFullpaths);
if all(output_exist_mat)
    is_done_flag = ~is_done_flag;
    return;
end

trial_counter = zeros(nF, 1);
if parseCluster
    job_ids = -ones(nF, 1);
    job_status_mat = -ones(nF, 2); % current and previous status
end

loop_counter = 0;
nF_done = 0;
ts = tic;
while ~all(is_done_flag | trial_counter >= maxTrialNum, 'all')
    if parseCluster
        lastP = find(~is_done_flag & trial_counter < maxTrialNum, 1, 'last');
        nB = ceil(nF / taskBatchNum);
    else
        % For no cluster computing, choose the first unfinished one, to avoid 
        % waiting time in each iteration. 
        lastP = find(~is_done_flag & trial_counter < maxTrialNum, 1, 'first');
        nB = nF;
    end
    output_exist_mat(~is_done_flag) = batch_file_exist(outputFullpaths(~is_done_flag));
    
    if parseCluster
        task_ids = rem(1 : nB, 5000);
        job_status_mat(~is_done_flag, 2) = job_status_mat(~is_done_flag, 1);
        job_status_mat(~is_done_flag, 1) = check_batch_slurm_jobs_status(job_ids(~is_done_flag), task_ids(~is_done_flag));
    end
        
    fsnames = cell(1, nF);
    for b = 1 : nB
        fs = (b - 1) * taskBatchNum + 1 : min(b * taskBatchNum, nF);
        task_id = rem(b, 5000);
        
        % check output exist and job status every 1000 batches
        if loop_counter > 0 && rem(b, 1000) == 0
            output_exist_mat(~is_done_flag) = batch_file_exist(outputFullpaths(~is_done_flag));

            if parseCluster
                job_status_mat(~is_done_flag, 2) = job_status_mat(~is_done_flag, 1);
                job_status_mat(~is_done_flag, 1) = check_batch_slurm_jobs_status(job_ids(~is_done_flag), task_ids(~is_done_flag));
            end
        end
        
        for f = fs
            inputFullpath_i = inputFullpaths{f};

            if iscell(inputFullpath_i)
                for fi = 1 : numel(inputFullpath_i)
                    inputFullpath = inputFullpath_i{fi};
                    if ~(exist(inputFullpath, 'file') || exist(inputFullpath, 'dir'))
                        sprintf('%s does not exist, skip it!', inputFullpath);
                        is_done_flag(f) = true;
                        % continue;
                    end
                end
            else
                if ~input_exist_mat(f)
                    inputFullpath = inputFullpath_i;                    
                    sprintf('%s does not exist, skip it!', inputFullpath);
                    is_done_flag(f) = true;
                    % continue;
                end
            end
            
            outputFullpath = outputFullpaths{f};            
            [outputDir, fsname, ext] = fileparts(outputFullpath);
            if isempty(fsname) && ~isempty(ext)
                fsname = ext;
            end
            fsnames{f} = fsname;
            if ~parseCluster
                if isempty(tmpDir)
                    tmpFullpath = sprintf('%s/%s.tmp', outputDir, fsname);
                else
                    tmpFullpath = sprintf('%s/%s.tmp', tmpDir, fsname);            
                end
            end
            
            % if exist(outputFullpath, 'file') || exist(outputFullpath, 'dir')
            if output_exist_mat(f)
                is_done_flag(f) = true;
                if ~parseCluster && exist(tmpFullpath, 'file')
                    delete(tmpFullpath);
                end
                
                % kill new pending jobs
                if  parseCluster && job_ids(f) > 0 && job_status_mat(f, 1) == 0
                    % job_status = check_slurm_job_status(job_ids(f), task_id); 
                    system(sprintf('scancel %d_%d', job_ids(f), task_id), '-echo');
                end
            end
        end
    
        if all(is_done_flag(fs) | trial_counter(fs) >= maxTrialNum)
            continue;
        end
        
        f = fs(end);
        % if parseCluster
            % job_status = check_slurm_job_status(job_ids(f), task_id);
            % job_status_mat(fs, 2) = job_status_mat(fs, 1);
            % job_status_mat(fs, 1) = job_status;
            
            % use master job to run the first pending job to avoid loop
            % through all jobs.
            % if job_status_mat(fs, 1) == 0
                % lastP = b;
            % end            
        % end
                
        % set parameter to skip job submission step in case of reaching max job
        % number and masterCompute is true
        skip_job_submission = false;
        if parseCluster && sum(job_status_mat(~is_done_flag(1 : taskBatchNum : nF), 1) >= 0) >= maxJobNum
            skip_job_submission = true;
        end
        
        func_str = strjoin(funcStrs(fs), ';');
        if ~skip_job_submission && (parseCluster || exist(tmpFullpath, 'file'))
            if parseCluster
                % kill the first pending job and use master node do the computing.
                if job_status_mat(fs, 1) == 0 && (masterCompute && b == lastP)
                    system(sprintf('scancel %d_%d', job_ids(f), task_id), '-echo');
                    trial_counter(fs) = trial_counter(fs) - 1;
                    job_status_mat(fs, 1) = -1;
                    job_status_mat(fs, 2) = -1;
                end

                % if the job is still running, skip it. 
                if job_status_mat(fs, 1) == 1 
                    continue;
                end
                
                % wait some time for the status change
                if loop_counter > 0 && job_status_mat(f, 1) < job_status_mat(f, 2)
                    pause(0.2);
                end

                % If there is no job, submit a job
                if job_status_mat(fs, 1) == -1 && job_status_mat(f, 2) == -1 && ~(masterCompute && f == lastP)
                    if rem(b, 50) == 0 || b == 1
                        fprintf('Task % 4d:    Process %s with function %s... \n', b, strjoin(fsnames(fs), ', '), func_str); 
                    else
                        fprintf('Task % 4d:    Process %s ... \n', b, strjoin(fsnames(fs), ', '));    
                    end
                    if strcmpi(language, 'matlab')
                        matlab_setup_str = 'setup([],true)';

                        matlab_cmd = sprintf('%s;t0_=tic;%s;toc(t0_)', matlab_setup_str, func_str);
                        process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                            '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask, SlurmParam, ...
                            slurm_constraint_str, matlab_cmd, process_cmd);
                    elseif strcmpi(language, 'bash')
                        % process_cmd = func_str;
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                            '--wrap="echo $PWD; echo bash command:  \\\"%s\\\"; %s"'], ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask, SlurmParam, ...
                            slurm_constraint_str, func_str, func_str);
                    end
                    [status, cmdout] = system(cmd, '-echo');
                    if (isempty(cmdout) || isempty(regexp(cmdout, 'Submitted batch job (\d+)\n', 'match'))) ...
                            && strcmpi(language, 'matlab')
                        fprintf('Unable to run matlab code, save the func str to disk and load it to run\n');
                        func_str_dir = sprintf('%s/func_strs/', outputDir);
                        if ~exist(func_str_dir, 'dir')
                            mkdir(func_str_dir);
                        end
                        func_str_fn = sprintf('%s/func_str_f%04d_%s_%s.mat', func_str_dir, f, fsnames{f}, uuid);
                        save('-v7.3', func_str_fn, 'func_str');
                        
                        matlab_cmd = sprintf('%s;t0_=tic;load(''%s'',''func_str'');func_str,feval(str2func([''@()'',func_str]));toc(t0_)', ...
                            matlab_setup_str, func_str_fn);
                        process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                            '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask, SlurmParam, ...
                            slurm_constraint_str, matlab_cmd, process_cmd);
                        [status, cmdout] = system(cmd, '-echo');
                    end
                    
                    job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                    job_id = str2double(job_id{1}{1});
                    job_ids(fs) = job_id;
                    job_status_mat(fs, 1) = 0;
                    trial_counter(fs) = trial_counter(f) + 1;                                
                end
            else
                temp_file_info = dir(tmpFullpath);
                if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 < unitWaitTime
                    continue; 
                else
                    fclose(fopen(tmpFullpath, 'w'));
                end
            end
        else
            if ~parseCluster
                fclose(fopen(tmpFullpath, 'w'));
            end
        end

        if ~parseCluster || (parseCluster && masterCompute && f == lastP)
            fprintf('Process %s with function %s... \n', strjoin(fsnames(fs), ', '), func_str); 
            if parseCluster && loop_counter > 0 && job_status_mat(f, 1) ~= job_status_mat(f, 2)
                pause(1);
            end
            if strcmpi(language, 'matlab')
                % tic; feval(str2func(['@()', func_str])); toc;
                try 
                    t0=tic; feval(str2func(['@()', func_str])); toc(t0);
                catch ME
                    disp(ME)
                    t0=tic; eval(func_str); toc(t0);
                end
            elseif strcmpi(language, 'bash')
                t0=tic; [status, cmdout] = system(func_str, '-echo'); toc(t0)
            end
            fprintf('Done!\n');
        end
    end
    
    if ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
        pause(1);
    end
    if nF_done < sum(is_done_flag)
        nF_done = sum(is_done_flag);
        sprintf('Time %d s: %d / %d (%0.3f) are finished!\n', toc(ts), nF_done, nF, nF_done / nF);
    end
    loop_counter = loop_counter + 1;
end

end


