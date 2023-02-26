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
% xruan (05/21/2022): check job status based on how many pending jobs, rather than each loop
% xruan (08/25/2022): add support for time limit (in hour)


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
ip.addParameter('cpuOnlyNodes', ~true, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);
ip.addParameter('maxJobNum', inf, @isnumeric); % submit limited number of jobs (pending/running)
ip.addParameter('taskBatchNum', 1, @isnumeric); % aggragate several tasks together
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2022a; matlab -nodisplay -nosplash -nodesktop -nojvm -r', @ischar);
ip.addParameter('BashLaunchStr', '', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);
ip.addParameter('jobTimeLimit', 24, @isnumeric); % in hour, [] means no limit
ip.addParameter('language', 'matlab', @ischar); % support matlab, bash

ip.parse(inputFullpaths, outputFullpaths, funcStrs, varargin{:});

% move to the root path
paths = split(which('slurm_cluster_generic_computing_wrapper'), 'LLSM5DTools');
cd(paths{1});
% use the setup within LLSM5DTools
if ~exist('setup.m', 'file')
    cd('LLSM5DTools');
end

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
jobTimeLimit = pr.jobTimeLimit;
MatlabLaunchStr = pr.MatlabLaunchStr;
BashLaunchStr = pr.BashLaunchStr;
language = pr.language;

if isempty(uuid)
    uuid = get_uuid();
end
if isempty(BashLaunchStr)
    BashLaunchStr = 'echo ';
end

[dataPath, ~] = fileparts(inputFullpaths{1});
[outPath, ~] = fileparts(outputFullpaths{1});

% check if a slurm-based computing cluster exist
if parseCluster 
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
    time_str = '';
    if ~isempty(jobTimeLimit) && ~(contains(SlurmParam, ' -t ') || contains(SlurmParam, '--time'))
        % only round to minutes (minimum 1 minute);
        jobTimeLimit = max(jobTimeLimit, 1 / 60);
        h = floor(jobTimeLimit);
        m = round((jobTimeLimit - h) * 60);
        time_str = sprintf(' -t %d:%d:00 ', h, m);
    end
end

nF = numel(inputFullpaths);
is_done_flag = false(nF, 1);
outputFullpaths = strip(outputFullpaths, 'right', filesep);
outputFullpaths = strip(outputFullpaths, 'right', '/');
output_exist_mat = batch_file_exist(outputFullpaths, [], true);
if all(output_exist_mat)
    is_done_flag = ~is_done_flag;
    return;
else
    is_done_flag = output_exist_mat;
end
input_exist_mat = output_exist_mat;
input_exist_mat(~output_exist_mat) = batch_file_exist(inputFullpaths(~output_exist_mat), [], true);

trial_counter = zeros(nF, 1);
if parseCluster
    job_ids = -ones(nF, 1);
    job_status_mat = -ones(nF, 2); % current and previous status
    job_timestamp_mat = zeros(nF, 1);

    nB = ceil(nF / taskBatchNum);        
    task_ids = ones(taskBatchNum, 1) * (1 : nB);
    task_ids = task_ids(:)';
    task_ids = task_ids(1 : nF);
    task_ids = rem(task_ids, 5000);
end

loop_counter = 0;
nF_done = 0;
n_status_check = 10000;
start_time = datetime('now');
ts = tic;
while ~all(is_done_flag | trial_counter >= maxTrialNum, 'all')
    if parseCluster
        job_status_mat(~is_done_flag, 2) = job_status_mat(~is_done_flag, 1);
        job_status_mat(~is_done_flag, 1) = check_batch_slurm_jobs_status(job_ids(~is_done_flag), task_ids(~is_done_flag));
        timestamp = seconds(datetime('now') - start_time);
        job_timestamp_mat(job_status_mat(:, 1) >= 0) = timestamp;
    end
        
    if loop_counter > 0
        output_exist_mat(~is_done_flag) = batch_file_exist(outputFullpaths(~is_done_flag), [], true);
        is_done_flag(~is_done_flag) = output_exist_mat(~is_done_flag);        
    end
    
    if parseCluster
        lastP = find(~is_done_flag & trial_counter < maxTrialNum & job_status_mat(:, 1) < 1, 1, 'last');
        if isempty(lastP)
            lastP = -1;
        else
            lastP = ceil(lastP / taskBatchNum);
        end
    else
        % For no cluster computing, choose the first unfinished one, to avoid 
        % waiting time in each iteration. 
        lastP = find(~is_done_flag & trial_counter < maxTrialNum, 1, 'first');
        nB = nF;
    end    
    fsnames = cell(1, nF);
    for b = 1 : nB
        fs = (b - 1) * taskBatchNum + 1 : min(b * taskBatchNum, nF);
        task_id = rem(b, 5000);
        
        % check output exist and job status every 10000 batches (except the
        % last small bacth (< 0.5 * n_status_check))
        if rem(b, n_status_check) == 0 && (b + n_status_check * 0.5 < nB)      
            if loop_counter == 0 && parseCluster
                numJobs = sum(job_status_mat(~is_done_flag, 1) >= 0);
            end
            if parseCluster && (loop_counter > 0 || (loop_counter == 0 && numJobs >= maxJobNum))
                job_status_mat(:, 2) = job_status_mat(:, 1);
                job_inds = ~is_done_flag & job_status_mat(:, 1) > -1;

                % only check the status of running/pending jobs. 
                output_exist_mat(job_inds) = batch_file_exist(outputFullpaths(job_inds), [], true);
                is_done_flag(job_inds) = output_exist_mat(job_inds);                
                
                job_status_mat(job_inds, 1) = check_batch_slurm_jobs_status(job_ids(job_inds), task_ids(job_inds));
                numJobs = sum(job_status_mat(~is_done_flag, 1) >= 0);
                timestamp = seconds(datetime('now') - start_time);                
                job_timestamp_mat(job_inds) = timestamp;
                if masterCompute
                    if b + n_status_check * 1.5 > nB
                        check_inds = b : nB;
                    else
                        check_inds = b : min(nB, b + n_status_check - 1);
                    end
                    lastPb = find(~is_done_flag(check_inds) & trial_counter(check_inds) < maxTrialNum & job_status_mat(check_inds, 1) < 0, 1, 'last');
                    if ~isempty(lastPb)
                        lastP = lastPb + b - 1;
                    end
                end                
            elseif ~parseCluster
                output_exist_mat(~is_done_flag) = batch_file_exist(outputFullpaths(~is_done_flag), [], true);
                is_done_flag(~is_done_flag) = output_exist_mat(~is_done_flag);                
            end
        end

        if all(is_done_flag(fs) | trial_counter(fs) >= maxTrialNum)
            continue;
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
                if  parseCluster && (numel(fs) == 1 || (numel(fs) > 1 && all(output_exist_mat(fs)))) ...
                        && job_ids(f) > 0 && job_status_mat(f, 1) == 0
                    % job_status = check_slurm_job_status(job_ids(f), task_id); 
                    system(sprintf('scancel %d_%d', job_ids(f), task_id), '-echo');
                end
            end
        end
    
        if all(is_done_flag(fs) | trial_counter(fs) >= maxTrialNum)
            continue;
        end
        
        f = fs(end);
                
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
                pending_flag = false;
                if job_status_mat(f, 1) == 0 && (masterCompute && b == lastP)
                    system(sprintf('scancel %d_%d', job_ids(f), task_id), '-echo');
                    pending_flag = true;
                    trial_counter(fs) = trial_counter(fs) - 1;
                    job_status_mat(fs, 1) = -1;
                    job_status_mat(fs, 2) = -1;
                end

                % if the job is still running, skip it. 
                if job_status_mat(f, 1) == 1 
                    continue;
                end
                
                % check if the job just finished, if masterCompute 30s, if not, 45s
                if loop_counter > 0 && trial_counter(f) > 0 && job_status_mat(f, 1) == -1 ...
                        && timestamp - job_timestamp_mat(f) < 45 - masterCompute * 15
                    continue;
                end

                % If there is no job, submit a job
                if job_status_mat(f, 1) == -1 && job_status_mat(f, 2) == -1 && ~(masterCompute && f == lastP)
                    if rem(b, 50) == 0 || b == 1
                        fprintf('Task % 4d:    Process %s with function %s... \n', b, strjoin(fsnames(fs), ', '), func_str); 
                    else
                        fprintf('Task % 4d:    Process %s ... \n', b, strjoin(fsnames(fs), ', '));    
                    end
                    if strcmpi(language, 'matlab')
                        matlab_setup_str = 'setup([],true)';

                        matlab_cmd = sprintf('%s;t0_=tic;%s;toc(t0_)', matlab_setup_str, func_str);
                        process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s %s ', ...
                            '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask, SlurmParam, ...
                            slurm_constraint_str, time_str, matlab_cmd, process_cmd);
                    elseif strcmpi(language, 'bash')
                        % process_cmd = func_str;
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s %s ', ...
                            '--wrap="echo $PWD; echo bash command:  \\\"%s\\\";  %s; %s"'], ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask, SlurmParam, ...
                            slurm_constraint_str, time_str, func_str, BashLaunchStr, func_str);
                    end
                    [status, cmdout] = system(cmd, '-echo');
                    if isempty(cmdout) || isempty(regexp(cmdout, 'Submitted batch job (\d+)\n', 'match'))
                        fprintf('Unable to run the code, save the func str to disk and load it to run\n');
                        func_str_dir = sprintf('%s/func_strs/', tmpDir);
                        if isempty(tmpDir)
                            func_str_dir = sprintf('%s/func_strs/', outPath);
                        end
                        if ~exist(func_str_dir, 'dir')
                            mkdir(func_str_dir);
                        end
                        
                        if strcmpi(language, 'matlab')
                            func_str_fn = sprintf('%s/func_str_f%04d_%s_%s.mat', func_str_dir, f, fsnames{f}, uuid);                            
                            save('-v7.3', func_str_fn, 'func_str');
                            
                            matlab_cmd = sprintf('%s;t0_=tic;load(''%s'',''func_str'');func_str,feval(str2func([''@()'',func_str]));toc(t0_)', ...
                                matlab_setup_str, func_str_fn);
                            process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                            cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                                '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], ...
                                task_id, job_log_fname, job_log_error_fname, cpusPerTask, SlurmParam, ...
                                slurm_constraint_str, matlab_cmd, process_cmd);                            
                        elseif strcmpi(language, 'bash')
                            func_str_fn = sprintf('%s/func_str_f%04d_%s_%s.sh', func_str_dir, f, fsnames{f}, uuid);     
                            fid = fopen(func_str_fn, 'w');
                            % escape '\ ' for writing to the disk
                            func_str = strrep(func_str, '\ ', '\\ ');
                            fprintf(fid, func_str);
                            fclose(fid);
                            
                            cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                                '--wrap="echo $PWD; echo bash command:  \\\"%s\\\"; %s; bash %s"'], ...
                                task_id, job_log_fname, job_log_error_fname, cpusPerTask, SlurmParam, ...
                                slurm_constraint_str, func_str_fn, BashLaunchStr, func_str_fn);
                        end
                        [status, cmdout] = system(cmd, '-echo');
                    end
                    
                    job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                    if isempty(job_id)
                        continue;
                    end
                    job_id = str2double(job_id{1}{1});
                    job_ids(fs) = job_id;
                    job_status_mat(fs, 1) = 0;
                    trial_counter(fs) = trial_counter(fs) + 1;                                
                end
            else
                temp_file_info = dir(tmpFullpath);
                if minutes(datetime('now') -  datetime(temp_file_info.date)) < unitWaitTime
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

        if ~parseCluster || (parseCluster && masterCompute && b == lastP)
            if parseCluster && loop_counter > 0
                % for nonpending killed jobs, wait a bit longer in case of just finished job.
                if ~pending_flag
                    pause(1);
                    if timestamp - job_timestamp_mat(f) < 30
                        continue;
                    end
                    if exist(outputFullpath, 'file') || exist(outputFullpath, 'dir')
                        continue;
                    end
                end
            end
            fprintf('Process %s with function %s... \n', strjoin(fsnames(fs), ', '), func_str);             
            if strcmpi(language, 'matlab')
                try 
                    t0=tic; feval(str2func(['@()', func_str])); t1=toc(t0);
                catch ME
                    disp(ME)
                    t0=tic; eval(func_str); t1=toc(t0);
                end
            elseif strcmpi(language, 'bash')
                t0=tic; [status, cmdout] = system(func_str, '-echo'); t1=toc(t0);
            end
            fprintf('Elapsed time is %f seconds.\n', t1);
            trial_counter(fs) = trial_counter(fs) + 1;
            if ~parseCluster && exist(outputFullpath, 'file') && exist(tmpFullpath, 'file')
                delete(tmpFullpath);
            end

            fprintf('Done!\n');
        end
    end
    
    if parseCluster && ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
        if masterCompute
            pause(0.1);
        else
            pause(1);
        end
    end
    if nF_done < sum(is_done_flag)
        nF_done = sum(is_done_flag);
        fprintf('Time %0.2f s: %d / %d (%0.2f%%) are finished!\n', toc(ts), nF_done, nF, nF_done / nF * 100);
    end
    loop_counter = loop_counter + 1;
end

end


