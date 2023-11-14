function [is_done_flag] = mcc_slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, funcStrs, varargin)
% generic slurm computing framework that use a function handle/string as input
% for the computing for mcc functions
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
% xruan (11/11/2023): add support for minimum query interval for jobs and file system


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('outputFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('functionStrs', @(x) iscell(x) || ischar(x));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('tmpDir', '', @ischar);
ip.addParameter('maxCPUNum', 24, @isnumeric);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('GPUJob', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);
ip.addParameter('maxJobNum', inf, @isnumeric); % submit limited number of jobs (pending/running)
ip.addParameter('taskBatchNum', 1, @isnumeric); % aggragate several tasks together
ip.addParameter('BashLaunchStr', '', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);
ip.addParameter('SlurmConstraint', '', @ischar);
ip.addParameter('MCRCacheRoot', '/tmp/', @ischar);
ip.addParameter('MCRParam', '/usr/local/MATLAB/R2022b', @ischar);
ip.addParameter('MCCMasterStr', '/home/xruan/Projects/XR_Repository/mcc/run_mccMaster.sh', @ischar);
ip.addParameter('jobTimeLimit', 24, @isnumeric); % in hour, [] means no limit
ip.addParameter('queryInterval', 3, @isnumeric); % in second, by default 3 s. 
ip.addParameter('language', 'bash', @ischar); % support matlab, bash
ip.addParameter('GNUparallel', false, @islogical); % support matlab, bash
ip.addParameter('paraJobNum', 1, @isnumeric); % support matlab, bash
ip.addParameter('masterParaFactor', 1, @isnumeric); % master job parallel job percentage

ip.parse(inputFullpaths, outputFullpaths, funcStrs, varargin{:});

% move to the root path
funcFn = which(mfilename);
paths = split(funcFn, 'LLSM5DTools');
cd(paths{1});
setupFn = [paths{1}, 'setup.m'];
% use the setup within LLSM5DTools
if ismcc || isdeployed
    if ~ispc
        cd('/tmp/');
    else
        if ~exist(setupFn, 'file')
            cd('../');
        end
    end
else
    if ~exist(setupFn, 'file')
        paths = split(funcFn, 'utils/generic_computing_framework');
        cd(paths{1});    
    end
end

pr = ip.Results;
 % Resolution = pr.Resolution;
jobLogDir = pr.jobLogDir;
tmpDir = pr.tmpDir;
maxCPUNum = pr.maxCPUNum;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
GPUJob = pr.GPUJob;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
maxJobNum = pr.maxJobNum;
taskBatchNum = pr.taskBatchNum;
uuid = pr.uuid;
SlurmParam = pr.SlurmParam;
SlurmConstraint = pr.SlurmConstraint;
MCRCacheRoot = pr.MCRCacheRoot;
MCRParam = pr.MCRParam;
MCCMasterStr = pr.MCCMasterStr;
jobTimeLimit = pr.jobTimeLimit;
queryInterval = pr.queryInterval;
BashLaunchStr = pr.BashLaunchStr;
GNUparallel = pr.GNUparallel;
paraJobNum = pr.paraJobNum;
masterParaFactor = pr.masterParaFactor;

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
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath, jobLogDir);
    time_str = '';
    if ~isempty(jobTimeLimit) && ~(contains(SlurmParam, ' -t ') || contains(SlurmParam, '--time'))
        % only round to minutes (minimum 1 minute);
        jobTimeLimit = max(jobTimeLimit, 1 / 60);
        h = floor(jobTimeLimit);
        m = round((jobTimeLimit - h) * 60);
        time_str = sprintf(' -t %d:%02d:00 ', h, m);
    end
end

nF = numel(inputFullpaths);
is_done_flag = false(nF, 1);
outputFullpaths = strip(outputFullpaths, 'right', filesep);
outputFullpaths = strip(outputFullpaths, 'right', '/');
output_exist_mat = batch_file_exist(outputFullpaths, [], true);
if all(output_exist_mat)
    is_done_flag = ~is_done_flag;
    fprintf('All output files (%d / %d) already exists!\n', nF, nF);        
    return;
else
    is_done_flag = output_exist_mat;
end
input_exist_mat = output_exist_mat;
input_exist_mat(~output_exist_mat) = batch_file_exist(inputFullpaths(~output_exist_mat), [], true);

% write mcc strings to disk
dt = char(datetime('now', 'Format', 'yyyyMMdd_HHmmSS'));

if isempty(tmpDir)
    funcInputDir = sprintf('%s/tmp/%s/', dataPath, dt);
else
    funcInputDir = sprintf('%s/%s/', tmpDir, dt);            
end
mkdir_recursive(funcInputDir);

batchSize = taskBatchNum;
nB = ceil(nF / batchSize);

for b = 1 : nB
    inputFn = sprintf('%s/input_%d.txt', funcInputDir, b);
    
    s = (b - 1) * batchSize + 1;
    t = min(b * batchSize, nF);
    if all(is_done_flag(s : t))
        continue;
    end

    fid = fopen(inputFn, 'w');

    for f = s : t
        func_str = funcStrs{f};
        [func_name, var_str] = convert_function_string_to_mcc_string(func_str);
        % tline = sprintf('%s %s %s %s \n', MCCMasterStr, MCRParam, func_name, var_str);
        % check output file in bash to avoid waste of time in loading mcc program if the output file exists. 
        tline = sprintf('if [ ! -f %s ]; then MCR_CACHE_ROOT=%s %s %s %s %s ; else echo output %s already exists; fi\n', ...
            outputFullpaths{f}, MCRCacheRoot, MCCMasterStr, MCRParam, func_name, var_str, outputFullpaths{f});
        fprintf(fid, tline);
    end
    fclose(fid);
end

% setup jobs
trial_counter = zeros(nF, 1);
if parseCluster
    job_ids = -ones(nF, 1);
    job_status_mat = -ones(nF, 2); % current and previous status
    job_timestamp_mat = zeros(nF, 1);

    task_ids = ones(taskBatchNum, 1) * (1 : nB);
    task_ids = task_ids(:);
    task_ids = task_ids(1 : nF);
    task_ids = rem(task_ids, 5000);
end

fprintf('Task number : %d, task batch size : %d, task parallel number : %d, job number : %d \n', ...
    nF, taskBatchNum, paraJobNum, nB);

GNUparallel = GNUparallel & (paraJobNum > 1);

loop_counter = 0;
nF_done = 0;
n_status_check = 10000;
start_time = datetime('now');
ts = tic;
qts = ts;
while (~parseCluster && ~all(is_done_flag | trial_counter >= maxTrialNum, 'all')) || ...
        (parseCluster && ~all(is_done_flag | (trial_counter >= maxTrialNum & job_status_mat(:, 1) < 0), 'all'))
    querySystem = ~parseCluster || (parseCluster && (loop_counter == 0 || toc(qts) > queryInterval));    
    if parseCluster && querySystem
        job_status_mat(~is_done_flag, 2) = job_status_mat(~is_done_flag, 1);
        job_status_mat(~is_done_flag, 1) = check_batch_slurm_jobs_status(job_ids(~is_done_flag), task_ids(~is_done_flag));
        timestamp = seconds(datetime('now') - start_time);
        job_timestamp_mat(job_status_mat(:, 1) >= 0) = timestamp;
        qts = tic;        
    end
        
    if loop_counter > 0 && querySystem
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
        if parseCluster
            job_status_mat(fs, :) = repmat(min(job_status_mat(fs, :)), numel(fs), 1);
        end
        
        % check output exist and job status every 10000 batches (except the
        % last small bacth (< 0.5 * n_status_check))
        if rem(b, n_status_check) == 0 && (b + n_status_check * 0.5 < nB)
            querySystem = ~parseCluster || (parseCluster &&  toc(qts) > queryInterval);            
            if parseCluster && querySystem
                job_status_mat(:, 2) = job_status_mat(:, 1);
                job_inds = ~is_done_flag & job_status_mat(:, 1) > -1;

                % only check the status of running/pending jobs. 
                output_exist_mat(job_inds) = batch_file_exist(outputFullpaths(job_inds), [], true);
                is_done_flag(job_inds) = output_exist_mat(job_inds);                
                
                job_status_mat(job_inds, 1) = check_batch_slurm_jobs_status(job_ids(job_inds), task_ids(job_inds));
                timestamp = seconds(datetime('now') - start_time);                
                job_timestamp_mat(job_inds) = timestamp;
                qts = tic;                
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
        inputFn = sprintf('%s/input_%d.txt', funcInputDir, b);

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
                if GPUJob
                    cpusPerTask_f = min(maxCPUNum, cpusPerTask * (trial_counter(f) + 1));
                else
                    cpusPerTask_f = min(maxCPUNum, cpusPerTask);
                end

                % If there is no job, submit a job
                if job_status_mat(f, 1) == -1 && job_status_mat(f, 2) == -1 && ~(masterCompute && b == lastP)
                    if rem(b, 50) == 0 || b == 1
                        fprintf('Task % 4d:    Process %s with function %s... \n', b, strjoin(fsnames(fs), ', '), func_str); 
                    else
                        fprintf('Task % 4d:    Process %s ... \n', b, strjoin(fsnames(fs), ', '));    
                    end
                    % if tried twice, still fail, not use parallel computing
                    if ~GNUparallel || trial_counter(f) > 1
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s %s ', ...
                            '--wrap="echo $PWD; echo bash command:  \\\"%s\\\";  %s; bash %s"'], ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask_f, SlurmParam, ...
                            SlurmConstraint, time_str, inputFn, BashLaunchStr, inputFn);
                    else
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d ', ...
                            '--ntasks=1 %s %s %s --wrap="echo $PWD; echo bash command: \\\"%s\\\"; ', ...
                            '%s; parallel --delay 0.2 \\\"srun --exclusive -c %d bash -c {}\\\" < %s"'], ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask_f, SlurmParam, ...
                            SlurmConstraint, time_str, inputFn, BashLaunchStr, floor(cpusPerTask_f / paraJobNum), inputFn);
                        % in case of some jobs fail because of memory issue, directly use parallel for computing
                        if trial_counter(f) <= 1
                            paraJobNum_f = max(1, round(paraJobNum / (trial_counter(f) + 1)));
                            cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d ', ...
                                '--ntasks=1 %s %s %s --wrap="echo $PWD; echo bash command: \\\"%s\\\"; ', ...
                                '%s; parallel --ungroup --jobs %d --delay 0.2 < %s"'], ...
                                task_id, job_log_fname, job_log_error_fname, cpusPerTask_f, SlurmParam, ...
                                SlurmConstraint, time_str, inputFn, BashLaunchStr, paraJobNum_f, inputFn);
                            if GPUJob
                                cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d ', ...
                                    '--ntasks=1 %s %s %s --wrap="echo $PWD; echo bash command: \\\"%s\\\"; ', ...
                                    '%s; parallel --ungroup --jobs %d --delay 0.2 CUDA_VISIBLE_DEVICES=''\\$(({%%} - 1))'' eval {} < %s"'], ...
                                    task_id, job_log_fname, job_log_error_fname, cpusPerTask_f, SlurmParam, ...
                                    SlurmConstraint, time_str, inputFn, BashLaunchStr, paraJobNum, inputFn);
                            end
                        end
                    end
                    [status, cmdout] = system(cmd, '-echo');
                    if isempty(cmdout) || isempty(regexp(cmdout, 'Submitted batch job (\d+)\n', 'match'))
                        fprintf('Unable to run the code, save the func str to disk and load it to run\n');
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
                % change the wait time to the maximum of 30s and half of computing time
                if ~pending_flag
                    pause(1);
                    if timestamp - job_timestamp_mat(f) < min(max(30, t1 * 0.5), 180)
                        continue;
                    end
                    if all(is_done_flag(fs))
                        continue;
                    end
                end
            end
            fprintf('Process %s with function %s... \n', strjoin(fsnames(fs), ', '), func_str);             
            if ~GNUparallel
                % process_cmd = func_str;
                cmd = sprintf(['%s; bash %s'], BashLaunchStr, inputFn);
            else
                if GPUJob
                    paraJobNum_f = max(1, min(paraJobNum - 1, round(paraJobNum * masterParaFactor)));
                else
                    paraJobNum_f = max(1, min(paraJobNum - 1, round(paraJobNum * masterParaFactor / (trial_counter(f) + 1))));
                end
                cmd = sprintf(['%s; parallel --jobs %d --delay 0.2 < %s'], BashLaunchStr, paraJobNum_f, inputFn);
                if ismcc || isdeployed
                    % reduce the load of master job in case of crash due to oom
                    cmd = sprintf(['%s; parallel --ungroup --jobs %d --delay 0.2 < %s'], BashLaunchStr, paraJobNum_f, inputFn);
                end
            end
            t0=tic; [status, cmdout] = system(cmd, '-echo'); t1=toc(t0);
            fprintf('Elapsed time is %f seconds.\n', t1);
            trial_counter(fs) = trial_counter(fs) + 1;
            if ~parseCluster && exist(outputFullpath, 'file') && exist(tmpFullpath, 'file')
                delete(tmpFullpath);
            end

            fprintf('Done!\n');
        end
    end
    
    if parseCluster && ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
        wt = max(0.001, queryInterval - toc(qts) - 0.001);
        % if all jobs are still pending, wait another query interval        
        if ~masterCompute && ~any(job_status_mat(~is_done_flag, 1), 'all')
            wt = wt + queryInterval;
        end
        pause(wt);
    end
    if nF_done < sum(is_done_flag)
        nF_done = sum(is_done_flag);
        fprintf('Time %0.2f s: %d / %d (%0.2f%%) are finished!\n', toc(ts), nF_done, nF, nF_done / nF * 100);
    end
    loop_counter = loop_counter + 1;
end

if all(is_done_flag)
    fprintf('All output files (%d / %d) are finished!\n', nF, nF);    
end

end

