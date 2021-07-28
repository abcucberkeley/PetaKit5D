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


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('outputFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('functionStrs', @(x) iscell(x) || ischar(x));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @isstr);
ip.addParameter('tmpDir', '', @isstr);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @isstr);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2020b; matlab -nodisplay -nosplash -nodesktop -nojvm -r', @ischar);
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
trial_counter = zeros(nF, 1);
if parseCluster
    job_ids = -ones(nF, 1);
    job_status_mat = -ones(nF, 2); % current and previous status
end

loop_counter = 0;
while ~all(is_done_flag | trial_counter >= maxTrialNum, 'all')
    if parseCluster
        lastP = find(~is_done_flag & trial_counter < maxTrialNum, 1, 'last');
    else
        % For no cluster computing, choose the first unfinished one, to avoid 
        % waiting time in each iteration. 
        lastP = find(~is_done_flag & trial_counter < maxTrialNum, 1, 'first');
    end
    
    for f = 1 : nF
        inputFullpath_i = inputFullpaths{f};
        task_id = rem(f, 5000);
        
        if iscell(inputFullpath_i)
            for fi = 1 : numel(inputFullpath_i)
                inputFullpath = inputFullpath_i{fi};
                if ~(exist(inputFullpath, 'file') || exist(inputFullpath, 'dir'))
                    sprintf('%s does not exist, skip it!', inputFullpath);
                    is_done_flag(f) = true;
                    continue;
                end
            end
        else
            inputFullpath = inputFullpath_i;
            if ~(exist(inputFullpath, 'file') || exist(inputFullpath, 'dir'))
                sprintf('%s does not exist, skip it!', inputFullpath);
                is_done_flag(f) = true;
                continue;
            end
        end
        if is_done_flag(f) || trial_counter(f) >= maxTrialNum 
            continue;
        end
        outputFullpath = outputFullpaths{f};
        if strcmp(outputFullpath(end), filesep)
            outputFullpath = outputFullpath(1 : end - 1);
        end
        [outputDir, fsname, ext] = fileparts(outputFullpath);
        if isempty(fsname) && ~isempty(ext)
            fsname = ext;
        end
        if isempty(tmpDir)
            tmpFullpath = sprintf('%s/%s.tmp', outputDir, fsname);
        else
            tmpFullpath = sprintf('%s/%s.tmp', tmpDir, fsname);            
        end
        if exist(outputFullpath, 'file') || exist(outputFullpath, 'dir')
            is_done_flag(f) = true;
            if exist(tmpFullpath, 'file')
                delete(tmpFullpath);
            end
            
            % kill new pending jobs
            if parseCluster && job_ids(f) > 0
                job_status = check_slurm_job_status(job_ids(f), task_id); 
                if job_status ~= 1
                    system(sprintf('scancel %d_%d', job_ids(f), task_id), '-echo');
                end
            end

            continue;
        end
        
        func_str = funcStrs{f};
        if exist(tmpFullpath, 'file') || parseCluster
            if parseCluster
                job_status = check_slurm_job_status(job_ids(f), task_id);
                job_status_mat(f, 2) = job_status_mat(f, 1);
                job_status_mat(f, 1) = job_status;

                % kill the first pending job and use master node do the computing.
                if job_status == 0.5 && (masterCompute && f == lastP)
                    system(sprintf('scancel %d_%d', job_ids(f), task_id), '-echo');
                    trial_counter(f) = trial_counter(f) - 1;
                end

                % if the job is still running, skip it. 
                if job_status == 1 
                    continue;
                end
                
                % wait some time for the status change
                if loop_counter > 0 && job_status_mat(f, 1) < job_status_mat(f, 2)
                    pause(5);
                end

                % If there is no job, submit a job
                if job_status == -1 && job_status_mat(f, 2) == -1 && ~(masterCompute && f == lastP)
                    if rem(f, 50) == 0 || f == 1
                        fprintf('Task % 4d:    Process %s with function %s... \n', f, fsname, func_str); 
                    else
                        fprintf('Task % 4d:    Process %s ... \n', f, fsname);    
                    end
                    if strcmpi(language, 'matlab')
                        matlab_setup_str = 'setup([],true)';

                        matlab_cmd = sprintf('%s;tic;%s;toc', matlab_setup_str, func_str);
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
                        func_str_fn = sprintf('%s/func_str_f%04d_%s_%s.mat', func_str_dir, f, fsname, uuid);
                        save('-v7.3', func_str_fn, 'func_str');
                        
                        matlab_cmd = sprintf('%s;tic;load(''%s'',''func_str'');func_str,feval(str2func([''@()'',func_str]));toc', matlab_setup_str, func_str_fn);
                        process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                            '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask, SlurmParam, ...
                            slurm_constraint_str, matlab_cmd, process_cmd);
                        [status, cmdout] = system(cmd, '-echo');
                    end
                    
                    job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                    job_id = str2double(job_id{1}{1});
                    job_ids(f) = job_id;
                    trial_counter(f) = trial_counter(f) + 1;                                
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
            fclose(fopen(tmpFullpath, 'w'));
        end

        if ~parseCluster || (parseCluster && masterCompute && f == lastP)
            fprintf('Process %s with function %s... \n', fsname, func_str); 
            if loop_counter > 0
                pause(5);
            end
            if strcmpi(language, 'matlab')
                % tic; feval(str2func(['@()', func_str])); toc;
                try 
                    tic; feval(str2func(['@()', func_str])); toc;
                catch ME
                    disp(ME)
                    tic; eval(func_str); toc;
                end
            elseif strcmpi(language, 'bash')
                tic; [status, cmdout] = system(func_str, '-echo'); toc
            end
            fprintf('Done!\n');
        end
        % toc
        if exist(outputFullpath, 'file') || exist(outputFullpath, 'dir')
            is_done_flag(f) = true;
            if exist(tmpFullpath, 'file')
                delete(tmpFullpath);
            end
        end
    end
    
    if ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
        pause(5);
    end
    loop_counter = loop_counter + 1;
end


end
