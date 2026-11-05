function [job_id, job_status, submit_status] = generic_single_job_submit_wrapper(funcStr, job_id, task_id, varargin)
% generic computing framework that use a function handle/string as input
% for the computing for integration of all computing schemes, including
% slurm, slurm mcc, and their corresponding single job submission. We will
% integrate more frameworks in the future
% 
% Author Xiongtao Ruan (03/09/2023)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('funcStr', @(x) ischar(x));
ip.addRequired('job_id', @(x) isnumeric(x) || ischar(x));
ip.addRequired('task_id', @(x) isnumeric(x) || ischar(x));
ip.addParameter('clusterType', 'slurm', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', false, @islogical); % master node participate in the task computing. 
ip.addParameter('lastFile', false, @islogical); % last file for master job
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('jobLogFname', '../job_logs/job_%A_%a.out', @ischar);
ip.addParameter('jobErrorFname', '../job_logs/job_%A_%a.err', @ischar);
ip.addParameter('tmpDir', '', @ischar);
ip.addParameter('maxCPUNum', 24, @isnumeric);
ip.addParameter('minCPUNum', 1, @isnumeric);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('memPerCPU', 20.9, @isnumeric);
ip.addParameter('memAllocate', [], @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);
ip.addParameter('maxJobNum', inf, @isnumeric); % submit limited number of jobs (pending/running)
ip.addParameter('taskBatchNum', 1, @isnumeric); % aggragate several tasks together
% -batch (non-interactive): with -r, MATLAB R2024a and later installs and starts
% the MathWorks Service Host in the job (in ~/.MathWorks, ~600 MB per node), and
% exiting took from 6 s to over 15 minutes on the cluster
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2024b; matlab -nodisplay -nosplash -nodesktop -nojvm -batch', @ischar);
ip.addParameter('BashLaunchStr', '', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);
ip.addParameter('SlurmConstraint', '', @ischar);
ip.addParameter('MCRCacheRoot', '/tmp/', @ischar);
ip.addParameter('MCRParam', '/usr/local/MATLAB/R2023a', @ischar);
ip.addParameter('MCCMasterStr', '/home/xruan/Projects/XR_Repository/mcc/run_mccMaster.sh', @ischar);
ip.addParameter('jobTimeLimit', 24, @isnumeric); % in hour, [] means no limit
ip.addParameter('language', 'matlab', @ischar); % support matlab, bash
ip.addParameter('GNUparallel', false, @islogical); % support matlab, bash
ip.addParameter('paraJobNum', 1, @isnumeric); % support matlab, bash
ip.addParameter('ConfigFile', '', @ischar); % cluster configuration file that override the default parameters
ip.addParameter('taskMetrics', false, @islogical); % record the resources of the task in <job log folder>/task_metrics/ (see task_metrics)
ip.addParameter('taskRunner', true, @islogical); % after its task, a job's process runs unstarted tasks of other single jobs of the same kind (see task_runner_run); off if the caller or the config says so
ip.addParameter('SlurmAlternatives', [], @(x) isempty(x) || isstruct(x) || iscell(x)); % other places for the job (entries with SlurmParam, memPerCPU, maxCPUNum, minCPUNum, gpuHeadroom; see slurm_alternatives)
ip.addParameter('gpuHeadroom', [], @isnumeric); % with SlurmAlternatives: CPUs and GB per idle GPU that CPU-only jobs leave free at the main place ([]: the node's share per GPU)

ip.parse(funcStr, job_id, task_id, varargin{:});

% move to the root path
funcFn = which(mfilename);
paths = split(funcFn, 'PetaKit5D');
if numel(paths) > 1 && exist(paths{1}, 'dir')
    cd(paths{1});
end
setupFn = [paths{1}, 'setup.m'];
% use the setup within PetaKit5D
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
        if ispc
            paths = split(funcFn, 'utils\generic_computing_framework');
        else
            paths = split(funcFn, 'utils/generic_computing_framework');        
        end
        cd(paths{1});    
    end
end

pr = ip.Results;
ConfigFile = pr.ConfigFile;

persistent ConfigFile_orig confData confModDate;
if ~isempty(ConfigFile)
    dir_info = dir(ConfigFile);
    conf_date = dir_info.date;    
    if ~strcmp(ConfigFile, ConfigFile_orig) || isempty(confData) || isempty(confModDate) || any(confModDate ~= conf_date)
        ConfigFile_orig = ConfigFile;
        confModDate = conf_date;        
        fprintf('Set up cluster configuration according to %s...\n', ConfigFile);
        [~, ~, ext] = fileparts(ConfigFile);
        switch ext 
            case '.json'
                fid = fopen(ConfigFile);
                raw = fread(fid);
                fclose(fid);
                str = char(raw');
                confData = jsondecode(str);            
            case '.mat'
                a = load(ConfigFile);
                confData = a;
        end
        disp('Cluster configuration:');
        disp(confData);        
    end

    conf_keys = fields(confData);
    pr_keys = fields(pr);
    for i = 1 : numel(conf_keys)
        % keys match the parameter names regardless of case, as name-value
        % arguments do (e.g., MemPerCPU sets memPerCPU)
        key = pr_keys(strcmpi(pr_keys, conf_keys{i}));
        if isempty(key)
            key = conf_keys(i);
        end
        pr.(key{1}) = confData.(conf_keys{i});
    end
end

clusterType = pr.clusterType;
mccMode = pr.mccMode;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
lastFile = pr.lastFile;
jobLogDir = pr.jobLogDir;
jobLogFname = pr.jobLogFname;
jobErrorFname = pr.jobErrorFname;
% empty job log paths (as "jobLogDir": "" in a config file) are the default
% ones, so that no job gets "-o  -e" (its log went to a file named -e)
if isempty(jobLogDir)
    jobLogDir = '../job_logs';
end
if isempty(jobLogFname)
    jobLogFname = '../job_logs/job_%A_%a.out';
end
if isempty(jobErrorFname)
    jobErrorFname = '../job_logs/job_%A_%a.err';
end
tmpDir = pr.tmpDir;
maxCPUNum = pr.maxCPUNum;
minCPUNum = pr.minCPUNum;
cpusPerTask = pr.cpusPerTask;
memPerCPU = pr.memPerCPU;
memAllocate = pr.memAllocate;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
maxJobNum = pr.maxJobNum;
taskBatchNum = pr.taskBatchNum;
MatlabLaunchStr = pr.MatlabLaunchStr;
BashLaunchStr = pr.BashLaunchStr;
SlurmParam = pr.SlurmParam;
SlurmConstraint = pr.SlurmConstraint;
MCRCacheRoot = pr.MCRCacheRoot;
MCRParam = pr.MCRParam;
MCCMasterStr = pr.MCCMasterStr;
jobTimeLimit = pr.jobTimeLimit;
language = pr.language;
GNUparallel = pr.GNUparallel;
paraJobNum = pr.paraJobNum;
taskMetrics = pr.taskMetrics;
taskRunner = ip.Results.taskRunner && pr.taskRunner;
SlurmAlternatives = pr.SlurmAlternatives;
gpuHeadroom = pr.gpuHeadroom;

% the Slurm check and the job log paths, kept for the job log folder while
% Slurm is there (callers call this for each job in each round; sinfo and
% the folder's listing took 10 ms or more per call)
persistent logDirCheck;
if exist(jobLogDir, 'dir')
    dataPath = tempdir();
    if ispc
        dataPath = strrep(dataPath, '\', '/');
    end
    logDirKey = [pwd, '|', jobLogDir];
    if isstruct(logDirCheck) && strcmp(logDirCheck.key, logDirKey)
        jobLogFname = logDirCheck.jobLogFname;
        jobErrorFname = logDirCheck.jobErrorFname;
    else
        % (if sinfo fails, checkSlurmCluster gives no log paths: the job
        % keeps its own, else it got "-o  -e" and its log went to a file
        % named -e)
        [parseCluster, logFname, errorFname] = checkSlurmCluster(dataPath, jobLogDir);
        if parseCluster
            jobLogFname = logFname;
            jobErrorFname = errorFname;
            logDirCheck = struct('key', logDirKey, 'jobLogFname', jobLogFname, 'jobErrorFname', jobErrorFname);
        end
    end
end

% the kinds of single jobs (see below) this session submitted a task of
persistent kindsSeen;
% the places for the job (with SlurmAlternatives, made when it is submitted)
places = [];
submit_status = 0;
trial_counter = 0;
poolTask = [];
while submit_status == 0 && trial_counter < maxTrialNum
    switch clusterType
        case 'slurm'
            array_id = rem(task_id, 5000);
            % (a task without a job yet: squeue could only say so)
            job_status = -1;
            if job_id > 0
                job_status = check_slurm_job_status(job_id, array_id);
                % a job submitted to one of several places moves to another
                % one where it can start now (see slurm_alternatives)
                if job_status == 0 && ~(masterCompute && lastFile)
                    slurm_alternatives('move', job_id, array_id);
                elseif job_status < 0
                    slurm_alternatives('forget', job_id, array_id);
                end
            end

            % task runners: after the job of a task submitted to a pool (see
            % below) ended, a runner of another job may still run the task
            % (then it counts as running); before the master submits it
            % again (as before) or runs it itself (last file), it takes the
            % task back from the pool, so that no runner takes it. While the
            % job is pending, a task that a runner failed is taken back (the
            % job is cancelled) and submitted again, and one whose runner
            % ended without finishing it is left to the job
            % (none before this session's first pool task: a session
            % without runners doesn't load task_runner_run)
            poolTask = [];
            if job_status < 1 && single_job_pools()
                poolTask = task_runner_run('single_task', funcStr);
            end
            if ~isempty(poolTask) && job_status == 0 && ~(masterCompute && lastFile)
                if strcmp(poolTask.jobKey, sprintf('%d_%d', job_id, array_id)) ...
                        && task_runner_run('single_status', job_id, array_id, 0) == -1
                    job_status = -1;
                end
                poolTask = [];
            end
            if ~isempty(poolTask)
                [live, claimed] = task_runner_run('single_live', funcStr);
                if live
                    job_status = 1;
                    return;
                end
                % (a task added in this call is submitted again after a
                % failed submission)
                if ~strcmp(poolTask.state, 'new') || claimed || (masterCompute && lastFile)
                    if ~task_runner_run('single_take', funcStr)
                        job_status = 1;
                        return;
                    end
                    poolTask = [];
                end
            end

            % kill the last pending job and use master node do the computing.
            if job_status == 0 && (masterCompute && lastFile)
                system(sprintf('scancel %d_%d', job_id, array_id), '-echo');
            end
    
            % if the job is pending or running, skip it; the caller checks
            % it again in its next round (waiting here polled squeue
            % continuously while the job was pending)
            if job_status >= 0 || (masterCompute && lastFile)
                return;
            end
    
            % If there is no job, submit a job
            if job_status == -1         
                if ~isempty(memAllocate) && cpusPerTask * memPerCPU < memAllocate
                    cpusPerTask = ceil(memAllocate / memPerCPU);
                end
                cpusPerTask = min(maxCPUNum, max(minCPUNum, cpusPerTask));
                time_str = '';
                if ~isempty(jobTimeLimit) && ~(contains(SlurmParam, ' -t ') || contains(SlurmParam, '--time'))
                    % only round to minutes (minimum 1 minute);
                    jobTimeLimit = max(jobTimeLimit, 1 / 60);
                    h = floor(jobTimeLimit);
                    m = round((jobTimeLimit - h) * 60);
                    time_str = sprintf(' -t %d:%d:00 ', h, m);
                end
                % with other places (SlurmAlternatives): the first with room
                % for the job now, with its shape there (cpus_k: the CPUs a
                % config with only that place gives it; see slurm_alternatives)
                if isempty(places) && ~isempty(SlurmAlternatives)
                    places = slurm_alternatives('places', struct('SlurmParam', pr.SlurmParam, 'SlurmConstraint', ...
                        pr.SlurmConstraint, 'memPerCPU', memPerCPU, 'maxCPUNum', maxCPUNum, 'minCPUNum', minCPUNum, ...
                        'gpuHeadroom', gpuHeadroom), SlurmAlternatives);
                    cpus_k = arrayfun(@(p) min(p.maxCPUNum, max(p.minCPUNum, max(pr.cpusPerTask, ...
                        ceil(max([memAllocate, 0]) / p.memPerCPU)))), places);
                    cpus_k(1) = cpusPerTask;
                    timeStrs = slurm_alternatives('time_strs', places, pr.jobTimeLimit, ' -t %d:%d:00 ', time_str);
                end
                k = 1;
                if numel(places) > 1
                    [k, cpusPerTask, SlurmParam] = slurm_alternatives('choose', places, cpus_k, timeStrs);
                    SlurmConstraint = places(k).SlurmConstraint;
                    time_str = timeStrs{k};
                end
                
                % the metrics folder of this MATLAB's single jobs ('' without
                % metrics; then a task's own metrics settings don't reach its jobs)
                runDir = '';
                if taskMetrics || ~isempty(getenv('PETAKIT5D_TASK_METRICS'))
                    runDir = task_metrics_run('session_dir', fileparts(jobLogFname), taskMetrics);
                end
                isMcc = ismcc || isdeployed || mccMode;
                if isMcc
                    % the compiled master must accept the options of the task
                    check_mcc_master_build(MCCMasterStr);
                    [func_name, task_var_str] = convert_function_string_to_mcc_string(funcStr);
                    % handle for nested "
                    var_str = strrep(task_var_str, '"', '\"');
                    % a MATLAB Runtime cache per job (jobs starting at the same
                    % time with a shared cache sometimes fail to start); the $ is
                    % escaped for the job's shell, as the command is in --wrap="..."
                    mcrCacheStr = sprintf('%s/mcr_cache_\\${SLURM_JOB_ID:-0}', MCRCacheRoot);
                    % the task's metrics environment
                    task_envStr = '';
                    if ~isempty(runDir)
                        task_envStr = task_metrics_run('env', runDir, task_id, false);
                    end
                    envStr = strrep(task_envStr, '"', '\"');
                    process_cmd = sprintf('%sMCR_CACHE_ROOT=%s %s %s %s %s \n', envStr, mcrCacheStr, MCCMasterStr, MCRParam, func_name, var_str);
                else
                    matlab_setup_str = 'setup([])';
                    % the task between its metrics calls
                    task_str = funcStr;
                    if ~isempty(runDir)
                        task_str = task_metrics_run('matlab_task', runDir, task_id, funcStr, false);
                    end
                    matlab_cmd = sprintf('%s;t0_=tic;%s;toc(t0_)', matlab_setup_str, task_str);
                    process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                end

                % task runners (Linux, MATLAB jobs or a compiled master built
                % with them): a first submission of a task joins the pool of
                % this MATLAB's single jobs of the same kind (function,
                % resources, master or MATLAB launch, job log folder; see
                % task_runner_run), made with the second task of the kind;
                % its job runs the task with its process of before, and after
                % it unstarted tasks of the pool's other jobs, or it ends if
                % another job's runner has its task. The first task of a kind
                % (a session with one pays nothing) and resubmissions run as
                % before
                poolDir = fileparts(jobLogFname);
                if isempty(poolTask) && taskRunner && job_id <= 0 && isunix && ~ismac && strcmpi(language, 'matlab') ...
                        && startsWith(poolDir, '/') && isfolder(poolDir)
                    k = find(funcStr == '(', 1);
                    if isempty(k)
                        k = numel(funcStr) + 1;
                    end
                    if isMcc
                        launch = sprintf('mcc|%s|%s|%s', MCCMasterStr, MCRParam, MCRCacheRoot);
                    else
                        launch = sprintf('matlab|%s', MatlabLaunchStr);
                    end
                    poolKey = sprintf('%s|%s|%d|%s|%s|%s|%s', strtrim(funcStr(1 : k - 1)), launch, cpusPerTask, ...
                        SlurmParam, SlurmConstraint, time_str, poolDir);
                    if ~any(strcmp(kindsSeen, poolKey))
                        kindsSeen{end + 1} = poolKey;
                    elseif ~isMcc || task_runner_run('supported', MCCMasterStr)
                        single_job_pools(true);
                        if isMcc
                            % the command of before, with the job's shell's view of $
                            task_line = char(sprintf('%sMCR_CACHE_ROOT=%s/mcr_cache_${SLURM_JOB_ID:-0} %s %s %s %s', task_envStr, ...
                                MCRCacheRoot, MCCMasterStr, MCRParam, func_name, task_var_str));
                            poolTask = task_runner_run('single_add', funcStr, poolKey, poolDir, task_line, MCCMasterStr, '');
                        else
                            % the stand-in master gets the task string (single-quoted)
                            task_line = sprintf('"$D/argdump/run_mccMaster.sh" x ''%s''', strrep(task_str, '''', '''\'''''));
                            poolTask = task_runner_run('single_add', funcStr, poolKey, poolDir, task_line, '', MatlabLaunchStr);
                        end
                    end
                end
                if ~isempty(poolTask)
                    process_cmd = task_runner_run('single_cmd', funcStr);
                end
        
                cmd = sprintf('sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s %s --wrap="date; echo Matlab command:  \\\"%s\\\"; %s"', ...
                    rem(task_id, 5000), jobLogFname, jobErrorFname, cpusPerTask, SlurmParam, ...
                    SlurmConstraint, time_str, funcStr, process_cmd);
                disp(cmd);
                
                [status, cmdout] = system(cmd, '-echo');
        
                job_out_info = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                if ~isempty(job_out_info) && ~isempty(job_out_info{1})
                    job_id = str2double(job_out_info{1}{1});
                    submit_status = 1;
                    if numel(places) > 1
                        slurm_alternatives('submitted', job_id, array_id, places, cpus_k, k, timeStrs);
                    end
                    if ~isempty(poolTask)
                        task_runner_run('single_job', funcStr, job_id, array_id);
                    end
                    if ~isempty(runDir)
                        task_metrics_run('record', runDir, 'submit', struct('task', task_id, 'jobId', job_id, ...
                            'functionStr', funcStr, 'cpusPerTask', cpusPerTask, 'memAllocate', memAllocate, ...
                            'memPerCPU', memPerCPU, 'mccMode', isMcc, 'SlurmParam', SlurmParam));
                    end
                else
                    fprintf('Job submission failed, try again!\n')
                end
                trial_counter = trial_counter + 1;
            end
    end
end

% a pool task whose job couldn't be submitted: the master takes it back (the
% caller submits it again)
if submit_status == 0 && ~isempty(poolTask)
    task_runner_run('single_job', funcStr, [], []);
end

end

