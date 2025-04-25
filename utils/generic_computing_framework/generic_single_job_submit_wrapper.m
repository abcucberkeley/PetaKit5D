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
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2023a; matlab -nodisplay -nosplash -nodesktop -nojvm -r', @ischar);
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

ip.parse(funcStr, job_id, task_id, varargin{:});

% move to the root path
funcFn = which(mfilename);
paths = split(funcFn, 'PetaKit5D');
cd(paths{1});
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
    for i = 1 : numel(conf_keys)
        pr.(conf_keys{i}) = confData.(conf_keys{i});
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

if exist(jobLogDir, 'dir')
    dataPath = tempdir();
    if ispc
        dataPath = strrep(dataPath, '\', '/');
    end
    [parseCluster, jobLogFname, jobErrorFname] = checkSlurmCluster(dataPath, jobLogDir);
end

submit_status = 0;
trial_counter = 0;
while submit_status == 0 && trial_counter < maxTrialNum
    switch clusterType
        case 'slurm'
            array_id = rem(task_id, 5000);
            job_status = check_slurm_job_status(job_id, array_id);
        
            % kill the last pending job and use master node do the computing.
            if job_status == 0.5 && (masterCompute && lastFile)
                system(sprintf('scancel %d_%d', job_id, array_id), '-echo');
            end
    
            % if the job is still running, skip it. 
            if job_status == 1 || (masterCompute && lastFile)
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
                
                if ismcc || isdeployed || mccMode
                    [func_name, var_str] = convert_function_string_to_mcc_string(funcStr);
                    % handle for nested "
                    var_str = strrep(var_str, '"', '\"');                
                    process_cmd = sprintf('MCR_CACHE_ROOT=%s %s %s %s %s \n', MCRCacheRoot, MCCMasterStr, MCRParam, func_name, var_str);
                else
                    matlab_setup_str = 'setup([])';
                    matlab_cmd = sprintf('%s;t0_=tic;%s;toc(t0_)', matlab_setup_str, funcStr);
                    process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
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
                else
                    fprintf('Job submission failed, try again!\n')
                end
                trial_counter = trial_counter + 1;
            end
    end
end

end

