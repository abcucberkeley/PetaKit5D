classdef SlurmClusterTest < pktest.TestCase
    % the job framework in Slurm cluster mode (parseCluster=true): tasks are
    % submitted with sbatch, polled with squeue, retried with more CPUs and
    % cancelled at the end. These tests submit real jobs, so they only run
    % when PETAKIT5D_TEST_SLURM_CONFIG is set (Linux and macOS only):
    %
    %   PETAKIT5D_TEST_SLURM_CONFIG  : JSON cluster configuration passed to
    %       the framework as configFile (e.g. SlurmParam, MatlabLaunchStr,
    %       maxCPUNum, memPerCPU), see tests/cluster/github_actions_slurm.json.
    %       It must not set jobLogDir or SlurmConstraint, which the tests use.
    %   PETAKIT5D_TEST_SLURM_SCRATCH : folder on storage shared with the
    %       compute nodes for test data and job logs (default: tempdir, which
    %       only works when jobs run on this machine, as on the single-node
    %       CI cluster)
    %   PETAKIT5D_TEST_SLURM_MATLAB  : 'true' to also run the tests whose jobs
    %       start MATLAB (needs MATLAB and a license on the compute nodes)
    %
    % Tasks are small bash scripts that record the Slurm environment of the
    % job that ran them (array job id, CPUs per task, time limit).

    properties (SetAccess = private)
        ConfigFile
        ScratchRoot
        MaxCPUNum = 24
        JobName
        WatchdogFired = false
    end

    methods (TestClassSetup)
        function requireSlurm(testCase)
            testCase.assumeFalse(ispc, 'Slurm cluster mode is not supported on Windows.');
            configFile = getenv('PETAKIT5D_TEST_SLURM_CONFIG');
            testCase.assumeNotEmpty(configFile, ['Set PETAKIT5D_TEST_SLURM_CONFIG to a cluster ', ...
                'configuration file to run the Slurm tests (they submit real jobs).']);
            testCase.assertTrue(isfile(configFile), sprintf('Cluster configuration file %s does not exist.', configFile));

            testCase.applyFixture(pktest.RealSlurmFixture);
            [status, cmdout] = system('sinfo');
            testCase.assertEqual(status, 0, sprintf('PETAKIT5D_TEST_SLURM_CONFIG is set but sinfo fails:\n%s', cmdout));

            testCase.ConfigFile = configFile;
            conf = jsondecode(fileread(configFile));
            if isfield(conf, 'maxCPUNum')
                testCase.MaxCPUNum = conf.maxCPUNum;
            end
            scratch = getenv('PETAKIT5D_TEST_SLURM_SCRATCH');
            if isempty(scratch)
                scratch = tempdir;
            end
            testCase.ScratchRoot = [strip(strrep(scratch, '\', '/'), 'right', '/'), '/'];
        end
    end

    methods (TestMethodSetup)
        function prepareJobs(testCase)
            % the framework caches the configuration file in persistent variables
            clear generic_computing_frameworks_wrapper
            uuid = get_uuid();
            testCase.JobName = ['pk5dtest_', uuid(1 : 8)];
            testCase.WatchdogFired = false;
            % cancel anything a failed test left in the queue
            testCase.addTeardown(@() system(sprintf('scancel --name=%s', testCase.JobName)));
        end
    end

    methods (Test, TestTags = {'Cluster', 'Unit'})
        function jobStatusFollowsSqueue(testCase)
            % check_batch_slurm_jobs_status parses the default squeue output
            dataPath = testCase.createSharedFolder();
            conf = jsondecode(fileread(testCase.ConfigFile));
            [status, cmdout] = system(sprintf(['sbatch --hold --array=3 -J %s -o %sjob_%%A_%%a.out %s -t 0:05:00 ', ...
                '--wrap="sleep 120"'], testCase.JobName, dataPath, conf.SlurmParam));
            testCase.assertEqual(status, 0, cmdout);
            jobId = str2double(regexp(cmdout, 'Submitted batch job (\d+)', 'tokens', 'once'));

            testCase.verifyEqual(check_batch_slurm_jobs_status(jobId, 3), 0, 'A held job should be pending (0).');
            system(sprintf('scontrol release %d_3', jobId));
            testCase.verifyEqual(testCase.waitForJobStatus(jobId, 3, 1, 120), 1, 'The released job should be running (1).');
            system(sprintf('scancel %d', jobId));
            testCase.verifyEqual(testCase.waitForJobStatus(jobId, 3, -1, 60), -1, 'A cancelled job should be gone (-1).');
        end
    end

    methods (Test, TestTags = {'Cluster', 'Integration'})
        function tasksRunAsSlurmJobs(testCase)
            dataPath = testCase.createSharedFolder();
            [inputs, outputs, funcStrs] = testCase.bashTasks(dataPath, 3);
            isDone = testCase.runFramework(dataPath, inputs, outputs, funcStrs);

            testCase.verifyEqual(isDone(:)', true(1, 3));
            jobs = testCase.readTaskField(outputs, 'array_job');
            testCase.verifyFalse(any(cellfun(@isempty, jobs)), 'Some tasks did not run as Slurm jobs.');
            testCase.verifyEqual(numel(unique(jobs)), 3, 'Each task should run as its own job.');
            testCase.verifyNotEmpty(dir([dataPath, 'job_logs/job_*.out']), 'No job logs in jobLogDir.');
        end

        function tasksInABatchShareAJob(testCase)
            dataPath = testCase.createSharedFolder();
            [inputs, outputs, funcStrs] = testCase.bashTasks(dataPath, 4);
            isDone = testCase.runFramework(dataPath, inputs, outputs, funcStrs, 'taskBatchNum', 2);

            testCase.verifyEqual(isDone(:)', true(1, 4));
            jobs = testCase.readTaskField(outputs, 'array_job');
            testCase.verifyFalse(any(cellfun(@isempty, jobs)), 'Some tasks did not run as Slurm jobs.');
            testCase.verifyEqual(jobs{1}, jobs{2}, 'Tasks 1 and 2 should share a job.');
            testCase.verifyEqual(jobs{3}, jobs{4}, 'Tasks 3 and 4 should share a job.');
            testCase.verifyNotEqual(jobs{1}, jobs{3}, 'Each batch should run as its own job.');
        end

        function existingOutputsAreNotSubmitted(testCase)
            dataPath = testCase.createSharedFolder();
            [inputs, outputs, funcStrs] = testCase.bashTasks(dataPath, 3);
            fid = fopen(outputs{2}, 'w');
            fprintf(fid, 'existing\n');
            fclose(fid);
            isDone = testCase.runFramework(dataPath, inputs, outputs, funcStrs);

            testCase.verifyEqual(isDone(:)', true(1, 3));
            testCase.verifyEqual(strtrim(fileread(outputs{2})), 'existing', 'A task with an existing output was run.');
            jobs = testCase.readTaskField(outputs([1, 3]), 'array_job');
            testCase.verifyFalse(any(cellfun(@isempty, jobs)), 'The other tasks did not run as Slurm jobs.');
        end

        function masterComputesPartOfTheTasks(testCase)
            % with masterCompute the MATLAB session runs a batch itself
            % instead of waiting for its job
            dataPath = testCase.createSharedFolder();
            [inputs, outputs, funcStrs] = testCase.bashTasks(dataPath, 3);
            isDone = testCase.runFramework(dataPath, inputs, outputs, funcStrs, 'masterCompute', true);

            testCase.verifyEqual(isDone(:)', true(1, 3));
            jobs = testCase.readTaskField(outputs, 'array_job');
            testCase.verifyTrue(any(cellfun(@isempty, jobs)), 'No task was run by the MATLAB session itself.');
        end

        function requestedResourcesReachTheJob(testCase)
            testCase.assumeGreaterThanOrEqual(testCase.MaxCPUNum, 2, 'The cluster configuration allows only 1 CPU per task.');
            dataPath = testCase.createSharedFolder();
            [inputs, outputs, funcStrs] = testCase.bashTasks(dataPath, 2);
            isDone = testCase.runFramework(dataPath, inputs, outputs, funcStrs, 'cpusPerTask', 2, 'jobTimeLimit', 0.25);

            testCase.verifyEqual(isDone(:)', true(1, 2));
            testCase.verifyEqual(testCase.readTaskField(outputs, 'cpus'), {'2', '2'});
            testCase.verifyEqual(testCase.readTaskField(outputs, 'time_limit'), {'00:15:00', '00:15:00'});
        end

        function finalOutputEndsTheRun(testCase)
            % the first task writes the final output; the others would run
            % for minutes, but the framework stops once the final output exists
            dataPath = testCase.createSharedFolder();
            [inputs, outputs, funcStrs] = testCase.bashTasks(dataPath, 3, 'sleepSeconds', [0, 600, 600]);
            finalFn = [dataPath, 'final.txt'];
            funcStrs{1} = sprintf('%s && cp %s %s', funcStrs{1}, outputs{1}, finalFn);
            t0 = tic;
            isDone = testCase.runFramework(dataPath, inputs, outputs, funcStrs, 'finalOutFullpath', finalFn);

            testCase.verifyLessThan(toc(t0), 300, 'The framework kept waiting after the final output existed.');
            testCase.verifyEqual(isDone(:)', true(1, 3));
        end

        function finalOutputCancelsRemainingJobs(testCase)
            % once the final output exists, the framework returns and
            % cancels the jobs still running (it used to return before the
            % scancel, leaving them running)
            dataPath = testCase.createSharedFolder();
            [inputs, outputs, funcStrs] = testCase.bashTasks(dataPath, 3, 'sleepSeconds', [0, 600, 600]);
            finalFn = [dataPath, 'final.txt'];
            funcStrs{1} = sprintf('%s && cp %s %s', funcStrs{1}, outputs{1}, finalFn);
            testCase.runFramework(dataPath, inputs, outputs, funcStrs, 'finalOutFullpath', finalFn);

            % cancelled jobs can take a few seconds to leave the queue
            t0 = tic;
            [~, cmdout] = system(sprintf('squeue -h --name=%s', testCase.JobName));
            while ~isempty(strtrim(cmdout)) && toc(t0) < 60
                pause(2);
                [~, cmdout] = system(sprintf('squeue -h --name=%s', testCase.JobName));
            end
            testCase.verifyEmpty(strtrim(cmdout), sprintf('Jobs still queued after the final output existed:\n%s', cmdout));
        end
    end

    methods (Test, TestTags = {'Cluster', 'Integration', 'Slow'})
        function failedTaskIsRetriedWithMoreCPUs(testCase)
            % a job that ends without its output is resubmitted (after a
            % grace period) with cpusPerTask * (trial + 1) CPUs, up to maxTrialNum
            testCase.assumeGreaterThanOrEqual(testCase.MaxCPUNum, 2, 'The cluster configuration allows only 1 CPU per task.');
            dataPath = testCase.createSharedFolder();
            [inputs, outputs, funcStrs] = testCase.bashTasks(dataPath, 2, 'fail', [false, true]);
            isDone = testCase.runFramework(dataPath, inputs, outputs, funcStrs, 'maxTrialNum', 2);

            testCase.verifyEqual(isDone(:)', [true, false]);
            attempts = splitlines(strtrim(fileread([dataPath, 'attempts_2.txt'])));
            testCase.verifyEqual(attempts(:)', {'cpus=1', 'cpus=2'}, 'The failed task should be tried twice, with 1 and then 2 CPUs.');
        end
    end

    methods (Test, TestTags = {'Cluster', 'Integration'})
        function matlabTasksRunAsJobs(testCase)
            % jobs start MATLAB with MatlabLaunchStr, run setup.m from the
            % repository root and then the task string
            testCase.assumeMatlabJobs();
            dataPath = testCase.createSharedFolder();
            logFn = [dataPath, 'calls.log'];
            inputs = cell(1, 2);
            outputs = cell(1, 2);
            funcStrs = cell(1, 2);
            for i = 1 : 2
                inputs{i} = sprintf('%sin_%d.txt', dataPath, i);
                outputs{i} = sprintf('%sout_%d.txt', dataPath, i);
                fid = fopen(inputs{i}, 'w');
                fprintf(fid, 'task %d\n', i);
                fclose(fid);
                funcStrs{i} = sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'')', inputs{i}, outputs{i}, logFn);
            end
            isDone = testCase.runFramework(dataPath, inputs, outputs, funcStrs, 'language', 'matlab');

            testCase.verifyEqual(isDone(:)', true(1, 2));
            for i = 1 : 2
                testCase.verifyEqual(fileread(outputs{i}), fileread(inputs{i}));
            end
            testCase.verifyNumElements(splitlines(strtrim(fileread(logFn))), 2, 'Each task should run once.');
            % the outputs alone would also appear if the framework fell back
            % to running the tasks in this session
            testCase.verifyGreaterThanOrEqual(testCase.countJobLogs({[dataPath, 'job_logs']}, dataPath), 2, ...
                'The tasks did not run as Slurm jobs.');
        end

        function deskewRotateWrapperRunsAsJobs(testCase)
            % a data wrapper end to end on the cluster gives the same result
            % as the in-memory function
            testCase.assumeMatlabJobs();
            dataPath = testCase.createSharedFolder();
            fsn = 'Scan_Iter_0000_0000_Cam%s_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t';
            cams = {'A', 'B'};
            frames = cell(1, 2);
            for i = 1 : 2
                frames{i} = uint16(pktest.beads([32, 48, 40], [16, 20, 12; 10, 30, 28], 2, 1000, 100) + 50 * i);
                writetiff(frames{i}, [dataPath, sprintf(fsn, cams{i}), '.tif']);
            end
            XR_deskew_rotate_data_wrapper({dataPath}, 'deskew', true, 'rotate', true, 'DSRCombined', true, ...
                'xyPixelSize', 0.108, 'dz', 0.2, 'skewAngle', 32.45, 'reverse', true, 'channelPatterns', {'CamA', 'CamB'}, ...
                'save16bit', true, 'parseCluster', true, 'masterCompute', false, 'configFile', testCase.ConfigFile, ...
                'jobLogDir', [dataPath, 'job_logs']);

            for i = 1 : 2
                expected = deskewRotateFrame3D(frames{i}, 32.45, 0.2, 0.108, 'reverse', true, 'save16bit', true);
                testCase.verifyEqual(readtiff([dataPath, 'DSR/', sprintf(fsn, cams{i}), '.tif']), expected);
            end
            % the job logs go to jobLogDir (the wrapper used to drop it)
            testCase.verifyGreaterThanOrEqual(testCase.countJobLogs({[dataPath, 'job_logs']}, dataPath), 2, ...
                'The frames were not processed as Slurm jobs, or their logs are not in jobLogDir.');
        end
    end

    methods (Access = private)
        function [dataPath] = createSharedFolder(testCase)
            % on storage the compute nodes can reach, removed after the
            % test once its jobs are gone (they may still write job logs)
            uuid = get_uuid();
            dataPath = sprintf('%spk5d_slurm_test_%s/', testCase.ScratchRoot, uuid(1 : 8));
            mkdir([dataPath, 'job_logs']);
            testCase.addTeardown(@() testCase.removeFolderAfterJobs(dataPath));
        end

        function [] = removeFolderAfterJobs(testCase, dataPath)
            system(sprintf('scancel --name=%s', testCase.JobName));
            t0 = tic;
            [~, cmdout] = system(sprintf('squeue -h --name=%s', testCase.JobName));
            while ~isempty(strtrim(cmdout)) && toc(t0) < 60
                pause(2);
                [~, cmdout] = system(sprintf('squeue -h --name=%s', testCase.JobName));
            end
            testCase.removeFolder(dataPath);
        end

        function [inputs, outputs, funcStrs] = bashTasks(~, dataPath, nTasks, varargin)
            % each task is "bash task_i.sh <output>"; the script writes the
            % Slurm environment of its job to the output (via a temporary
            % file, so the output appears only when complete), or records
            % the attempt and exits without output if it should fail
            ip = inputParser;
            ip.addParameter('sleepSeconds', zeros(1, nTasks), @isnumeric);
            ip.addParameter('fail', false(1, nTasks), @islogical);
            ip.parse(varargin{:});
            pr = ip.Results;

            inputs = cell(1, nTasks);
            outputs = cell(1, nTasks);
            funcStrs = cell(1, nTasks);
            for i = 1 : nTasks
                inputs{i} = sprintf('%sin_%d.txt', dataPath, i);
                outputs{i} = sprintf('%sout_%d.txt', dataPath, i);
                scriptFn = sprintf('%stask_%d.sh', dataPath, i);
                fclose(fopen(inputs{i}, 'w'));
                if pr.fail(i)
                    body = sprintf('echo "cpus=${SLURM_CPUS_PER_TASK:-}" >> %sattempts_%d.txt\nexit 1\n', dataPath, i);
                else
                    body = sprintf(['sleep %d\n', ...
                        'limit=$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | grep -o "TimeLimit=[^ ]*" | cut -d= -f2)\n', ...
                        'printf "array_job=%%s\\ncpus=%%s\\ntime_limit=%%s\\n" "${SLURM_ARRAY_JOB_ID:-}" "${SLURM_CPUS_PER_TASK:-}" "$limit" > "$1.part"\n', ...
                        'mv "$1.part" "$1"\n'], pr.sleepSeconds(i));
                end
                fid = fopen(scriptFn, 'w');
                fprintf(fid, '#!/bin/bash\n%s', body);
                fclose(fid);
                funcStrs{i} = sprintf('bash %s %s', scriptFn, outputs{i});
            end
        end

        function [isDone] = runFramework(testCase, dataPath, inputs, outputs, funcStrs, varargin)
            testCase.startWatchdog(10);
            isDone = generic_computing_frameworks_wrapper(inputs, outputs, funcStrs, 'parseCluster', true, ...
                'masterCompute', false, 'language', 'bash', 'jobLogDir', [dataPath, 'job_logs'], ...
                'ConfigFile', testCase.ConfigFile, 'maxTrialNum', 2, 'jobTimeLimit', 0.25, 'queryInterval', 2, ...
                'SlurmConstraint', ['-J ', testCase.JobName], varargin{:});
            testCase.verifyFalse(testCase.WatchdogFired, 'Jobs did not finish within 10 minutes and were cancelled.');
        end

        function [] = startWatchdog(testCase, minutes)
            % if jobs hang (e.g., pending forever), keep cancelling them so
            % the framework runs out of trials and returns
            t = timer('StartDelay', minutes * 60, 'Period', 20, 'ExecutionMode', 'fixedRate', ...
                'TimerFcn', @(~, ~) testCase.cancelJobs());
            start(t);
            testCase.addTeardown(@() testCase.deleteTimer(t));
        end

        function [] = cancelJobs(testCase)
            testCase.WatchdogFired = true;
            system(sprintf('scancel --name=%s', testCase.JobName));
        end

        function [n] = countJobLogs(~, logDirs, dataPath)
            % number of Slurm job logs in logDirs whose command mentions
            % dataPath (each job echoes its command to its log)
            n = 0;
            for i = 1 : numel(logDirs)
                dir_info = dir(fullfile(logDirs{i}, 'job_*.out'));
                for j = 1 : numel(dir_info)
                    if contains(fileread(fullfile(dir_info(j).folder, dir_info(j).name)), dataPath)
                        n = n + 1;
                    end
                end
            end
        end

        function [values] = readTaskField(~, outputs, key)
            % value of key=value lines in each task output ('' if missing)
            values = cell(1, numel(outputs));
            for i = 1 : numel(outputs)
                values{i} = '';
                if isfile(outputs{i})
                    % in MATLAB regexp, . also matches newlines
                    tok = regexp(fileread(outputs{i}), ['(?m)^', key, '=([^\r\n]*)$'], 'tokens', 'once');
                    if ~isempty(tok)
                        values{i} = strtrim(tok{1});
                    end
                end
            end
        end

        function [status] = waitForJobStatus(~, jobId, arrayId, target, timeoutSeconds)
            t0 = tic;
            status = check_batch_slurm_jobs_status(jobId, arrayId);
            while status ~= target && toc(t0) < timeoutSeconds
                pause(2);
                status = check_batch_slurm_jobs_status(jobId, arrayId);
            end
        end

        function [] = assumeMatlabJobs(testCase)
            testCase.assumeTrue(any(strcmpi(getenv('PETAKIT5D_TEST_SLURM_MATLAB'), {'true', '1'})), ...
                'Set PETAKIT5D_TEST_SLURM_MATLAB=true to run tests whose jobs start MATLAB.');
        end
    end

    methods (Static, Access = private)
        function [] = removeFolder(dataPath)
            % network file systems can take a moment to release files that
            % jobs just closed, so retry a few times
            for i = 1 : 5
                if ~exist(dataPath, 'dir')
                    return;
                end
                try
                    rmdir(dataPath, 's');
                catch
                    pause(2);
                end
            end
        end

        function [] = deleteTimer(t)
            stop(t);
            delete(t);
        end
    end
end
