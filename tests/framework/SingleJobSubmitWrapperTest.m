classdef SingleJobSubmitWrapperTest < pktest.TestCase
    % generic_single_job_submit_wrapper builds the sbatch command for one
    % task. pktest.NoSlurmFixture (shared by all tests) makes squeue and
    % sbatch fail, so nothing is submitted; the command is read from the
    % console output. The job-state tests use pktest.FakeSlurmFixture.

    properties (TestParameter)
        % name-value arguments and the expected --cpus-per-task
        cpuCase = struct( ...
            'cpusAboveMax', {{{'cpusPerTask', 40}, 24}}, ...
            'minCPUNum', {{{'minCPUNum', 4}, 4}}, ...
            'memAllocate', {{{'memAllocate', 100}, 5}}, ...
            'memAllocateAboveMax', {{{'memAllocate', 1000}, 24}})
    end

    methods (Test, TestTags = {'Unit'})
        function cpusPerTaskIsClampedToLimits(testCase, cpuCase)
            % fb6a98b regression: user-defined cpusPerTask/minCPUNum were
            % only clamped to [minCPUNum, maxCPUNum] when memAllocate was set
            cmd = testCase.submitCommand(cpuCase{1}{:});
            cpus = regexp(cmd, '--cpus-per-task=(\d+)', 'tokens', 'once');
            testCase.assertNotEmpty(cpus, 'No sbatch command was printed.');
            testCase.verifyEqual(str2double(cpus{1}), cpuCase{2});
        end

        function commandWrapsTheTaskString(testCase)
            cmd = testCase.submitCommand('SlurmParam', '-p test_partition -n1');
            testCase.verifySubstring(cmd, 'sbatch --array=7 ');
            testCase.verifySubstring(cmd, '-p test_partition -n1');
            testCase.verifySubstring(cmd, 'disp(12345)');
            testCase.verifySubstring(cmd, 'matlab/r2024b');
            % MATLAB's non-interactive mode (with -r it starts the MathWorks
            % Service Host and was very slow to exit on the cluster)
            testCase.verifySubstring(cmd, '-nojvm -batch ');
        end

        function configFileKeysIgnoreCase(testCase)
            % the cluster configs in use spell MemPerCPU; config keys used to
            % be matched case-sensitively, so it was silently ignored
            configFn = [testCase.createTempFolder(), 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(struct('MemPerCPU', 5)));
            fclose(fid);
            cmd = testCase.submitCommand('ConfigFile', configFn, 'memAllocate', 50);
            testCase.verifySubstring(cmd, '--cpus-per-task=10');
        end

        function emptyJobLogPathsInAConfigFileAreTheDefaults(testCase)
            % empty job log paths in a cluster config ("jobLogDir": "",
            % "jobLogFname": "") gave the job "-o  -e", so its log went to a
            % file named -e; they are now the default ones
            configFn = [testCase.createTempFolder(), 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(struct('jobLogDir', '', 'jobLogFname', '', 'jobErrorFname', '')));
            fclose(fid);
            cmd = testCase.submitCommand('ConfigFile', configFn);
            testCase.assertNotEmpty(cmd, 'No sbatch command was printed.');
            paths = regexp(cmd, ' -o (\S+) -e (\S+) ', 'tokens', 'once');
            testCase.assertNotEmpty(paths, cmd);
            testCase.verifyEqual(paths, {'../job_logs/job_%A_%a.out', '../job_logs/job_%A_%a.err'});
        end

        function pendingJobReturnsWithoutWaiting(testCase)
            % the wrapper used to call squeue in a loop until a pending job
            % started (74-94 times a second on the cluster); callers check
            % the job again in their next sweep, as for running jobs
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('123_7', {'PD', 'PD', 'PD', 'R'});
            [jobId, jobStatus, submitStatus] = testCase.callWrapper(123, 7);
            testCase.verifyEqual([jobId, jobStatus, submitStatus], [123, 0, 0]);
            testCase.verifyEqual(slurm.callCount('squeue'), 1);
            testCase.verifyEqual(slurm.callCount('sbatch'), 0);
        end

        function runningJobReturnsWithoutWaiting(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('123_7', {'R'});
            [jobId, jobStatus, submitStatus] = testCase.callWrapper(123, 7);
            testCase.verifyEqual([jobId, jobStatus, submitStatus], [123, 1, 0]);
            testCase.verifyEqual(slurm.callCount('squeue'), 1);
        end

        function jobNotInTheQueueIsSubmitted(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [jobId, ~, submitStatus] = testCase.callWrapper(-1, 7);
            testCase.verifyEqual([jobId, submitStatus], [1000, 1]);
            testCase.verifyEqual(slurm.callCount('sbatch'), 1);
        end

        function taskWithoutAJobIsSubmittedWithoutAQueueCheck(testCase)
            % squeue (8-10 ms a call on the cluster) ran also for a task
            % without a job (job id -1), where it could only say so; a job
            % is still checked
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [jobId, ~, submitStatus] = testCase.callWrapper(-1, 7);
            testCase.verifyEqual([jobId, submitStatus], [1000, 1]);
            testCase.verifyEqual(slurm.callCount('squeue'), 0);
            slurm.setJobStates('1000_7', {'R'});
            [~, jobStatus] = testCase.callWrapper(1000, 7);
            testCase.verifyEqual(jobStatus, 1);
            testCase.verifyEqual(slurm.calls('squeue'), {'-j 1000_7'});
        end

        function slurmIsCheckedOncePerJobLogFolder(testCase)
            % callers call the wrapper for each job in each round, and each
            % call ran sinfo and listed the job log folder (10 ms or more);
            % both are kept for the folder while Slurm is there, with the
            % same job log paths
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            clear generic_single_job_submit_wrapper task_runner_run single_job_pools
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            logDir = testCase.createTempFolder();
            otherLogDir = testCase.createTempFolder();
            testCase.callWrapper(-1, 7, 'jobLogDir', logDir);
            testCase.callWrapper(-1, 8, 'jobLogDir', logDir);
            testCase.verifyEqual(nnz(startsWith(slurm.history(), 'sinfo')), 1);
            testCase.callWrapper(-1, 9, 'jobLogDir', otherLogDir);
            testCase.verifyEqual(nnz(startsWith(slurm.history(), 'sinfo')), 2);
            logs = regexp(slurm.calls('sbatch'), '-o (\S+) -e (\S+)', 'tokens', 'once');
            testCase.assertNumElements(logs, 3);
            testCase.verifyEqual(logs{2}, logs{1});
            testCase.verifyEqual(logs{1}, {[logDir, 'job_%A_%a.out'], [logDir, 'job_%A_%a.err']});
            testCase.verifyEqual(logs{3}, {[otherLogDir, 'job_%A_%a.out'], [otherLogDir, 'job_%A_%a.err']});
        end

        function sessionWithoutRunnersDoesNotLoadThem(testCase)
            % with runners off, the wrapper doesn't look for pool tasks
            % before a task of the session joined a pool (the runner code's
            % first call took 60-80 ms)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            clear generic_single_job_submit_wrapper task_runner_run single_job_pools
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            logDir = testCase.createTempFolder();
            testCase.callWrapper(-1, 7, 'jobLogDir', logDir, 'taskRunner', false);
            testCase.callWrapper(-1, 8, 'jobLogDir', logDir, 'taskRunner', false);
            testCase.verifyEqual(slurm.callCount('sbatch'), 2);
            testCase.verifyFalse(ismember('task_runner_run', inmem), 'The runner code was loaded.');
        end

        function slurmIsCheckedInEachCallWithoutSlurm(testCase)
            % without Slurm (sinfo fails) nothing is kept: each call checks
            % again, as before
            testCase.assumeFalse(ispc, 'The sinfo stub needs a POSIX shell.');
            clear generic_single_job_submit_wrapper task_runner_run single_job_pools
            stubDir = testCase.createTempFolder();
            logFn = [stubDir, 'sinfo.log'];
            fid = fopen([stubDir, 'sinfo'], 'w');
            fprintf(fid, '#!/bin/sh\necho x >> "%s"\nexit 1\n', logFn);
            fclose(fid);
            fileattrib([stubDir, 'sinfo'], '+x');
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('PATH', ...
                [stubDir, pathsep, getenv('PATH')]));
            logDir = testCase.createTempFolder();
            testCase.callWrapper(-1, 7, 'jobLogDir', logDir);
            testCase.callWrapper(-1, 8, 'jobLogDir', logDir);
            testCase.verifyNumElements(splitlines(strip(fileread(logFn))), 2);
        end

        function compiledJobsUseOwnRuntimeCache(testCase)
            % single compiled jobs shared one MATLAB Runtime cache, so jobs
            % starting at the same time on a node sometimes failed to start
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            mkdir([dataPath, 'mcr/bin']);
            logFn = [dataPath, 'cache.log'];
            masterFn = [dataPath, 'run_mccMaster.sh'];
            fid = fopen(masterFn, 'w');
            fprintf(fid, '#!/bin/bash\necho "$MCR_CACHE_ROOT" > "%s"\n', logFn);
            fclose(fid);
            fileattrib(masterFn, '+x');
            configFn = [dataPath, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(struct('MCRParam', [dataPath, 'mcr/'], 'MCCMasterStr', masterFn, ...
                'MCRCacheRoot', [dataPath, 'cache'])));
            fclose(fid);
            evalc(['generic_single_job_submit_wrapper(''XR_crop_frame(''''in'''',''''out'''',[1,1,1,2,2,2],''''pad'''',false)'', ', ...
                '-1, 7, ''maxTrialNum'', 1, ''mccMode'', true, ''ConfigFile'', configFn)']);
            calls = slurm.calls('sbatch');
            testCase.assertNumElements(calls, 1);
            % run the job's command as Slurm would, in the job's shell
            jobFn = [dataPath, 'job.sh'];
            fid = fopen(jobFn, 'w');
            fprintf(fid, '%s\n', extractAfter(calls{1}, '--wrap='));
            fclose(fid);
            [status, out] = system(sprintf('env SLURM_JOB_ID=42 bash %s', jobFn));
            testCase.assertEqual(status, 0, out);
            testCase.verifyEqual(strip(fileread(logFn)), [dataPath, 'cache/mcr_cache_42']);
        end

        function pendingLastJobGoesToTheMaster(testCase)
            % with masterCompute the master cancels the last file's pending
            % job and computes it itself
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('123_7', {'PD'});
            [~, jobStatus, submitStatus] = testCase.callWrapper(123, 7, 'masterCompute', true, 'lastFile', true);
            testCase.verifyEqual([jobStatus, submitStatus], [0, 0]);
            testCase.verifyEqual(slurm.calls('scancel'), {'123_7'});
        end
    end

    methods (Access = private)
        function [cmd] = submitCommand(~, varargin)
            out = evalc('generic_single_job_submit_wrapper(''disp(12345)'', 1, 7, ''maxTrialNum'', 1, varargin{:})');
            cmd = regexp(out, 'sbatch [^\n]*', 'match', 'once');
        end

        function [jobId, jobStatus, submitStatus] = callWrapper(~, jobId, taskId, varargin)
            % the sbatch command is printed; keep it out of the test log
            [~, jobId, jobStatus, submitStatus] = evalc(['generic_single_job_submit_wrapper(''disp(1)'', jobId, ', ...
                'taskId, ''maxTrialNum'', 1, varargin{:})']);
        end
    end
end
