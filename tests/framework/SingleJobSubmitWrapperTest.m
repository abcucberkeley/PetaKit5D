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
