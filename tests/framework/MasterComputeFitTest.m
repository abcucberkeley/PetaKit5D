classdef MasterComputeFitTest < pktest.TestCase
    % jobs that coordinate jobs: a master in a Slurm job runs tasks itself
    % (masterCompute) only if they fit its memory (master_task_fits, in the
    % job framework), and the request of a coordinating job is what it
    % holds itself (coordinator_resources, MIP_zarr_merge_memory,
    % cluster_config_value). The framework tests submit to a fake Slurm.

    methods (TestMethodSetup)
        function clearConfigCache(testCase)
            % (config files are cached by name and modification time)
            clear generic_computing_frameworks_wrapper
        end
    end

    methods (Test, TestTags = {'Unit'})
        function masterInAJobRunsNoTaskThatDoesNotFitIt(testCase)
            % a master in a Slurm job ran the last task (or a single task)
            % itself whatever its estimate: a user's 4-CPU, 64 GB
            % pipeline job ran an Imaris conversion and was killed for
            % memory. Now tasks that don't fit its memory all run as jobs
            % (the fake jobs end at once without results)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.setSlurmJob('4242');
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 16));
            for n = [1, 3]
                slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
                [logFn, out] = testCase.runFramework(n, 'memAllocate', 40);
                testCase.verifyEmpty(testCase.readLog(logFn), sprintf('%d task(s): the master ran a task.', n));
                testCase.verifyEqual(slurm.callCount('sbatch'), n, sprintf('%d task(s) were not all submitted.', n));
                testCase.verifySubstring(out, 'do not fit this job');
            end
        end

        function masterInAJobRunsTasksThatFitIt(testCase)
            % with the memory of a task, the master runs the last one as
            % before, and submits the others
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.setSlurmJob('4242');
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 64));
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [logFn, ~, outFns] = testCase.runFramework(3, 'memAllocate', 40);
            testCase.verifyEqual(testCase.readLog(logFn), outFns(3));
            testCase.verifyEqual(slurm.callCount('sbatch'), 2);
        end

        function tasksWithoutAnEstimateNeedTheMemoryOfTheirJobs(testCase)
            % without memAllocate, a task's job has the memory of its CPUs
            % (4 x 10 GB here), which is what the master needs
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.setSlurmJob('4242');
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 16));
            cases = {4, 0, 2; 1, 1, 1};
            for i = 1 : size(cases, 1)
                slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
                logFn = testCase.runFramework(2, 'cpusPerTask', cases{i, 1}, 'memPerCPU', 10);
                testCase.verifyNumElements(testCase.readLog(logFn), cases{i, 2}, sprintf('%d CPUs per task.', cases{i, 1}));
                testCase.verifyEqual(slurm.callCount('sbatch'), cases{i, 3}, sprintf('%d CPUs per task.', cases{i, 1}));
            end
        end

        function masterOutsideAJobRunsTasksAsBefore(testCase)
            % the user's MATLAB (no Slurm job) has all of the machine: it
            % runs the last task whatever the estimate, as before
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.setSlurmJob('');
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 1));
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [logFn, ~, outFns] = testCase.runFramework(3, 'memAllocate', 400);
            testCase.verifyEqual(testCase.readLog(logFn), outFns(3));
            testCase.verifyEqual(slurm.callCount('sbatch'), 2);
        end

        function localRunsRunAllTasksAsBefore(testCase)
            % without a cluster the master runs every task, also in a job
            testCase.setSlurmJob('4242');
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 1));
            [logFn, ~, outFns] = testCase.runFramework(3, 'memAllocate', 400, 'parseCluster', false);
            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort(outFns));
        end

        function fitIsTheTasksMemoryAgainstTheJobs(testCase)
            testCase.setSlurmJob('4242');
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 16));
            % memAllocate, up to 5% over (memPerCPU is often rounded)
            testCase.verifyTrue(master_task_fits(16.5, 1, 20.9, 24, 1));
            testCase.verifyFalse(master_task_fits(17, 1, 20.9, 24, 1));
            % at most what a job gets (maxCPUNum x memPerCPU)
            testCase.verifyTrue(master_task_fits(100, 1, 10, 1, 1));
            % without memAllocate, the memory of the job's CPUs (at least
            % minCPUNum, at most maxCPUNum)
            testCase.verifyFalse(master_task_fits([], 2, 10, 24, 1));
            testCase.verifyTrue(master_task_fits([], 2, 10, 1, 1));
            testCase.verifyFalse(master_task_fits([], 1, 10, 24, 2));
            [~, taskMem, masterMem] = master_task_fits([], 3, 5, 24, 1);
            testCase.verifyEqual([taskMem, masterMem], [15, 16]);
        end

        function unknownMemoryDoesNotFit(testCase)
            % (getSystemMemory gives NaN when it can't find the memory)
            testCase.setSlurmJob('4242');
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', NaN));
            testCase.verifyFalse(master_task_fits(1, 1, 20.9, 24, 1));
        end

        function everyTaskFitsOutsideAJob(testCase)
            testCase.setSlurmJob('');
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 1));
            [tf, ~, masterMem] = master_task_fits(1000, 24, 20.9, 24, 1);
            testCase.verifyTrue(tf);
            testCase.verifyEqual(masterMem, inf);
        end

        function coordinatorAsksForWhatItHolds(testCase)
            % MATLAB (2 GB) and the merges with one CPU; a merge needs up to
            % half of the job's memory (XR_MIP_zarr)
            [cpus, mem] = coordinator_resources(0.5, false, 24, 40);
            testCase.verifyEqual([cpus, mem], [1, 2.5]);
            [cpus, mem] = coordinator_resources(3, false, 24, 40);
            testCase.verifyEqual([cpus, mem], [1, 6]);
            % one of its tasks too, if it computes them
            [cpus, mem] = coordinator_resources(0.5, true, 24, 40);
            testCase.verifyEqual([cpus, mem], [24, 40]);
            [cpus, mem] = coordinator_resources(0.5, true, 24, []);
            testCase.verifyEqual([cpus, mem], [24, 2.5]);
            [cpus, mem] = coordinator_resources(30, true, 2, 40);
            testCase.verifyEqual([cpus, mem], [2, 60]);
        end

        function configValuesIgnoreCase(testCase)
            folder = testCase.createTempFolder();
            fn = [folder, 'config.json'];
            fid = fopen(fn, 'w');
            fprintf(fid, '{"MasterCompute": false, "memPerCPU": 7.6}');
            fclose(fid);
            testCase.verifyEqual(cluster_config_value(fn, 'masterCompute', true), false);
            testCase.verifyEqual(cluster_config_value(fn, 'MemPerCPU', 20.9), 7.6);
            testCase.verifyEqual(cluster_config_value(fn, 'maxCPUNum', 24), 24);
            testCase.verifyEqual(cluster_config_value('', 'masterCompute', true), true);
            testCase.verifyEqual(cluster_config_value([folder, 'none.json'], 'masterCompute', true), true);
            matFn = [folder, 'config.mat'];
            masterCompute = false;
            save(matFn, 'masterCompute');
            testCase.verifyEqual(cluster_config_value(matFn, 'masterCompute', true), false);
        end

        function mergeMemoryCoversTheMergeOfTheBlocks(testCase)
            % the estimate is at least what MIP_zarr_merge_block_mips holds
            % for the blocks (computed here from their boxes as it does),
            % and not much more: decon-like grids (also not aligned with
            % the 2048 MIP batches) and DSR-like y batches across x and z
            cases = {[4096, 4096, 2048], [1024, 1024, 1024], [256, 256, 256], 'uint16'
                [5000, 3000, 2100], [768, 768, 768], [256, 256, 256], 'uint16'
                [2200, 300, 64], [768, 300, 64], [128, 128, 64], 'single'
                [8192, 2673, 1095], [1024, 2673, 1095], [256, 256, 256], 'uint16'
                [100, 120, 1], [64, 64, 1], [16, 16, 1], 'uint16'};
            for i = 1 : size(cases, 1)
                [imSize, regionSize, blockSize, dtype] = cases{i, :};
                est = MIP_zarr_merge_memory(imSize, regionSize, blockSize, dtype);
                exact = testCase.exactMergeMemory(imSize, regionSize, MIP_zarr_batch_size(imSize, blockSize, dtype), dtype);
                testCase.verifyGreaterThanOrEqual(est, exact, sprintf('Case %d.', i));
                testCase.verifyLessThanOrEqual(est, exact * 1.5, sprintf('Case %d.', i));
            end
        end
    end

    methods (Access = private)
        function [] = setSlurmJob(testCase, jobId)
            % SLURM_JOB_ID for the test ('' unsets it)
            orig = getenv('SLURM_JOB_ID');
            if isempty(orig)
                testCase.addTeardown(@unsetenv, 'SLURM_JOB_ID');
            else
                testCase.addTeardown(@setenv, 'SLURM_JOB_ID', orig);
            end
            if isempty(jobId)
                unsetenv('SLURM_JOB_ID');
            else
                setenv('SLURM_JOB_ID', jobId);
            end
        end

        function [logFn, out, outFns] = runFramework(testCase, n, varargin)
            % the framework in cluster mode with masterCompute on n copy
            % tasks, one try; varargin adds or overrides parameters
            dataPath = testCase.createTempFolder();
            mkdir([dataPath, 'out']);
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : n, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout/out_%d.txt', dataPath, x), 1 : n, 'unif', 0);
            for f = 1 : n
                fid = fopen(inFns{f}, 'w');
                fprintf(fid, 'task %d', f);
                fclose(fid);
            end
            logFn = [dataPath, 'calls.log'];
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''copy'')', x, y, logFn), ...
                inFns, outFns, 'unif', 0);
            args = [{'parseCluster', true, 'masterCompute', true, 'maxTrialNum', 1, 'queryInterval', 0, ...
                'jobLogDir', [dataPath, 'logs'], 'taskRunner', false}, varargin];
            out = evalc('generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, args{:})');
        end

        function [lines] = readLog(~, logFn)
            % output paths of all task calls the master made, in order
            lines = cell(1, 0);
            if exist(logFn, 'file')
                lines = splitlines(strip(fileread(logFn)))';
            end
        end

        function [mem] = exactMergeMemory(~, imSize, regionSize, mipBatchSize, dtype)
            % MIP_zarr_merge_block_mips' memory for the blocks of a grid of
            % regionSize, the largest axis
            nR = ceil(imSize ./ regionSize);
            [y, x, z] = ndgrid(1 : nR(1), 1 : nR(2), 1 : nR(3));
            starts = ([y(:), x(:), z(:)] - 1) .* regionSize + 1;
            boxes = [starts, min(starts + regionSize - 1, imSize)];
            byteNum = dataTypeToByteNumber(dtype);
            mem = 0;
            for a = 1 : 3
                slabSize = imSize;
                slabSize(a) = ceil(imSize(a) / mipBatchSize(a));
                k1 = floor((boxes(:, a) - 1) / mipBatchSize(a)) + 1;
                k2 = floor((boxes(:, a + 3) - 1) / mipBatchSize(a)) + 1;
                mipSizes = boxes(:, 4 : 6) - boxes(:, 1 : 3) + 1;
                mipSizes(:, a) = k2 - k1 + 1;
                mem = max(mem, (prod(slabSize) * (byteNum + 1) + sum(prod(mipSizes, 2)) * byteNum) / 2^30);
            end
        end
    end
end
