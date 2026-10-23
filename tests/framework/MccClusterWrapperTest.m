classdef MccClusterWrapperTest < pktest.TestCase
    % mcc_slurm_cluster_generic_computing_wrapper, the job framework's
    % compiled mode, with a fake Slurm (pktest.FakeSlurmFixture) and a fake
    % compiled master that records its environment. The jobs never run
    % (a fake job can only leave empty outputs); the tests run the task
    % lines the framework writes for them.

    methods (TestMethodSetup)
        function clearConfigCache(testCase)
            % config file contents are cached in persistent variables keyed
            % on the file name and its modification time
            clear generic_computing_frameworks_wrapper
        end
    end

    methods (Test, TestTags = {'Unit'})
        function tasksUseOwnRuntimeCachePerSlot(testCase)
            % tasks of a job that started at the same time with a shared
            % MATLAB Runtime cache sometimes failed to start
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.applyFixture(pktest.FakeSlurmFixture);
            [lines, logFn, cacheRoot] = testCase.writeTaskLines();
            testCase.assertNotEmpty(lines, 'No task lines were written.');

            testCase.runTaskLine(lines{1}, 'SLURM_JOB_ID=42 PARALLEL_JOBSLOT=3');
            testCase.runTaskLine(lines{1}, 'SLURM_JOB_ID=42 PARALLEL_JOBSLOT=4');
            caches = splitlines(strip(fileread(logFn)));
            testCase.verifyEqual(caches, {sprintf('%s/mcr_cache_42_3', cacheRoot); sprintf('%s/mcr_cache_42_4', cacheRoot)});
        end

        function submitsEachTaskOnceWithMaxJobNum(testCase)
            % with more jobs than maxJobNum the framework holds tasks back
            % until earlier jobs end (it counts its active jobs for every
            % batch); the fake jobs end at once without results
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.runFramework(3, 'maxJobNum', 1);
            arrays = cellfun(@(x) str2double(regexp(x, '--array=(\d+)', 'tokens', 'once')), slurm.calls('sbatch'));
            testCase.verifyEqual(sort(arrays(:))', 1 : 3, 'Each task should be submitted once.');
        end

        function taskFilesThatCannotBeRemovedDoNotFailTheRun(testCase)
            % on NFS, a job that still has its task script open keeps a .nfs
            % file in the task folder, so removing the folder at the end
            % failed and the run errored although its outputs were complete
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            slurm.setJobStates('1000_1', {'R', 'done'}, {[dataPath, 'out_1.txt']});
            rmdirs = testCase.applyFixture(pktest.FunctionSpyFixture('rmdir', ...
                @() error('pktest:rmdir', 'Directory not empty (a .nfs file)')));

            testCase.runFrameworkIn(dataPath, 1, 'queryInterval', 0);

            testCase.verifyNotEmpty(rmdirs.getCalls(), 'The task files were not removed.');
            testCase.verifyTrue(isfile([dataPath, 'out_1.txt']));
        end

        function countsJobUntilAllItsTasksAreDone(testCase)
            % with taskBatchNum > 1 a job runs several tasks; it counts
            % against maxJobNum until all of them are done, not just one
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            % the job of batch 1 (tasks 1 and 2) writes the output of task
            % 1 at the first poll and is gone at the fourth
            slurm.setJobStates('1000_1', {'R', 'R', 'R', 'done'}, {[dataPath, 'out_1.txt']});

            testCase.runFrameworkIn(dataPath, 4, 'taskBatchNum', 2, 'maxJobNum', 1, 'queryInterval', 0);

            cmds = slurm.history();
            submit2 = find(startsWith(cmds, 'sbatch') & contains(cmds, '--array=2 '));
            testCase.assertNumElements(submit2, 1, 'Batch 2 should be submitted once.');
            polls = find(startsWith(cmds, 'squeue') & contains(cmds, '1000_1'));
            testCase.assertGreaterThanOrEqual(numel(polls), 4, 'The job of batch 1 was not polled until it ended.');
            testCase.verifyGreaterThan(submit2, polls(4), ...
                'Batch 2 was submitted while the job of batch 1 was still running.');
        end

        function finishedShortLastBatchFreesItsJobSlot(testCase)
            % a batch whose tasks are all done stops counting against
            % maxJobNum even while its job is still listed (e.g., exiting).
            % Here the batch is the last one and shorter than taskBatchNum;
            % the count pads it to taskBatchNum tasks, which must count as done.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            % Slurm rejects batch 1 (tasks 1 and 2), so batch 2 (task 3) is
            % submitted first; its job writes the output at the first poll
            % and stays listed
            slurm.rejectSubmissions(1);
            slurm.setJobStates('1000_2', {'R'}, {[dataPath, 'out_3.txt']});
            % a framework that holds batch 1 back forever is stopped by a
            % watchdog that writes the final output
            finalFn = [dataPath, 'final.txt'];
            watchdog = timer('StartDelay', 20, 'TimerFcn', @(~, ~) fclose(fopen(finalFn, 'w')));
            testCase.addTeardown(@() delete(watchdog));
            testCase.addTeardown(@() stop(watchdog));
            start(watchdog);

            testCase.runFrameworkIn(dataPath, 3, 'taskBatchNum', 2, 'maxJobNum', 1, 'queryInterval', 0, ...
                'finalOutFullpath', finalFn);

            testCase.verifyFalse(exist(finalFn, 'file') == 2, ...
                'The finished batch 2 kept counting as a job until the watchdog stopped the framework.');
            arrays = cellfun(@(x) str2double(regexp(x, '--array=(\d+)', 'tokens', 'once')), slurm.calls('sbatch'));
            testCase.verifyEqual(arrays(:)', [1, 2, 1], 'Batch 1 should be submitted again once batch 2 is done.');
        end

        function tasksOutsideGNUParallelUseOwnRuntimeCache(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.applyFixture(pktest.FakeSlurmFixture);
            [lines, logFn, cacheRoot] = testCase.writeTaskLines();
            testCase.assertNotEmpty(lines, 'No task lines were written.');

            testCase.runTaskLine(lines{1}, '');
            cache = strip(fileread(logFn));
            testCase.verifyMatches(cache, ['^', regexptranslate('escape', cacheRoot), '/mcr_cache_0_\d+$']);
        end
    end

    methods (Access = private)
        function [lines, logFn, cacheRoot] = writeTaskLines(testCase)
            % runs the framework in mcc mode on one task; the fake Slurm never
            % runs the job, so the framework gives up after one try and
            % leaves the task file behind
            [dataPath, logFn, cacheRoot] = testCase.runFramework(1);
            taskFiles = dir([dataPath, 'tmp/*/input_*.txt']);
            lines = {};
            if ~isempty(taskFiles)
                lines = splitlines(strip(fileread(fullfile(taskFiles(1).folder, taskFiles(1).name))));
            end
        end

        function [dataPath, logFn, cacheRoot] = runFramework(testCase, nTasks, varargin)
            % runFrameworkIn in a new folder
            dataPath = testCase.createTempFolder();
            [logFn, cacheRoot] = testCase.runFrameworkIn(dataPath, nTasks, varargin{:});
        end

        function [logFn, cacheRoot] = runFrameworkIn(~, dataPath, nTasks, varargin)
            % the framework in mcc mode on nTasks tasks in dataPath (inputs
            % in_<k>.txt, outputs out_<k>.txt), with a fake MATLAB Runtime
            % and a fake compiled master that logs MCR_CACHE_ROOT
            mcrPath = [dataPath, 'mcr/'];
            mkdir([mcrPath, 'bin']);
            logFn = [dataPath, 'caches.log'];
            masterFn = [dataPath, 'run_mccMaster.sh'];
            fid = fopen(masterFn, 'w');
            fprintf(fid, '#!/bin/bash\necho "$MCR_CACHE_ROOT" >> "%s"\n', logFn);
            fclose(fid);
            fileattrib(masterFn, '+x');
            cacheRoot = [dataPath, 'cache'];
            configFn = [dataPath, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(struct('MCRParam', mcrPath, 'MCCMasterStr', masterFn, ...
                'MCRCacheRoot', cacheRoot, 'BashLaunchStr', 'echo')));
            fclose(fid);
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : nTasks, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout_%d.txt', dataPath, x), 1 : nTasks, 'unif', 0);
            cellfun(@(x) fclose(fopen(x, 'w')), inFns);

            % task strings as the data wrappers write them (positional
            % arguments, then name-value pairs); the fake master ignores them
            funcStrs = cellfun(@(x, y) sprintf('XR_crop_frame(''%s'',''%s'',[1,1,1,2,2,2],''pad'',false)', x, y), ...
                inFns, outFns, 'unif', 0);
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''mccMode'', true, ''maxTrialNum'', 1, ''jobLogDir'', dataPath, ', ...
                '''ConfigFile'', configFn, varargin{:})']);
        end

        function runTaskLine(testCase, line, envStr)
            % run a task line as a job would, with bash and the given environment
            lineFn = [tempname, '.sh'];
            fid = fopen(lineFn, 'w');
            fprintf(fid, '%s\n', line);
            fclose(fid);
            testCase.addTeardown(@() delete(lineFn));
            [status, out] = system(sprintf('env -u SLURM_JOB_ID -u PARALLEL_JOBSLOT %s bash %s', envStr, lineFn));
            testCase.assertEqual(status, 0, out);
        end
    end
end
