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

        function outdatedBuildStopsBeforeSubmitting(testCase)
            % a compiled master built from code with other task interfaces
            % rejects the tasks' options; the run stops at once with a clear
            % error instead of every task failing (here a build from before
            % XR_crop_block, the tasks' function, got mipBatchSize)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            manifest = strrep(mcc_parser_manifest(fileparts(which('mccMaster'))), ...
                'Overwrite,mipBatchSize,uuid', 'Overwrite,uuid');
            fid = fopen([dataPath, 'parsers_manifest.txt'], 'w');
            fprintf(fid, '%s', manifest);
            fclose(fid);

            testCase.verifyError(@() testCase.runFrameworkIn(dataPath, 2), 'mcc:outdatedBuild');
            testCase.verifyEqual(slurm.callCount('sbatch'), 0, 'Tasks were submitted to the outdated build.');
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

        function failedQueryKeepsJobStates(testCase)
            % a failed squeue (e.g., a timeout of a busy controller) read as
            % all jobs gone: they stopped counting against maxJobNum (and
            % were submitted again after 45 s); they keep their states now
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            % the first two polls fail; the job of batch 1 then runs and
            % writes its output
            slurm.failQueries(2);
            slurm.setJobStates('1000_1', {'R', 'R', 'done'}, {[dataPath, 'out_1.txt']});

            testCase.runFrameworkIn(dataPath, 2, 'maxJobNum', 1, 'queryInterval', 0);

            cmds = slurm.history();
            submit2 = find(startsWith(cmds, 'sbatch') & contains(cmds, '--array=2 '));
            testCase.assertNumElements(submit2, 1, 'Batch 2 should be submitted once.');
            polls = find(startsWith(cmds, 'squeue') & contains(cmds, '1000_1'));
            testCase.assertGreaterThanOrEqual(numel(polls), 3, 'The job of batch 1 was not polled until it ran.');
            testCase.verifyGreaterThan(submit2, polls(3), ...
                'Batch 2 was submitted while squeue failed to answer about the job of batch 1.');
        end

        function lostSubmissionAnswerFindsTheJob(testCase)
            % sbatch failed with "Socket timed out" after Slurm made the job,
            % and the framework submitted the batch again (also past
            % maxJobNum): two jobs, one unknown to it. It now looks for the
            % job by its comment and uses it, and until then the batch
            % counts against maxJobNum
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            slurm.loseSubmissionAnswers(1);
            slurm.setJobStates('1000_1', {'R'}, {[dataPath, 'out_1.txt']});

            % (rounds without a query between the queries)
            testCase.runFrameworkIn(dataPath, 2, 'maxJobNum', 1, 'queryInterval', 1);

            cmds = slurm.history();
            arrays = cellfun(@(x) str2double(regexp(x, '--array=(\d+)', 'tokens', 'once')), slurm.calls('sbatch'));
            testCase.verifyEqual(arrays(:)', [1, 2], strjoin(cmds, newline));
            lookup = find(startsWith(cmds, 'squeue -u'), 1);
            testCase.assertNotEmpty(lookup, 'The job was not looked for.');
            submit2 = find(startsWith(cmds, 'sbatch') & contains(cmds, '--array=2 '));
            testCase.verifyGreaterThan(submit2, lookup, 'Batch 2 was submitted past maxJobNum.');
        end

        function lostSubmissionWaitsWhileSqueueFails(testCase)
            % a batch submitted without an answer isn't submitted again
            % while its job can't be looked for (squeue fails, here the
            % first two times)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            slurm.loseSubmissionAnswers(1);
            slurm.failQueries(2);
            slurm.setJobStates('1000_1', {'R'}, {[dataPath, 'out_1.txt']});

            testCase.runFrameworkIn(dataPath, 1, 'queryInterval', 0);

            testCase.verifyNumElements(slurm.calls('sbatch'), 1, strjoin(slurm.history(), newline));
            testCase.verifyTrue(isfile([dataPath, 'out_1.txt']), 'The job found did not run the task.');
        end

        function runExtraTasksResubmitsBatchOne(testCase)
            % with runExtraTasks and master compute, the last waiting batch
            % was left out of the submissions also when the master doesn't
            % take it (batch 1, whose job runs other batches' tasks too): a
            % batch 1 whose job failed was never submitted again, and the
            % call didn't end (here a watchdog ends it after 60 s)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            % the first job of batch 1 ends without its output, the second
            % one writes it; the master runs batch 2 (its fake compiled
            % master writes the output, its fourth argument)
            slurm.setJobStates('1000_1', {'done'});
            slurm.setJobStates('1001_1', {'R'}, {[dataPath, 'out_1.txt']});
            finalFn = [dataPath, 'final.txt'];
            watchdog = timer('StartDelay', 60, 'TimerFcn', @(~, ~) fclose(fopen(finalFn, 'w')));
            testCase.addTeardown(@() delete(watchdog));
            testCase.addTeardown(@() stop(watchdog));
            start(watchdog);

            flags = testCase.runCompiledMaster(dataPath, 2, {'#!/bin/bash', 'echo x > "$4"'}, ...
                'runExtraTasks', true, 'finalOutFullpath', finalFn, 'queryInterval', 1);

            testCase.verifyFalse(isfile(finalFn), 'Batch 1 was not submitted again (the watchdog ended the call).');
            testCase.verifyEqual(flags, true(2, 1));
            arrays = cellfun(@(x) str2double(regexp(x, '--array=(\d+)', 'tokens', 'once')), slurm.calls('sbatch'));
            testCase.verifyEqual(arrays(:)', [1, 1]);
        end

        function stoppedCallCancelsItsJobs(testCase)
            % a call stopped while it waits for its jobs (Ctrl+C; here an
            % error of pause) or ended by an error left them in the queue;
            % here also the job of a submission whose answer was lost
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            slurm.setJobStates('1000_1', {'PD'});
            slurm.setJobStates('1001_2', {'PD'});
            slurm.loseSubmissionAnswers(1);
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() error('pktest:stopped', 'Stopped.')));

            testCase.verifyError(@() testCase.runFrameworkIn(dataPath, 2, 'maxTrialNum', 2), 'pktest:stopped');

            testCase.verifyEqual(slurm.calls('scancel'), {'1000 1001'}, 'The jobs were left in the queue.');
        end

        function pausesWhileOnlyLastTriesRun(testCase)
            % while only jobs of last tries were listed, the framework looped
            % without a pause (a busy core) until they ended
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            slurm.setJobStates('1000_1', {'R', 'R', 'R', 'done'});
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() builtin('pause', 0.05)));

            testCase.runFrameworkIn(dataPath, 1, 'queryInterval', 0.2);

            testCase.verifyGreaterThanOrEqual(numel(pauses.getCalls()), 3, ...
                'The framework did not pause while the job of the last try ran.');
        end

        function gpuSlotsUseTheJobsGPUs(testCase)
            % GNU parallel gave slot k GPU k - 1, another job's GPU where
            % Slurm doesn't renumber the GPUs of a job (it doesn't where it
            % doesn't constrain them); slot k takes the k-th GPU in the job's
            % CUDA_VISIBLE_DEVICES now, a job with fewer GPUs than slots runs
            % fewer slots, and without a list slot k takes GPU k - 1 as
            % before. The job script runs here with GPU lists of a test (no
            % GPU is used)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.assumeEqual(system('command -v parallel > /dev/null'), 0, 'GNU parallel is not installed.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            gpuFn = [dataPath, 'gpus.log'];
            % each task logs its GPU and slot (the tasks overlap: each slot
            % runs one)
            testCase.runFrameworkIn(dataPath, 2, 'GPUJob', true, 'GNUparallel', true, 'paraJobNum', 2, ...
                'taskBatchNum', 2, 'masterLines', {'#!/bin/bash', ...
                sprintf('echo "${CUDA_VISIBLE_DEVICES-unset} ${PARALLEL_JOBSLOT-}" >> "%s"', gpuFn), 'sleep 1'});
            sbatch = slurm.calls('sbatch');
            testCase.assertNumElements(sbatch, 1);
            jobFn = [dataPath, 'job.sh'];
            fid = fopen(jobFn, 'w');
            fprintf(fid, '%s\n', extractAfter(sbatch{1}, '--wrap='));
            fclose(fid);

            % GPU lists, and the GPUs and slots of the tasks
            lists = {'CUDA_VISIBLE_DEVICES=4,6', {'4'; '6'}, 2; 'CUDA_VISIBLE_DEVICES=5', {'5'; '5'}, 1; ...
                '-u CUDA_VISIBLE_DEVICES', {'0'; '1'}, 2};
            for i = 1 : size(lists, 1)
                if isfile(gpuFn)
                    delete(gpuFn);
                end
                [status, out] = system(sprintf('env -u SLURM_JOB_ID -u PARALLEL_JOBSLOT %s bash %s 2>&1', lists{i, 1}, jobFn));
                testCase.assertEqual(status, 0, out);
                runs = split(splitlines(strip(fileread(gpuFn))));
                msg = [lists{i, 1}, newline, out];
                testCase.verifyEqual(sort(runs(:, 1)), lists{i, 2}, msg);
                testCase.verifyNumElements(unique(runs(:, 2)), lists{i, 3}, msg);
            end
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

        function [logFn, cacheRoot] = runFrameworkIn(testCase, dataPath, nTasks, varargin)
            % the framework in mcc mode on nTasks tasks in dataPath (inputs
            % in_<k>.txt, outputs out_<k>.txt), with a fake MATLAB Runtime
            % and a fake compiled master that logs MCR_CACHE_ROOT (or has
            % the lines of the option 'masterLines')
            logFn = [dataPath, 'caches.log'];
            masterLines = {'#!/bin/bash', sprintf('echo "$MCR_CACHE_ROOT" >> "%s"', logFn)};
            k = find(strcmp(varargin(1 : 2 : end), 'masterLines'), 1);
            if ~isempty(k)
                masterLines = varargin{2 * k};
                varargin(2 * k - 1 : 2 * k) = [];
            end
            [inFns, outFns, funcStrs, s] = testCase.compiledTasks(dataPath, nTasks, masterLines); %#ok<ASGLU>
            cacheRoot = s.MCRCacheRoot;
            configFn = [dataPath, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(s));
            fclose(fid);
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''mccMode'', true, ''maxTrialNum'', 1, ''jobLogDir'', dataPath, ', ...
                '''ConfigFile'', configFn, varargin{:})']);
        end

        function [flags] = runCompiledMaster(testCase, dataPath, nTasks, masterLines, varargin)
            % mcc_slurm_cluster_generic_computing_wrapper with master compute
            % (as in a compiled master) on the tasks of runFrameworkIn
            [inFns, outFns, funcStrs, s] = testCase.compiledTasks(dataPath, nTasks, masterLines); %#ok<ASGLU>
            args = [{'MCRParam', s.MCRParam, 'MCCMasterStr', s.MCCMasterStr, 'MCRCacheRoot', s.MCRCacheRoot, ...
                'BashLaunchStr', s.BashLaunchStr}, varargin]; %#ok<NASGU>
            evalc(['flags = mcc_slurm_cluster_generic_computing_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', true, ''maxTrialNum'', 3, ''jobLogDir'', dataPath, args{:});']);
        end

        function [inFns, outFns, funcStrs, s] = compiledTasks(~, dataPath, nTasks, masterLines)
            % nTasks tasks in dataPath (inputs in_<k>.txt, outputs
            % out_<k>.txt), a fake MATLAB Runtime and a fake compiled master
            % (masterLines), and the compiled mode's settings s
            mcrPath = [dataPath, 'mcr/'];
            mkdir([mcrPath, 'bin']);
            masterFn = [dataPath, 'run_mccMaster.sh'];
            fid = fopen(masterFn, 'w');
            fprintf(fid, '%s\n', masterLines{:});
            fclose(fid);
            fileattrib(masterFn, '+x');
            s = struct('MCRParam', mcrPath, 'MCCMasterStr', masterFn, 'MCRCacheRoot', [dataPath, 'cache'], ...
                'BashLaunchStr', 'echo');
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : nTasks, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout_%d.txt', dataPath, x), 1 : nTasks, 'unif', 0);
            cellfun(@(x) fclose(fopen(x, 'w')), inFns);

            % task strings as the data wrappers write them (positional
            % arguments, then name-value pairs); the fake master ignores them
            funcStrs = cellfun(@(x, y) sprintf('XR_crop_frame(''%s'',''%s'',[1,1,1,2,2,2],''pad'',false)', x, y), ...
                inFns, outFns, 'unif', 0);
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
