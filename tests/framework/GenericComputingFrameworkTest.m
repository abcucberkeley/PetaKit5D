classdef GenericComputingFrameworkTest < pktest.TestCase
    % generic_computing_frameworks_wrapper in local mode (parseCluster =
    % false): tasks run sequentially in this MATLAB session. Tasks call
    % pktest.frameworkTask, which logs every call so runs can be counted.
    % The clusterMode tests submit to a fake Slurm (pktest.FakeSlurmFixture).

    methods (TestMethodSetup)
        function clearConfigCache(testCase)
            % config file contents are cached in persistent variables keyed
            % on the file name and its modification time
            clear generic_computing_frameworks_wrapper
        end
    end

    methods (Test, TestTags = {'Unit'})
        function runsEveryTaskOnceAndWritesOutputs(testCase)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(3);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true);

            testCase.verifyEqual(flags, true(3, 1));
            for f = 1 : 3
                testCase.verifyEqual(fileread(outFns{f}), fileread(inFns{f}));
            end
            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort(outFns(:)), 'Each task should run exactly once.');
            testCase.verifyEmpty(dir([dataPath, 'out/*.tmp']), 'Task flag files should be removed after success.');
        end

        function skipsTasksWhoseOutputsExist(testCase)
            [inFns, outFns, logFn] = testCase.makeTasks(3);
            cellfun(@copyfile, inFns, outFns);
            % the tasks would throw if they ran
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'error');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true);

            testCase.verifyEqual(flags, true(3, 1));
            testCase.verifyEmpty(testCase.readLog(logFn), 'No task should run when all outputs exist.');
        end

        function runsOnlyTasksWithMissingOutputs(testCase)
            [inFns, outFns, logFn] = testCase.makeTasks(3);
            copyfile(inFns{1}, outFns{1});
            copyfile(inFns{2}, outFns{2});
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true);

            testCase.verifyEqual(flags, true(3, 1));
            testCase.verifyEqual(testCase.readLog(logFn), outFns(3));
        end

        function failingTaskIsRetriedOnceThenRethrown(testCase)
            [inFns, outFns, logFn] = testCase.makeTasks(1);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'error');

            testCase.verifyError(@() generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true), 'pktest:frameworkTask:failed');
            % one run in try, one more in catch, then the error propagates
            testCase.verifyNumElements(testCase.readLog(logFn), 2);
        end

        function taskWithoutOutputIsRetriedUpToMaxTrialNum(testCase)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop');

            % unitWaitTime = 0, otherwise each retry busy-waits for minutes
            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true, 'maxTrialNum', 2, 'unitWaitTime', 0);

            testCase.verifyEqual(flags, false(2, 1));
            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort([outFns(:); outFns(:)]), ...
                'Each task should run maxTrialNum times.');
            % the flag files of failed tasks stay behind; a rerun within
            % unitWaitTime minutes skips (busy-waits on) those tasks
            testCase.verifyNumElements(dir([dataPath, 'out/*.tmp']), 2);
        end

        function taskDoneInItsLastTryCountsAsDone(testCase)
            % a task that makes its output in its last try ran after its
            % round looked at the outputs; the framework looked no more, so
            % the last task run returned unfinished (and callers ran the
            % call again); it now looks once more at the end
            [inFns, outFns, logFn] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'secondRun');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true, 'maxTrialNum', 2, 'unitWaitTime', 0);

            testCase.verifyEqual(flags, true(2, 1), 'A task done in its last try counts as unfinished.');
            testCase.verifyTrue(all(cellfun(@isfile, outFns)));
            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort([outFns(:); outFns(:)]));
        end

        function staleTaskFlagFromCrashedRunIsIgnored(testCase)
            testCase.assumeTrue(usejava('jvm'), 'Needs Java to set the file modification time.');
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');
            % flag file left by a worker that died an hour ago
            tmpFn = [dataPath, 'out/out_2.tmp'];
            fclose(fopen(tmpFn, 'w'));
            java.io.File(tmpFn).setLastModified(java.lang.System.currentTimeMillis() - 3600 * 1000);

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true);

            testCase.verifyEqual(flags, true(2, 1));
            testCase.verifyTrue(exist(outFns{2}, 'file') == 2);
            testCase.verifyFalse(exist(tmpFn, 'file') == 2, 'The stale flag file should be removed.');
        end

        function missingInputIsSkippedAndReportedDone(testCase)
            % characterizes current behavior: tasks whose input is missing
            % are not run, and their done flag is true even though no
            % output exists
            [inFns, outFns, logFn] = testCase.makeTasks(2);
            delete(inFns{2});
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true);

            testCase.verifyEqual(flags, true(2, 1));
            testCase.verifyEqual(testCase.readLog(logFn), outFns(1));
            testCase.verifyFalse(exist(outFns{2}, 'file') == 2);
        end

        function existingFinalOutputSkipsAllTasks(testCase)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            finalFn = [dataPath, 'final.txt'];
            fclose(fopen(finalFn, 'w'));
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'error');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true, 'finalOutFullpath', finalFn);

            testCase.verifyEqual(flags, true(2, 1));
            testCase.verifyEmpty(testCase.readLog(logFn));
        end

        function finalOutputMadeByATaskEndsTheRun(testCase)
            % the first task writes the final output; the second never
            % writes its own and would be retried up to maxTrialNum times,
            % but the run ends once the final output exists (in local mode
            % after the first pass over the tasks)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            outFns{1} = [dataPath, 'final.txt'];
            funcStrs = [testCase.taskStrings(inFns(1), outFns(1), logFn, 'copy'), ...
                testCase.taskStrings(inFns(2), outFns(2), logFn, 'noop')];

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, 'parseCluster', false, ...
                'masterCompute', true, 'unitWaitTime', 0, 'finalOutFullpath', outFns{1});

            testCase.verifyEqual(flags, true(2, 1));
            testCase.verifyEqual(testCase.readLog(logFn), outFns(:), 'Each task should run once.');
        end

        function taskBatchesRunEachTaskOnce(testCase)
            [inFns, outFns, logFn] = testCase.makeTasks(4);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true, 'taskBatchNum', 2);

            testCase.verifyEqual(flags, true(4, 1));
            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort(outFns(:)));
        end

        function configFileOverridesParameters(testCase)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop');
            configFn = [dataPath, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(struct('maxTrialNum', 1)));
            fclose(fid);

            % maxTrialNum = 3 here, the config file sets it to 1
            generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, 'parseCluster', false, ...
                'masterCompute', true, 'maxTrialNum', 3, 'unitWaitTime', 0, 'ConfigFile', configFn);

            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort(outFns(:)), ...
                'With maxTrialNum = 1 from the config file, each task should run once.');
        end

        function configFileKeysIgnoreCase(testCase)
            % the cluster configs in use spell MemPerCPU; config keys used to
            % be matched case-sensitively, so it was silently ignored
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(1);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop');
            configFn = [dataPath, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(struct('MemPerCPU', 5)));
            fclose(fid);

            % 50 GB at 5 GB per CPU is 10 CPUs, rounded up to 12 to fit
            % 24-CPU nodes (the default 20.9 GB per CPU would give 3)
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''memAllocate'', 50, ''jobLogDir'', dataPath, ', ...
                '''ConfigFile'', configFn)']);

            sbatch = slurm.calls('sbatch');
            testCase.assertNotEmpty(sbatch, 'No job was submitted.');
            testCase.verifySubstring(sbatch{1}, '--cpus-per-task=12');
        end

        function emptyJobLogDirInAConfigFileIsTheDefaultFolder(testCase)
            % a cluster config with "jobLogDir": "" (the user's data-processing
            % configs have it) gave checkSlurmCluster no folder, so the jobs
            % got "-o  -e" and each one's log went to a file named -e in the
            % folder the framework submits from (the code folder); it is now
            % the default folder (../job_logs beside the code folder if it
            % is there, else <data folder>/job_logs). The fake jobs make
            % their log files as Slurm does when they start (only in the
            % test's folders)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.makeLogFiles();
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop'); %#ok<NASGU>
            configFn = [dataPath, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(struct('jobLogDir', '')));
            fclose(fid);
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(dataPath));

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''taskRunner'', false, ''ConfigFile'', configFn)']);

            [logDir, inData] = pktest.defaultJobLogDir(dataPath);
            sbatch = slurm.calls('sbatch');
            testCase.assertNotEmpty(sbatch, 'No job was submitted.');
            logs = regexp(sbatch, ' -o (\S+) -e (\S+) ', 'tokens', 'once');
            testCase.verifyTrue(all(cellfun(@(x) numel(x) == 2 && strcmp(pktest.fullPath(fileparts(x{1})), logDir) ...
                && endsWith(x{1}, '/job_%A_%a.out') && strcmp(x{2}, strrep(x{1}, '.out', '.err')), logs)), ...
                sprintf('Not in %s:\n%s', logDir, strjoin(sbatch, newline)));
            if inData
                testCase.verifyNumElements(dir([dataPath, 'job_logs/job_*.out']), numel(sbatch), ...
                    'The job logs are not in the default folder.');
            end
        end

        function clusterModeSubmitsEachTaskOnceWithMaxJobNum(testCase)
            % with more jobs than maxJobNum the framework holds tasks back
            % until earlier jobs end (it counts its active jobs for every
            % batch); the fake jobs end at once without results
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(3);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop');

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''maxJobNum'', 1, ''jobLogDir'', dataPath, ', ...
                '''taskRunner'', false)']);

            arrays = cellfun(@(x) str2double(regexp(x, '--array=(\d+)', 'tokens', 'once')), slurm.calls('sbatch'));
            testCase.verifyEqual(sort(arrays(:))', 1 : 3, 'Each task should be submitted once.');
            % the jobs start MATLAB in its non-interactive mode (with -r it
            % starts the MathWorks Service Host and was very slow to exit;
            % TaskRunnerTest checks the task runners' job script)
            testCase.verifyTrue(all(contains(slurm.calls('sbatch'), '-nojvm -batch ')), ...
                'The jobs do not start MATLAB with -batch.');
        end

        function clusterModeCountsJobUntilAllItsTasksAreDone(testCase)
            % with taskBatchNum > 1 a job runs several tasks; it counts
            % against maxJobNum until all of them are done, not just one
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(4);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop');
            % the job of batch 1 (tasks 1 and 2) writes the output of task
            % 1 at the first poll and is gone at the fourth
            slurm.setJobStates('1000_1', {'R', 'R', 'R', 'done'}, outFns(1));

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''maxJobNum'', 1, ''taskBatchNum'', 2, ', ...
                '''queryInterval'', 0, ''jobLogDir'', dataPath)']);

            cmds = slurm.history();
            submit2 = find(startsWith(cmds, 'sbatch') & contains(cmds, '--array=2 '));
            testCase.assertNumElements(submit2, 1, 'Batch 2 should be submitted once.');
            polls = find(startsWith(cmds, 'squeue') & contains(cmds, '1000_1'));
            testCase.assertGreaterThanOrEqual(numel(polls), 4, 'The job of batch 1 was not polled until it ended.');
            testCase.verifyGreaterThan(submit2, polls(4), ...
                'Batch 2 was submitted while the job of batch 1 was still running.');
        end

        function clusterModeFinishedShortLastBatchFreesItsJobSlot(testCase)
            % a batch whose tasks are all done stops counting against
            % maxJobNum even while its job is still listed (e.g., exiting).
            % Here the batch is the last one and shorter than taskBatchNum;
            % the count pads it to taskBatchNum tasks, which must count as done.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(3);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop');
            % Slurm rejects batch 1 (tasks 1 and 2) and the framework's
            % immediate retry, so batch 2 (task 3) is submitted first; its
            % job writes the output at the first poll and stays listed
            slurm.rejectSubmissions(2);
            slurm.setJobStates('1000_2', {'R'}, outFns(3));
            % a framework that holds batch 1 back forever is stopped by a
            % watchdog that writes the final output
            finalFn = [dataPath, 'final.txt'];
            watchdog = timer('StartDelay', 20, 'TimerFcn', @(~, ~) fclose(fopen(finalFn, 'w')));
            testCase.addTeardown(@() delete(watchdog));
            testCase.addTeardown(@() stop(watchdog));
            start(watchdog);

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''maxJobNum'', 1, ''taskBatchNum'', 2, ', ...
                '''queryInterval'', 0, ''jobLogDir'', dataPath, ''finalOutFullpath'', finalFn)']);

            testCase.verifyFalse(exist(finalFn, 'file') == 2, ...
                'The finished batch 2 kept counting as a job until the watchdog stopped the framework.');
            arrays = cellfun(@(x) str2double(regexp(x, '--array=(\d+)', 'tokens', 'once')), slurm.calls('sbatch'));
            testCase.verifyEqual(arrays(:)', [1, 1, 2, 1], ...
                'Batch 1 should be submitted again once batch 2 is done.');
        end

        function clusterModeRefusedTasksAreSavedAsTheyAre(testCase)
            % after a refused submission the framework saves the task to a
            % file and submits a job that runs it from there; a bash task
            % was written as fprintf's format (its % and \ taken as
            % conversions and escapes, e.g. of a file name with ' quoted as
            % '\''), and a MATLAB task with its double quotes escaped for
            % the shell
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns] = testCase.makeTasks(1);
            task = ['printf ''%d\n'' 5 > ''', outFns{1}, '''; echo ''it''\''''s'''];
            slurm.rejectSubmissions(1);
            testCase.runCluster(inFns, outFns, {task}, 'language', 'bash', 'maxTrialNum', 1, 'queryInterval', 0);
            taskFiles = dir([fileparts(outFns{1}), '/func_strs/*.sh']);
            testCase.assertNumElements(taskFiles, 1, 'The task was not saved.');
            testCase.verifyEqual(fileread(fullfile(taskFiles(1).folder, taskFiles(1).name)), task);
            calls = slurm.calls('sbatch');
            testCase.verifyNotEmpty(regexp(calls{end}, ['bash \S+/', taskFiles(1).name], 'once'), calls{end});

            slurm.rejectSubmissions(1);
            [inFns, outFns] = testCase.makeTasks(1);
            task = 'fprintf("%d\n", 5)';
            testCase.runCluster(inFns, outFns, {task}, 'maxTrialNum', 1, 'queryInterval', 0);
            taskFiles = dir([fileparts(outFns{1}), '/func_strs/*.mat']);
            testCase.assertNumElements(taskFiles, 1, 'The task was not saved.');
            a = load(fullfile(taskFiles(1).folder, taskFiles(1).name), 'func_str');
            testCase.verifyEqual(a.func_str, task);
        end

        function clusterModeFailedQueryKeepsJobStates(testCase)
            % a failed squeue (e.g., a timeout of a busy controller) read as
            % all jobs gone: they stopped counting against maxJobNum (and
            % were submitted again after 45 s); they keep their states now
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop');
            % the first two polls fail; the job of batch 1 then runs and
            % writes its output
            slurm.failQueries(2);
            slurm.setJobStates('1000_1', {'R', 'R', 'done'}, outFns(1));

            testCase.runCluster(inFns, outFns, funcStrs, 'maxTrialNum', 1, 'maxJobNum', 1, 'queryInterval', 0);

            cmds = slurm.history();
            submit2 = find(startsWith(cmds, 'sbatch') & contains(cmds, '--array=2 '));
            testCase.assertNumElements(submit2, 1, 'Batch 2 should be submitted once.');
            polls = find(startsWith(cmds, 'squeue') & contains(cmds, '1000_1'));
            testCase.assertGreaterThanOrEqual(numel(polls), 3, 'The job of batch 1 was not polled until it ran.');
            testCase.verifyGreaterThan(submit2, polls(3), ...
                'Batch 2 was submitted while squeue failed to answer about the job of batch 1.');
        end

        function clusterModeLostSubmissionAnswerFindsTheJob(testCase)
            % sbatch failed with "Socket timed out" after Slurm made the job,
            % and the framework submitted the batch again at once (with the
            % task in a file): two jobs, one unknown to it. It now looks for
            % the job by its comment and uses it; while squeue fails (here
            % the first two times) the batch waits
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn] = testCase.makeTasks(1);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop');
            slurm.loseSubmissionAnswers(1);
            slurm.failQueries(2);
            slurm.setJobStates('1000_1', {'R'}, outFns(1));

            flags = testCase.runCluster(inFns, outFns, funcStrs, 'maxTrialNum', 1, 'queryInterval', 0);

            testCase.verifyEqual(flags, true);
            testCase.verifyNumElements(slurm.calls('sbatch'), 1, 'The batch was submitted again.');
            cmds = slurm.history();
            lookup = find(startsWith(cmds, 'squeue -u') & contains(cmds, '%k'), 1);
            testCase.assertNotEmpty(lookup, strjoin(cmds, newline));
            testCase.verifyTrue(any(strcmp(cmds(lookup + 1 : end), 'squeue -j 1000_1')), ...
                'The job found was not polled.');
            % cancelled once at the end (in the background)
            pause(1);
            testCase.verifyEqual(slurm.calls('scancel'), {'1000'});
        end

        function clusterModeLostSubmissionWithoutJobIsSubmittedAgain(testCase)
            % a submission whose answer was lost and of which Slurm lists no
            % job (looked for until 30 s after it) is submitted again; until
            % then it counts against maxJobNum
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop');
            slurm.rejectSubmissions(1, 'sbatch: error: Batch job submission failed: Socket timed out on send/recv operation');
            slurm.setJobStates('1000_1', {'R'}, outFns(1));
            slurm.setJobStates('1001_2', {'R'}, outFns(2));

            t0 = tic;
            flags = testCase.runCluster(inFns, outFns, funcStrs, 'maxTrialNum', 1, 'maxJobNum', 1, 'queryInterval', 1);

            testCase.verifyEqual(flags, true(2, 1));
            cmds = slurm.history();
            arrays = cellfun(@(x) str2double(regexp(x, '--array=(\d+)', 'tokens', 'once')), slurm.calls('sbatch'));
            testCase.assertEqual(arrays(:)', [1, 1, 2], strjoin(cmds, newline));
            submits = find(startsWith(cmds, 'sbatch'));
            lookups = find(startsWith(cmds, 'squeue -u'));
            testCase.verifyTrue(any(lookups > submits(1) & lookups < submits(2)), ...
                'The batch was submitted again before its job was looked for.');
            testCase.verifyGreaterThanOrEqual(toc(t0), 30, 'The batch was submitted again before 30 s.');
        end

        function clusterModeFailingMasterCancelsItsJobs(testCase)
            % a master whose own batch fails twice raises; the call's jobs
            % stayed in the queue
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn] = testCase.makeTasks(2);
            funcStrs = [testCase.taskStrings(inFns(1), outFns(1), logFn, 'noop'), ...
                testCase.taskStrings(inFns(2), outFns(2), logFn, 'error')];
            slurm.setJobStates('1000_1', {'PD'});

            % (the master runs the last batch; memAllocate: its task fits the
            % test's job)
            testCase.verifyError(@() testCase.runCluster(inFns, outFns, funcStrs, 'masterCompute', true, ...
                'memAllocate', 1, 'queryInterval', 0), 'pktest:frameworkTask:failed');

            testCase.verifyEqual(slurm.calls('scancel'), {'1000'}, 'The job of batch 1 was left in the queue.');
        end

        function clusterModePausesWhileOnlyLastTriesRun(testCase)
            % while only jobs of last tries were listed, the framework looped
            % without a pause (a busy core) until they ended
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn] = testCase.makeTasks(1);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'noop');
            slurm.setJobStates('1000_1', {'R', 'R', 'R', 'done'});
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() builtin('pause', 0.05)));

            testCase.runCluster(inFns, outFns, funcStrs, 'maxTrialNum', 1, 'queryInterval', 0.2);

            testCase.verifyGreaterThanOrEqual(numel(pauses.getCalls()), 3, ...
                'The framework did not pause while the job of the last try ran.');
        end

        function missingConfigFileErrors(testCase)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');
            testCase.verifyError(@() generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'ConfigFile', [dataPath, 'missing.json']), ?MException);
            testCase.verifyEmpty(testCase.readLog(logFn));
        end

        function masterParaFactorAcceptsFractions(testCase)
            % masterParaFactor is the fraction of parallel jobs on the master
            % in mcc mode; it used to be validated with @islogical
            [inFns, outFns, logFn] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true, 'masterParaFactor', 0.5);

            testCase.verifyEqual(flags, true(2, 1));
        end

        function callLogWaitsOutRefusedFileOperations(testCase)
            % Windows refuses for a moment to make the call log's lock folder
            % (while another task removes it), to open the log and to remove
            % the lock; pktest.frameworkTask tries again (it took the lock
            % without making it and then failed to remove it, which failed
            % tasks on GitHub's Windows runners)
            [inFns, outFns, logFn] = testCase.makeTasks(1);
            stubPath = testCase.createTempFolder();
            testCase.writeRefusingStub(stubPath, 'mkdir', '.lock', ...
                '{false, ''Access is denied.'', ''MATLAB:MKDIR:OSError''}');
            testCase.writeRefusingStub(stubPath, 'fopen', '.log', '{-1, ''Permission denied''}');
            testCase.writeRefusingStub(stubPath, 'rmdir', '.lock', ...
                '{false, ''Access is denied.'', ''MATLAB:RMDIR:OSError''}');
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture('MATLAB:dispatcher:nameConflict'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(stubPath));

            pktest.frameworkTask(inFns{1}, outFns{1}, logFn, 'noop');

            testCase.verifyEqual(testCase.readLog(logFn), outFns(1));
            testCase.verifyFalse(isfolder([logFn, '.lock']), 'The lock folder was left behind.');
            for name = {'mkdir', 'fopen', 'rmdir'}
                testCase.verifyTrue(isfile([stubPath, name{1}, '_refused']), [name{1}, ' was not refused.']);
            end
        end
    end

    methods (Test, TestTags = {'Integration', 'Slow'})
        function parforWrapperRunsEveryTask(testCase)
            % starts a local pool with 2 workers. parseParfor used to fail
            % on the first line, which looked for the root path in
            % 'XiongtaoScripts'.
            [inFns, outFns, logFn] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');
            testCase.addTeardown(@() delete(gcp('nocreate')));

            flags = matlab_parfor_generic_computing_wrapper(inFns, outFns, funcStrs, 'GPUJob', false, 'nworker', 2);

            testCase.verifyTrue(all(flags));
            for i = 1 : 2
                testCase.verifyEqual(fileread(outFns{i}), fileread(inFns{i}));
            end
        end

        function parforGPUJobHasAWorkerPerUsableGPU(testCase)
            % GPUs newer than this MATLAB supports (e.g., B200 GPUs in
            % R2024b) weren't counted without forward compatibility, so a
            % GPU job asked for a pool of 0 workers. usable_gpu_num and
            % gpuDevice are stubs (2 GPUs; the workers select them too).
            testCase.applyFixture(pktest.FunctionStubFixture('usable_gpu_num', 2));
            testCase.applyFixture(pktest.FunctionStubFixture('gpuDevice', struct('Index', 1)));
            [inFns, outFns, logFn] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');
            delete(gcp('nocreate'));
            testCase.addTeardown(@() delete(gcp('nocreate')));

            out = evalc('flags = matlab_parfor_generic_computing_wrapper(inFns, outFns, funcStrs, ''GPUJob'', true);');

            testCase.verifyTrue(all(flags));
            pool = gcp('nocreate');
            testCase.assertNotEmpty(pool, 'No pool was started.');
            testCase.verifyEqual(pool.NumWorkers, 2);
            % each worker selected a GPU
            testCase.verifyEqual(numel(strfind(out, 'Using GPU')), 2, out);
            for i = 1 : 2
                testCase.verifyEqual(fileread(outFns{i}), fileread(inFns{i}));
            end
        end
    end

    methods (Test, TestTags = {'Integration', 'Slow'})
        function parforBashWorkersUseTheJobsGPUs(testCase)
            % the bash tasks of a parfor GPU job got GPU k - 1 for worker k,
            % other GPUs than a Slurm job's own where they aren't 0 to n - 1;
            % they get the k-th GPU in its CUDA_VISIBLE_DEVICES now, with no
            % more workers than it lists (usable_gpu_num and gpuDevice are
            % stubs of 2 GPUs; no GPU is used)
            testCase.assumeFalse(ispc, 'The tasks are POSIX shell commands.');
            testCase.applyFixture(pktest.FunctionStubFixture('usable_gpu_num', 2));
            testCase.applyFixture(pktest.FunctionStubFixture('gpuDevice', struct('Index', 1)));
            delete(gcp('nocreate'));
            testCase.addTeardown(@() delete(gcp('nocreate')));
            % GPU lists, the GPUs of the tasks and the workers
            lists = {'4,6', {'4', '6'}, 2; '5', {'5'}, 1};
            for i = 1 : size(lists, 1)
                testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('CUDA_VISIBLE_DEVICES', lists{i, 1}));
                [inFns, outFns] = testCase.makeTasks(2);
                funcStrs = cellfun(@(x) sprintf('sh -c ''echo $CUDA_VISIBLE_DEVICES > %s''', x), outFns, 'unif', 0); %#ok<NASGU>
                flags = []; %#ok<NASGU>
                evalc(['flags = matlab_parfor_generic_computing_wrapper(inFns, outFns, funcStrs, ''GPUJob'', true, ', ...
                    '''language'', ''bash'');']);
                testCase.verifyTrue(all(flags));
                gpus = cellfun(@(x) strtrim(fileread(x)), outFns(cellfun(@isfile, outFns)), 'unif', 0);
                msg = sprintf('CUDA_VISIBLE_DEVICES=%s: %s', lists{i, 1}, strjoin(gpus, ', '));
                testCase.verifyTrue(all(ismember(gpus, lists{i, 2})), msg);
                pool = gcp('nocreate');
                testCase.assertNotEmpty(pool, 'No pool was started.');
                testCase.verifyEqual(pool.NumWorkers, lists{i, 3}, msg);
                delete(pool);
            end
        end
    end

    methods (Access = private)
        function [inFns, outFns, logFn, dataPath] = makeTasks(testCase, n)
            % n input text files and the output paths in a separate folder
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
        end

        function [flags] = runCluster(~, inFns, outFns, funcStrs, varargin)
            % the framework in cluster mode with a fake Slurm (output
            % hidden), with the jobs' logs in the tasks' folder
            dataPath = [fileparts(fileparts(outFns{1})), '/'];
            evalc(['flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''jobLogDir'', dataPath, varargin{:});']);
        end

        function [funcStrs] = taskStrings(~, inFns, outFns, logFn, mode)
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''%s'')', ...
                x, y, logFn, mode), inFns, outFns, 'unif', 0);
        end

        function [] = writeRefusingStub(~, stubPath, name, suffix, refusal)
            % a stub of the built-in name that refuses its first call for a
            % path ending in suffix (returns refusal, or errors without
            % outputs, and leaves a file name_refused) and runs it otherwise
            lines = {
                sprintf('function varargout = %s(varargin)', name)
                'persistent refused'
                sprintf('if isempty(refused) && nargin > 0 && (ischar(varargin{1}) || isstring(varargin{1})) && endsWith(varargin{1}, ''%s'')', suffix)
                '    refused = true;'
                sprintf('    fclose(builtin(''fopen'', [fileparts(mfilename(''fullpath'')), ''/%s_refused''], ''w''));', name)
                '    if nargout == 0'
                sprintf('        error(''pktest:refused'', ''%s refused'');', name)
                '    end'
                sprintf('    varargout = %s;', refusal)
                '    return;'
                'end'
                sprintf('[varargout{1 : nargout}] = builtin(''%s'', varargin{:});', name)
                'end'};
            fid = fopen([stubPath, name, '.m'], 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
        end

        function [lines] = readLog(~, logFn)
            % output paths of all task calls, in call order
            if ~exist(logFn, 'file')
                lines = cell(0, 1);
                return;
            end
            lines = splitlines(strip(fileread(logFn)));
            lines = lines(~cellfun(@isempty, lines));
        end
    end
end
