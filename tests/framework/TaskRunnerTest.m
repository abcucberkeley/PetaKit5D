classdef TaskRunnerTest < pktest.TestCase
    % task runners of the job framework (task_runner, task_runner_run,
    % task_runner.sh): a process runs one task after another, its batch's
    % first, then unstarted tasks of the call's other batches, with claims
    % in the call folder. The job script runs with a stand-in compiled
    % master (a bash script that runs a fake task and, as a runner, claims
    % and runs the next tasks as task_runner does); the MATLAB runner runs
    % in this MATLAB with mccMaster replaced by a recorder; the framework
    % talks to a fake Slurm (pktest.FakeSlurmFixture).

    properties (Constant)
        EnvNames = {'PETAKIT5D_RUNNER', 'PETAKIT5D_TASK_METRICS', 'PETAKIT5D_TASK', 'PETAKIT5D_TASK_GPU', ...
            'PK_RUNNER_STATE', 'PETAKIT5D_RUNNER_MEM_GB', 'PETAKIT5D_RUNNER_MINUTES', 'SLURM_JOB_END_TIME', ...
            'SLURM_JOB_START_TIME'}
    end

    methods (TestMethodSetup)
        function clearCaches(testCase)
            % config files are cached in persistent variables
            clear generic_computing_frameworks_wrapper
        end

        function restoreEnvironment(testCase)
            for i = 1 : numel(testCase.EnvNames)
                testCase.addTeardown(@setenv, testCase.EnvNames{i}, getenv(testCase.EnvNames{i}));
                setenv(testCase.EnvNames{i}, '');
            end
        end
    end

    methods (Test, TestTags = {'Unit'})
        function runnersAreOnlyUsedWithBuildsThatHaveThem(testCase)
            % a build without task_runner.txt (made before runners) or with
            % another runner version runs as before; runners need Linux
            dataPath = testCase.createTempFolder();
            masterFn = [dataPath, 'run_mccMaster.sh'];
            testCase.writeLines(masterFn, {'#!/bin/sh'});
            testCase.verifyFalse(task_runner_run('supported', masterFn));
            testCase.writeLines([dataPath, 'task_runner.txt'], {sprintf('%d', task_runner('version') + 1)});
            testCase.verifyFalse(task_runner_run('supported', masterFn));
            testCase.writeLines([dataPath, 'task_runner.txt'], {sprintf('%d', task_runner('version'))});
            testCase.verifyEqual(task_runner_run('supported', masterFn), isunix && ~ismac);
        end
    end

    methods (Test, TestTags = {'Unit'})
        function jobRunsEachTaskOnceInFewerProcesses(testCase)
            testCase.assumeLinux();
            [c, b] = testCase.makeCall(8, 8, {});
            out = testCase.runJob(c, 'job 1 4 0 1 2');
            runs = testCase.taskRuns(b);
            testCase.verifyEqual(sort([runs.task]), 1 : 8, out);
            testCase.verifyTrue(all(cellfun(@isfile, c.outFns)));
            testCase.verifyLessThanOrEqual(numel(unique([runs.pid])), 4, 'More processes than slots.');
            testCase.verifyLessThan(numel(unique([runs.pid])), 8, 'Each task ran in a process of its own.');
            testCase.verifyEqual(sort(testCase.claimNames(c)), sort(arrayfun(@(x) sprintf('t%d', x), 1 : 8, 'unif', 0)));
        end

        function jobTakesUnstartedTasksOfOtherBatches(testCase)
            testCase.assumeLinux();
            [c, b] = testCase.makeCall(12, 4, {});
            testCase.runJob(c, 'job 1 2 0 1 2');
            runs = testCase.taskRuns(b);
            testCase.verifyEqual(sort([runs.task]), 1 : 12);
            % without stealing, only its own batch
            [c, b] = testCase.makeCall(12, 4, {});
            testCase.runJob(c, 'job 2 2 0 0 2');
            runs = testCase.taskRuns(b);
            testCase.verifyEqual(sort([runs.task]), 5 : 8);
        end

        function jobLeavesClaimedTasksResubmittedBatchesAndOutputs(testCase)
            % a task another job claimed, a task whose output exists, and a
            % batch the master resubmitted (its job runs it)
            testCase.assumeLinux();
            [c, b] = testCase.makeCall(12, 4, {});
            system(sprintf('ln -s "99_1:othernode:1:1" "%s/t3"', c.rc.claimsDir));
            testCase.writeLines(c.outFns{4}, {'x'});
            testCase.writeLines([c.rc.claimsDir, '/r2'], {});
            out = testCase.runJob(c, 'job 1 1 0 1 1');
            runs = testCase.taskRuns(b);
            testCase.verifyEqual(sort([runs.task]), [1, 2, 9 : 12]);
            testCase.verifySubstring(out, sprintf('output %s already exists', c.outFns{4}));
            testCase.verifyEqual(fileread(c.outFns{4}), sprintf('x\n'), 'The existing output was rewritten.');
        end

        function tasksWithoutInputsStayInTheirBatch(testCase)
            % a task whose input was missing when the call started runs only
            % in its own batch's job (as before), not in another job's runner
            testCase.assumeLinux();
            [c, b] = testCase.makeCall(4, 2, {}, '', [false; false; false; true]);
            testCase.verifyEqual(sscanf(fileread([c.rc.D, '/runnable.txt']), '%d')', 1 : 3);
            testCase.runJob(c, 'job 1 1 0 1 2');
            testCase.verifyEqual(sort([testCase.taskRuns(b).task]), 1 : 3);
        end

        function batchesWaitThirtySecondsBeforeRunnersTakeThem(testCase)
            % the master lists a batch for the runners once its tasks have
            % waited 30 s in the queue (or for their submission)
            [~, outFns] = testCase.makeTasks(4);
            evalc(['rc = task_runner_run(''setup_matlab'', [testCase.createTempFolder(), ''call''], outFns, ', ...
                'outFns, '''', false(4, 1), 2, ''x'', '''', false);']);
            waiting = [false; false; true; true];
            listed = @(rc) strtrim(fileread_or_empty([rc.D, '/waiting.txt']));
            rc = task_runner_run('waiting', rc, waiting, 100, 2);
            testCase.verifyEmpty(listed(rc));
            rc = task_runner_run('waiting', rc, waiting, 129, 2);
            testCase.verifyEmpty(listed(rc));
            rc = task_runner_run('waiting', rc, waiting, 131, 2);
            testCase.verifyEqual(listed(rc), '2');
            % a batch whose job started leaves the list
            rc = task_runner_run('waiting', rc, false(4, 1), 140, 2);
            testCase.verifyEmpty(listed(rc));
        end

        function jobTakesOnlyTasksOfBatchesWithItsKey(testCase)
            % a batch's key: the CPUs and GPU of its job's processes, or
            % 'master' for the batch the master runs itself; a runner takes
            % other batches' tasks only if they have the key of its own
            testCase.assumeLinux();
            testCase.verifyEqual(task_runner_run('keys', 3, 4, false, 3), {'4'; '4'; 'master'});
            testCase.verifyEqual(task_runner_run('keys', 2, 8, true, 0), {'8gpu'; '8gpu'});
            [c, b] = testCase.makeCall(16, 4, {});
            testCase.writeLines([c.rc.D, '/keys.txt'], {'4'; '4'; 'master'; '8'});
            testCase.runJob(c, 'job 1 1 0 1 2');
            testCase.verifyEqual(sort([testCase.taskRuns(b).task]), 1 : 8);
            % the master's batch: its own tasks only
            [c, b] = testCase.makeCall(16, 4, {});
            testCase.writeLines([c.rc.D, '/keys.txt'], {'4'; '4'; 'master'; '8'});
            testCase.runJob(c, 'job 3 1 0 1 2');
            testCase.verifyEqual(sort([testCase.taskRuns(b).task]), 9 : 12);
        end

        function failedTaskRunsAgainInItsOwnProcessAndIsMarked(testCase)
            % a task error and a crash: each task runs once more in a
            % process of its own (as GNU parallel --retries 2 did), then it
            % is marked failed; the other tasks run, and the slot takes no
            % more tasks of other batches
            testCase.assumeLinux();
            [c, b] = testCase.makeCall(12, 8, {2, 'error'; 5, 'crash'});
            testCase.runJob(c, 'job 1 1 0 1 2');
            runs = testCase.taskRuns(b);
            tasks = [runs.task];
            testCase.verifyEqual(sum(tasks == 2), 2);
            testCase.verifyEqual(sum(tasks == 5), 2);
            testCase.verifyEqual(sort(tasks(tasks ~= 2 & tasks ~= 5)), [1, 3, 4, 6, 7, 8]);
            names = testCase.claimNames(c);
            testCase.verifyTrue(all(ismember({'f2', 'f5'}, names)));
            testCase.verifyFalse(any(ismember({'f1', 'f3'}, names)));
            % retries 1 (sequential jobs): one try, then marked
            [c, b] = testCase.makeCall(4, 4, {2, 'error'});
            testCase.runJob(c, 'job 1 1 0 0 1');
            tasks = [testCase.taskRuns(b).task];
            testCase.verifyEqual(sum(tasks == 2), 1);
            testCase.verifyTrue(ismember('f2', testCase.claimNames(c)));
        end

        function jobWithNoTaskLeftStartsNoProcess(testCase)
            testCase.assumeLinux();
            [c, b] = testCase.makeCall(4, 4, {});
            for f = 1 : 4
                system(sprintf('ln -s "99_1:othernode:1:1" "%s/t%d"', c.rc.claimsDir, f));
            end
            out = testCase.runJob(c, 'job 1 2 0 1 2');
            testCase.verifyEmpty(testCase.taskRuns(b));
            testCase.verifySubstring(out, 'no task left');
        end

        function callFolderInTheHomeFolderAsTildeRunsItsTasks(testCase)
            % a call folder given as ~/... (data paths such as the demos'
            % ~/Downloads/...): the job script finds the claims and tasks
            % (bash doesn't expand ~ in the script's quoted paths)
            testCase.assumeLinux();
            home = testCase.createTempFolder();
            testCase.addTeardown(@setenv, 'HOME', getenv('HOME'));
            setenv('HOME', regexprep(home, '/$', ''));
            [c, b] = testCase.makeCall(4, 4, {}, '', false(4, 1), '~/data/tmp/call/');
            testCase.verifyEqual(c.rc.D, [home, 'data/tmp/call']);
            out = testCase.runJob(c, 'job 1 2 0 1 2');
            testCase.verifyEqual(sort([testCase.taskRuns(b).task]), 1 : 4, out);
            testCase.verifyTrue(all(cellfun(@isfile, c.outFns)));
            testCase.verifyEqual(sort(testCase.claimNames(c)), sort(arrayfun(@(x) sprintf('t%d', x), 1 : 4, 'unif', 0)));
        end

        function matlabJobWithNoTaskLeftStartsNoMatlab(testCase)
            % a MATLAB job starts MATLAB only for a task it claimed
            testCase.assumeLinux();
            [~, outFns] = testCase.makeTasks(2);
            for claimed = [false, true]
                evalc(['rc = task_runner_run(''setup_matlab'', [testCase.createTempFolder(), ''call''], outFns, ', ...
                    'outFns, '''', false(2, 1), 2, ''echo STARTED'', '''', false);']);
                if claimed
                    for f = 1 : 2
                        system(sprintf('ln -s "99_1:othernode:1:1" "%s/t%d"', rc.claimsDir, f));
                    end
                end
                [status, out] = system(sprintf('bash "%s" matlab 1 1', rc.script));
                testCase.assertEqual(status, 0, out);
                testCase.verifyEqual(contains(out, 'STARTED'), ~claimed, out);
            end
            testCase.verifySubstring(out, 'no task left');
        end

        function matlabJobMarksTasksItDidNotFinish(testCase)
            % a task a MATLAB job's process didn't finish is marked failed
            % (the master resubmits its batch with more CPUs), and the job
            % then takes no tasks of other batches
            testCase.assumeLinux();
            [~, outFns] = testCase.makeTasks(4);
            evalc(['rc = task_runner_run(''setup_matlab'', [testCase.createTempFolder(), ''call''], outFns, ', ...
                'outFns, '''', false(4, 1), 2, ''true'', '''', false);']);
            [status, out] = system(sprintf('bash "%s" matlab 1 1', rc.script));
            testCase.assertEqual(status, 0, out);
            names = dir(rc.claimsDir);
            names = setdiff({names.name}, {'.', '..'});
            testCase.verifyTrue(all(ismember({'f1', 'f2'}, names)), strjoin(names, ' '));
            testCase.verifyFalse(any(ismember({'t3', 't4'}, names)), 'The job took tasks of another batch after a failure.');
        end

        function jobStartsNoMoreSlotsThanTasks(testCase)
            % as GNU parallel started one process per task line, a job starts
            % no more slots than it has tasks to take (each slot starts 0.2 s
            % after the one before, and the job ends after the last one)
            testCase.assumeLinux();
            testCase.assumeEqual(system('command -v parallel > /dev/null'), 0, 'GNU parallel is not installed.');
            [c, b] = testCase.makeCall(2, 2, {});
            t0 = tic;
            testCase.runJob(c, 'job 1 12 0 0 2');
            t = toc(t0);
            testCase.verifyEqual(sort([testCase.taskRuns(b).task]), 1 : 2);
            testCase.verifyLessThan(t, 1.5, 'The job started all 12 slots for 2 tasks.');
        end

        function processWithManyMoreOpenFilesIsRecycled(testCase)
            % before a runner's next task, the job script compares its open
            % files with those after its first task: 64 more end the process
            testCase.assumeLinux();
            c = testCase.makeCall(4, 4, {});
            L = testCase.createTempFolder();
            % a process that opens 100 files after a second
            [~, pid] = system(['bash -c ''sleep 1; for i in $(seq 100); do eval "exec $((i + 10))< /dev/null"; done; ', ...
                'sleep 30'' > /dev/null 2>&1 & echo $!']);
            pid = strtrim(pid);
            testCase.addTeardown(@() system(['kill ', pid, ' 2> /dev/null']));
            claim = sprintf('PK_PID=%s bash "%s" claim 1 0 1 0 "%ss1" "%s" "%snext"', pid, c.rc.script, L, L(1 : end - 1), L);
            system(claim);
            testCase.verifyEqual(strtrim(fileread([L, 'next'])), '1');
            pause(2);
            system(claim);
            testCase.verifySubstring(fileread([L, 'next']), 'recycle');
        end

        function gpuSlotsGetOneDeviceEach(testCase)
            testCase.assumeLinux();
            testCase.assumeEqual(system('command -v parallel > /dev/null'), 0, 'GNU parallel is not installed.');
            [c, b] = testCase.makeCall(4, 4, {1, 'slow'; 2, 'slow'; 3, 'slow'; 4, 'slow'});
            testCase.runJob(c, 'job 1 2 1 0 2');
            runs = testCase.taskRuns(b);
            testCase.verifyEqual(sort([runs.task]), 1 : 4);
            testCase.verifyEqual(sort(unique({runs.cuda})), {'0', '1'});
            % each process sees the device of its slot
            for p = unique([runs.pid])
                testCase.verifyNumElements(unique({runs([runs.pid] == p).cuda}), 1);
            end
        end

        function stealingStopsNearTheJobsEnd(testCase)
            % a job takes tasks of other batches only if it has the time for
            % them (twice the longest task so far and a minute); its own
            % tasks run anyway
            testCase.assumeLinux();
            [c, b] = testCase.makeCall(8, 4, {});
            now_ = floor(posixtime(datetime('now', 'TimeZone', 'local')));
            setenv('SLURM_JOB_START_TIME', sprintf('%d', now_ - 3600));
            setenv('SLURM_JOB_END_TIME', sprintf('%d', now_ + 30));
            testCase.runJob(c, 'job 1 1 0 1 2');
            testCase.verifyEqual(sort([testCase.taskRuns(b).task]), 1 : 4);
            setenv('SLURM_JOB_END_TIME', sprintf('%d', now_ + 3600));
            [c, b] = testCase.makeCall(8, 4, {});
            testCase.runJob(c, 'job 1 1 0 1 2');
            testCase.verifyEqual(sort([testCase.taskRuns(b).task]), 1 : 8);
        end

        function claimsAreTakenOnceUnderConcurrency(testCase)
            % many runners claiming at once: each task once
            testCase.assumeLinux();
            nT = 400;
            nP = 16;
            [c, b] = testCase.makeCall(nT, nT, {});
            L = testCase.createTempFolder();
            cmd = sprintf(['for s in $(seq %d); do (for k in $(seq %d); do ', ...
                'bash "%s" claim 1 0 $s 0 "%ss$s" "%s" "%snext_$s" > /dev/null; ', ...
                'f=$(cat "%snext_$s"); [ -z "$f" ] && break; echo "$f" >> "%sgot_$s"; done) & done; wait'], ...
                nP, nT, c.rc.script, L, L(1 : end - 1), L, L, L);
            [status, out] = system(cmd);
            testCase.assertEqual(status, 0, out);
            got = [];
            for s = 1 : nP
                fn = sprintf('%sgot_%d', L, s);
                if isfile(fn)
                    got = [got; sscanf(fileread(fn), '%d')]; %#ok<AGROW>
                end
            end
            testCase.verifyEqual(sort(got)', 1 : nT, 'Each task should be claimed exactly once.');
            testCase.verifyEmpty(testCase.taskRuns(b), 'Claiming ran tasks.');
        end

        function localRunsSkipTasksThatOtherCallsProcess(testCase)
            % local runs (no cluster) keep the .tmp files with which
            % framework calls on the same data share their tasks: a task
            % whose .tmp is younger than the wait time is left alone, an
            % older one is taken; .tmp files of finished tasks are removed
            testCase.assumeLinux();
            [c, b] = testCase.makeCall(4, 4, {}, 'local');
            testCase.writeLines(c.tmpFns{2}, {});
            testCase.writeLines(c.tmpFns{3}, {});
            system(sprintf('touch -d "-10 minutes" "%s"', c.tmpFns{3}));
            testCase.runJob(c, 'job 1 1 0 2 1');
            testCase.verifyEqual(sort([testCase.taskRuns(b).task]), [1, 3, 4]);
            testCase.verifyEqual(cellfun(@isfile, c.tmpFns(:))', [false, true, false, false]);
            testCase.verifyFalse(ismember('t2', testCase.claimNames(c)), 'The task left to another call is claimed.');
        end

        function localRunsRunATaskARunnerGaveBack(testCase)
            % local runs: a task a runner gives back (e.g. its memory guard)
            % loses the .tmp of its claim and runs in a new process at once
            testCase.assumeLinux();
            [c, b] = testCase.makeCall(4, 4, {2, 'giveback'}, 'local');
            testCase.runJob(c, 'job 1 1 0 2 1');
            testCase.verifyEqual(sort([testCase.taskRuns(b).task]), 1 : 4);
            testCase.verifyTrue(all(cellfun(@isfile, c.outFns)));
            testCase.verifyFalse(any(cellfun(@isfile, c.tmpFns)), 'A .tmp file is left.');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function runnerGivesMccMasterTheArgumentsOfTheTaskLine(testCase)
            % the runner must give mccMaster what run_mccMaster.sh gives
            % the compiled mccMaster for the same task line: paths with
            % spaces, cell arrays, empty strings, numbers, non-ASCII
            % characters, final outputs
            testCase.assumeLinux();
            dataPath = testCase.createTempFolder();
            [lines, taskFn, master] = testCase.frameworkTaskLines(dataPath, 'mixed');
            nT = numel(lines);
            % the arguments the build's run_mccMaster.sh gives a stand-in
            % mccMaster that writes them (NUL separated)
            build = testCase.argsBuild([dataPath, 'argsbuild/']);
            expected = cell(nT, 1);
            for i = 1 : nT
                argFn = sprintf('%sargs_%d.bin', dataPath, i);
                line = strrep(lines{i}, master, build.master);
                testCase.writeLines([dataPath, 'line.sh'], {line});
                [status, out] = system(sprintf('PK_ARGS_OUT="%s" bash "%sline.sh"', argFn, dataPath));
                testCase.assertEqual(status, 0, out);
                expected{i} = testCase.readNul(argFn);
            end
            % the runner, with mccMaster replaced by a recorder
            recorded = testCase.runInProcess(taskFn, master, 1 : nT);
            testCase.assertNumElements(recorded, nT);
            for i = 1 : nT
                testCase.verifyEqual(recorded{i}, expected{i}, sprintf('Arguments of task %d differ.', i));
            end
        end

        function runnerResetsTheProcessBetweenTasks(testCase)
            % each task sees the state the process had before its first one
            testCase.assumeLinux();
            dataPath = testCase.createTempFolder();
            n = 4;
            outFns = arrayfun(@(x) sprintf('%sstate_%d.mat', dataPath, x), 1 : n, 'unif', 0);
            master = testCase.stateTaskLines(dataPath, outFns, repmat({'ok'}, 1, n));
            nThreads = maxNumCompThreads;
            testCase.addTeardown(@maxNumCompThreads, nThreads);
            testCase.addTeardown(@() evalin('base', 'clear pkRunnerStateBase'));
            testCase.addTeardown(@() clear('global', 'pkRunnerStateGlobal'));
            % a new process: default random numbers, no last warning
            rng('default');
            lastwarn('', '');
            testCase.runInProcess([dataPath, 'tasks.txt'], master, 1 : n, 'state');
            states = cellfun(@load, outFns);
            testCase.assertNumElements(states, n);
            first = states(1);
            testCase.verifyEqual(first.persistentCalls, 1);
            % (a task that starts compiled processes must not make them runners)
            testCase.verifyEmpty(first.runnerEnv, 'The task saw the runner settings.');
            for i = 2 : n
                s = states(i);
                testCase.verifyEqual(s, first, sprintf('Task %d saw another state than the first task.', i));
            end
        end

        function runnerRecyclesAProcessThatHoldsMemory(testCase)
            % after a task the process holds more memory than at its start
            % (beyond the margin): it ends, and a new process runs the next
            testCase.assumeLinux();
            testCase.addTeardown(@() testCase.removeHeldMemory());
            dataPath = testCase.createTempFolder();
            outFns = arrayfun(@(x) sprintf('%sstate_%d.mat', dataPath, x), 1 : 3, 'unif', 0);
            master = testCase.stateTaskLines(dataPath, outFns, {'hold', 'ok', 'ok'});
            setenv('PETAKIT5D_RUNNER_MEM_GB', '0.25');
            [~, status] = testCase.runInProcess([dataPath, 'tasks.txt'], master, 1 : 3, 'state');
            testCase.verifyTrue(isfile(outFns{1}));
            testCase.verifyFalse(isfile(outFns{2}), 'The process ran another task.');
            testCase.verifySubstring(status, 'x recycle');
        end

        function runnerRecyclesAProcessWhoseLaterTaskHoldsMemory(testCase)
            % a later task: compared with the process after its first task
            % and the reset (with the libraries the tasks load)
            testCase.assumeLinux();
            testCase.addTeardown(@() testCase.removeHeldMemory());
            dataPath = testCase.createTempFolder();
            outFns = arrayfun(@(x) sprintf('%sstate_%d.mat', dataPath, x), 1 : 3, 'unif', 0);
            master = testCase.stateTaskLines(dataPath, outFns, {'ok', 'hold', 'ok'});
            setenv('PETAKIT5D_RUNNER_MEM_GB', '0.25');
            [~, status] = testCase.runInProcess([dataPath, 'tasks.txt'], master, 1 : 3, 'state');
            testCase.verifyTrue(isfile(outFns{2}));
            testCase.verifyFalse(isfile(outFns{3}), 'The process ran another task.');
            testCase.verifySubstring(status, 'x recycle');
        end

        function runnerEndsAfterItsTime(testCase)
            % a process runs tasks for PETAKIT5D_RUNNER_MINUTES (10) at most
            testCase.assumeLinux();
            testCase.addTeardown(@() evalin('base', 'clear pkRunnerStateBase'));
            testCase.addTeardown(@() clear('global', 'pkRunnerStateGlobal'));
            dataPath = testCase.createTempFolder();
            outFns = arrayfun(@(x) sprintf('%sstate_%d.mat', dataPath, x), 1 : 2, 'unif', 0);
            master = testCase.stateTaskLines(dataPath, outFns, {'ok', 'ok'});
            setenv('PETAKIT5D_RUNNER_MINUTES', '0');
            [~, status] = testCase.runInProcess([dataPath, 'tasks.txt'], master, 1 : 2, 'state');
            testCase.verifyFalse(isfile(outFns{2}), 'The process ran another task.');
            testCase.verifyTrue(endsWith(strtrim(status), 'x recycle'), status);
        end

        function runnerEndsAfterATaskThatChangesTheThreadCount(testCase)
            % all tasks of a process run with its thread count (MATLAB fixes
            % at the first FFT how a process's FFTs use threads): a task that
            % changes it ends the process, and a new one runs the next task
            testCase.assumeLinux();
            nThreads = maxNumCompThreads;
            testCase.addTeardown(@maxNumCompThreads, nThreads);
            testCase.addTeardown(@() evalin('base', 'clear pkRunnerStateBase'));
            testCase.addTeardown(@() clear('global', 'pkRunnerStateGlobal'));
            dataPath = testCase.createTempFolder();
            outFns = arrayfun(@(x) sprintf('%sstate_%d.mat', dataPath, x), 1 : 3, 'unif', 0);
            master = testCase.stateTaskLines(dataPath, outFns, {'ok', 'threads', 'ok'});
            [~, status] = testCase.runInProcess([dataPath, 'tasks.txt'], master, 1 : 3, 'state');
            testCase.verifyTrue(isfile(outFns{2}));
            testCase.verifyFalse(isfile(outFns{3}), 'The process ran another task after the thread count changed.');
            testCase.verifyTrue(contains(status, 'd 2') && endsWith(strtrim(status), 'x recycle'), status);
        end

        function runnerStopsAfterAFailedTask(testCase)
            testCase.assumeLinux();
            dataPath = testCase.createTempFolder();
            outFns = arrayfun(@(x) sprintf('%sstate_%d.mat', dataPath, x), 1 : 3, 'unif', 0);
            master = testCase.stateTaskLines(dataPath, outFns, {'ok', 'error', 'ok'});
            [~, status] = testCase.runInProcess([dataPath, 'tasks.txt'], master, 1 : 3, 'state');
            testCase.verifyTrue(isfile(outFns{2}));
            testCase.verifyFalse(isfile(outFns{3}), 'The process ran another task after a failed one.');
            testCase.verifyTrue(contains(status, sprintf('e 2\nx failed')), status);
        end

        function matlabRunnerRunsItsBatchThenOtherBatches(testCase)
            % a MATLAB job's runner: its batch, then the other batches' tasks
            % from the last one; the base workspace keeps its variables
            testCase.assumeLinux();
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(6);
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''copy'')', x, y, logFn), ...
                inFns, outFns, 'unif', 0);
            evalc(['rc = task_runner_run(''setup_matlab'', [dataPath, ''call''], funcStrs, outFns, '''', ', ...
                'false(6, 1), 2, ''matlab -batch'', '''', false);']);
            testCase.assertTrue(rc.on);
            L = [testCase.createTempFolder(), 'local'];
            mkdir(L);
            st = [L, '/s1_1'];
            system(sprintf('ln -s "0_0:test:1:1" "%s/t1"', rc.claimsDir));
            testCase.writeLines(st, {'c 1'});
            assignin('base', 'pkRunnerKeep', 1);
            testCase.addTeardown(@() evalin('base', 'clear pkRunnerKeep'));
            setenv('PETAKIT5D_RUNNER_MEM_GB', '100');
            evalc('task_runner(''matlab'', rc.D, ''1'', 1, st, L, 1)');
            testCase.verifyEqual(testCase.readLog(logFn), outFns([1, 2, 6, 5, 4, 3])');
            testCase.verifyEqual(cellfun(@fileread, outFns, 'unif', 0), cellfun(@fileread, inFns, 'unif', 0));
            testCase.verifyTrue(evalin('base', 'exist(''pkRunnerKeep'', ''var'')') == 1, 'A base variable was cleared.');
            testCase.verifySubstring(fileread(st), 'x none');
        end

        function matlabRunnerRecordsTaskMetrics(testCase)
            % with task metrics, each task a MATLAB runner runs has its
            % record; the job script starts MATLAB with the launch string
            testCase.assumeLinux();
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(3);
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''copy'')', x, y, logFn), ...
                inFns, outFns, 'unif', 0);
            runDir = [dataPath, 'metrics/run_1/'];
            mkdir(runDir);
            launch = 'module load matlab/r2024b; matlab -nodisplay -nosplash -nodesktop -nojvm -batch';
            evalc(['rc = task_runner_run(''setup_matlab'', [dataPath, ''call''], funcStrs, outFns, '''', ', ...
                'false(3, 1), 3, launch, runDir, false);']);
            testCase.verifySubstring(fileread(rc.script), sprintf('MLAUNCH=''%s''', launch));
            L = [testCase.createTempFolder(), 'local'];
            mkdir(L);
            system(sprintf('ln -s "0_0:test:1:1" "%s/t1"', rc.claimsDir));
            testCase.writeLines([L, '/s1_1'], {'c 1'});
            setenv('PETAKIT5D_RUNNER_MEM_GB', '100');
            evalc('task_runner(''matlab'', rc.D, ''1'', 0, [L, ''/s1_1''], L, 1)');
            t = read_task_metrics(runDir);
            testCase.assertEqual(height(t), 3);
            testCase.verifyEqual(sort(t.task)', 1 : 3);
            testCase.verifyEqual(t.mode, repmat("slurm_matlab", 3, 1));
            testCase.verifyEqual(t.status, repmat("ok", 3, 1));
            testCase.verifyEqual(numel(unique(t.pid)), 1);
        end
    end

    methods (Test, TestTags = {'Unit'})
        function compiledJobsRunTheJobScript(testCase)
            % with a build that has task runners, a first submission runs the
            % batch with the call's job script (one slot per GNU parallel
            % slot, stealing, two tries as --retries 2), with the job's
            % resources as before; a sequential job has one slot and one try;
            % a GPU job gives each slot its GPU
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [dataPath, cmds] = testCase.compiledCall(slurm, 8, {}); %#ok<ASGLU>
            testCase.assertNumElements(cmds, 2);
            for b = 1 : 2
                testCase.verifyNotEmpty(regexp(cmds{b}, sprintf(['--cpus-per-task=4 --ntasks=1 .*echo; bash \\S+/run_tasks\\.sh ', ...
                    'job %d 4 0 1 2; echo; printf'], b), 'once'), cmds{b});
            end
            testCase.verifyEmpty(regexp([cmds{:}], 'parallel ', 'once'), 'A job runs the task lines with GNU parallel.');

            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [~, cmds] = testCase.compiledCall(slurm, 2, {'GNUparallel', false});
            testCase.verifyNotEmpty(regexp(cmds{1}, 'run_tasks\.sh job 1 1 0 1 1;', 'once'), cmds{1});
            testCase.verifyEmpty(regexp(cmds{1}, '--ntasks=1', 'once'), 'The sequential job has other resources.');

            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [~, cmds] = testCase.compiledCall(slurm, 2, {'GPUJob', true, 'GPUNum', 2, 'maxCPUNum', 8});
            testCase.verifyNotEmpty(regexp(cmds{1}, 'run_tasks\.sh job 1 2 1 1 2;', 'once'), cmds{1});
        end

        function taskRunnerOffOrOldBuildsRunAsBefore(testCase)
            % taskRunner false (option or config key), or a build without
            % task runners: the job commands of before and no runner files
            testCase.assumeLinux();
            expected = 'echo; parallel --ungroup --retries 2 --jobs 4 --delay 0.2 < ';
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [dataPath, cmds] = testCase.compiledCall(slurm, 4, {}, 'taskRunner', false);
            testCase.verifySubstring(cmds{1}, expected);
            testCase.verifyEmpty(dir([dataPath, 'tmp/*/run_tasks.sh']));
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [dataPath, cmds] = testCase.compiledCall(slurm, 4, {'taskRunner', false});
            testCase.verifySubstring(cmds{1}, expected);
            testCase.verifyEmpty(dir([dataPath, 'tmp/*/run_tasks.sh']));
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [dataPath, cmds] = testCase.compiledCall(slurm, 4, {}, 'oldBuild', true);
            testCase.verifySubstring(cmds{1}, expected);
            testCase.verifyEmpty(dir([dataPath, 'tmp/*/run_tasks.sh']));
        end

        function compiledEmptyJobLogDirIsTheDefaultFolder(testCase)
            % compiled mode: an empty jobLogDir (as from a config file) gave
            % the jobs "-o  -e", so their logs went to a file named -e in the
            % folder the framework submits from (the code folder); it is now
            % the default folder (../job_logs beside the code folder if it
            % is there, else <data folder>/job_logs; the fake jobs make
            % their log files only in the test's folders)
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.makeLogFiles();
            dataPath = testCase.createTempFolder();
            b = testCase.recordingBuild([dataPath, 'build/']); %#ok<NASGU>
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : 2, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout_%d.txt', dataPath, x), 1 : 2, 'unif', 0);
            cellfun(@(x) fclose(fopen(x, 'w')), inFns);
            funcStrs = cellfun(@(x, y) sprintf('XR_crop_frame(''%s'',''%s'',[1,1,1,2,2,2])', x, y), inFns, outFns, 'unif', 0); %#ok<NASGU>
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(dataPath));
            evalc(['mcc_slurm_cluster_generic_computing_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''taskBatchNum'', 1, ''maxTrialNum'', 1, ''taskRunner'', false, ', ...
                '''MCCMasterStr'', b.master, ''MCRParam'', b.mcr, ''MCRCacheRoot'', [dataPath, ''cache''], ', ...
                '''BashLaunchStr'', ''echo'', ''jobLogDir'', '''')']);
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

        function matlabJobsRunTheJobScript(testCase)
            % MATLAB jobs: the job script starts MATLAB (with the launch string
            % of before) if tasks are left; off: the command of before
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''noop'')', x, y, logFn), ...
                inFns, outFns, 'unif', 0);
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''jobLogDir'', dataPath)']);
            cmds = slurm.calls('sbatch');
            testCase.assertNumElements(cmds, 2);
            for b = 1 : 2
                testCase.verifyNotEmpty(regexp(cmds{b}, sprintf('; bash \\S+/run_tasks\\.sh matlab %d 1; echo;', b), 'once'), cmds{b});
                testCase.verifyEmpty(regexp(cmds{b}, '-batch ', 'once'), 'The job starts MATLAB itself.');
            end
            % the call's files are removed at the end
            testCase.verifyEmpty(dir([dataPath, 'task_runner_*']));

            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''jobLogDir'', dataPath, ''taskRunner'', false)']);
            cmds = slurm.calls('sbatch');
            testCase.verifySubstring(cmds{1}, sprintf('-nojvm -batch "setup([]);t0_=tic;%s;toc(t0_)"', funcStrs{1}));
        end

        function compiledResubmissionsRunAsBefore(testCase)
            % compiled mode: a batch resubmitted after its job ended without
            % its output runs the job command of before
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            b = testCase.recordingBuild([dataPath, 'build/']);
            configFn = [dataPath, 'config.json'];
            testCase.writeLines(configFn, {jsonencode(struct('MCRParam', b.mcr, 'MCCMasterStr', b.master, ...
                'MCRCacheRoot', [dataPath, 'cache'], 'BashLaunchStr', 'echo', 'GNUparallel', true, ...
                'wholeNodeJob', true, 'maxCPUNum', 4))});
            inFns = {[dataPath, 'in_1.txt']};
            outFns = {[dataPath, 'out_1.txt']};
            fclose(fopen(inFns{1}, 'w'));
            funcStrs = {sprintf('XR_crop_frame(''%s'',''%s'',[1,1,1,2,2,2])', inFns{1}, outFns{1})}; %#ok<NASGU>
            slurm.setJobStates('1000_1', {'R', 'done'});
            slurm.setJobStates('1001_1', {'R', 'done'}, outFns);
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''mccMode'', true, ''maxTrialNum'', 2, ''queryInterval'', 1, ', ...
                '''jobLogDir'', dataPath, ''ConfigFile'', configFn)']);
            cmds = slurm.calls('sbatch');
            testCase.assertNumElements(cmds, 2);
            testCase.verifySubstring(cmds{1}, 'run_tasks.sh job 1 ');
            testCase.verifyEmpty(strfind(cmds{2}, 'run_tasks.sh'), cmds{2});
        end

        function matlabResubmissionsRunAsBefore(testCase)
            % a batch resubmitted after its job ended without its output
            % runs the MATLAB command of before (no runner, no stealing)
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(1);
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''noop'')', x, y, logFn), ...
                inFns, outFns, 'unif', 0);
            slurm.setJobStates('1000_1', {'R', 'done'});
            slurm.setJobStates('1001_1', {'R', 'done'}, outFns(1));
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 2, ''queryInterval'', 1, ''jobLogDir'', dataPath)']);
            cmds = slurm.calls('sbatch');
            testCase.assertNumElements(cmds, 2);
            testCase.verifySubstring(cmds{1}, 'run_tasks.sh matlab 1 1;');
            testCase.verifyEmpty(strfind(cmds{2}, 'run_tasks.sh'), cmds{2});
            testCase.verifySubstring(cmds{2}, '-nojvm -batch');
        end

        function masterLeavesTasksToRunningRunnersAndResubmitsStaleOnes(testCase)
            % a task another job's runner took is not submitted while that
            % job runs; when the job is gone without its output, the claim
            % is removed and the task's batch submitted
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(3);
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''noop'')', x, y, logFn), ...
                inFns, outFns, 'unif', 0);
            % jobs: batch 1 (1000) runs for 6 polls with its output, batch 2
            % (1001) is done with its output, batch 3 (1002) once submitted
            slurm.setJobStates('1000_1', [repmat({'R'}, 1, 6), {'done'}], outFns(1));
            slurm.setJobStates('1001_2', {'R', 'done'}, outFns(2));
            slurm.setJobStates('1002_3', {'R'}, outFns(3));
            % job 1000's runner takes task 3 before the master submits it
            claimer = timer('StartDelay', 0.5, 'TimerFcn', @(~, ~) testCase.claimInCall(dataPath, 3, '1000_1'));
            testCase.addTeardown(@() delete(claimer));
            start(claimer);

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 2, ''maxJobNum'', 2, ''queryInterval'', 1, ''jobLogDir'', dataPath)']);

            cmds = slurm.history();
            submit3 = find(startsWith(cmds, 'sbatch') & contains(cmds, '--array=3 '));
            testCase.assertNumElements(submit3, 1, 'Batch 3 should be submitted once.');
            running = find(startsWith(cmds, 'squeue') & contains(cmds, '1000_1'));
            testCase.assertGreaterThanOrEqual(numel(running), 6);
            testCase.verifyGreaterThan(submit3, running(6), 'Batch 3 was submitted while the job that took its task ran.');
            testCase.verifyTrue(all(cellfun(@isfile, outFns)));
        end

        function claimStatesFollowTheirOwners(testCase)
            % a claim is live while its owner's job is in the queue (running
            % or pending), or (no job) while its process runs on this node;
            % failed claims, the master's own and those of ended jobs are not
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('77_1', {'R'});
            slurm.setJobStates('79_1', {'PD'});
            [~, outFns] = testCase.makeTasks(7);
            evalc('rc = task_runner_run(''setup_matlab'', [testCase.createTempFolder(), ''call''], outFns, outFns, '''', false(7, 1), 1, ''x'', '''', false);');
            % the node name as the job script writes it
            node = getenv('SLURMD_NODENAME');
            if isempty(node)
                [~, node] = system('uname -n');
                node = strtrim(node);
            end
            [~, pid] = system('sleep 30 > /dev/null & echo $!');
            testCase.addTeardown(@() system(['kill ', strtrim(pid)]));
            owners = {'77_1:n1:5:1', '78_1:n1:5:1', sprintf('0_0:%s:%s:0', node, strtrim(pid)), ...
                sprintf('0_0:%s:%d:0', node, feature('getpid')), '77_1:n1:5:1', '79_1:n2:6:2'};
            for f = 1 : 6
                system(sprintf('ln -s "%s" "%s/t%d"', owners{f}, rc.claimsDir, f));
            end
            system(sprintf('ln -s "77_1:n1:5:1" "%s/f5"', rc.claimsDir));
            [claimed, live] = task_runner_run('claims', rc, 1 : 7);
            testCase.verifyEqual(claimed', [true, true, true, true, true, true, false]);
            testCase.verifyEqual(live', [true, false, true, false, false, true, false]);
            % release removes the claims and failed marks, and marks a batch resubmitted
            task_runner_run('release', rc, [2, 5], 3);
            names = dir(rc.claimsDir);
            names = setdiff({names.name}, {'.', '..'});
            testCase.verifyEqual(sort(names), sort({'t1', 't3', 't4', 't6', 'r3'}));
        end

        function oneTaskPerSlotRunsTaskLinesAsBefore(testCase)
            % with one task per slot and no batch waiting in the queue, the
            % task lines run as before (no runner; a slot that is done takes
            % another batch's task the same way); with a batch that waits,
            % the processes are runners and take its tasks
            testCase.assumeLinux();
            testCase.assumeEqual(system('command -v parallel > /dev/null'), 0, 'GNU parallel is not installed.');
            slow = [num2cell((1 : 4)'), repmat({'slow'}, 4, 1)];
            [c, b] = testCase.makeCall(8, 4, slow);
            out = testCase.runJob(c, 'job 1 4 0 1 2');
            runs = testCase.taskRuns(b);
            testCase.verifyEqual(sort([runs.task]), 1 : 8);
            testCase.verifyEqual([runs.runner], zeros(1, 8), 'A task line ran with a runner.');
            testCase.verifyEqual(numel(unique([runs.pid])), 8, 'Tasks shared a process.');
            testCase.verifyEmpty(strfind(out, 'did not finish'), out);
            [c, b] = testCase.makeCall(8, 4, slow);
            testCase.writeLines([c.rc.D, '/waiting.txt'], {'2'});
            testCase.runJob(c, 'job 1 4 0 1 2');
            runs = testCase.taskRuns(b);
            testCase.verifyEqual(sort([runs.task]), 1 : 8);
            testCase.verifyTrue(all([runs.runner] == 1));
        end

        function compiledRunnerTasksGetTheEnvironmentOfTheirLines(testCase)
            % the metrics environment a task line sets (PETAKIT5D_TASK...)
            % is set for the task in the runner too
            testCase.assumeLinux();
            dataPath = testCase.createTempFolder();
            b = testCase.recordingBuild([dataPath, 'build/']);
            lines = cell(3, 1);
            for f = 1 : 3
                lines{f} = sprintf(['OUT="%sout_%d"; if [ ! -f "$OUT" ] && [ ! -d "$OUT" ]; then ', ...
                    'PETAKIT5D_TASK_METRICS="%sm/" PETAKIT5D_TASK=%d MCR_CACHE_ROOT=/tmp/x %s "%s" pktest.runnerEnvTask %sout_%d ; ', ...
                    'else echo output "$OUT" already exists; fi'], dataPath, f, dataPath, f, b.master, b.mcr, dataPath, f);
            end
            evalc('rc = task_runner_run(''setup'', dataPath, lines, b.master, 3, [], false(3, 1));');
            testCase.assertTrue(rc.on);
            testCase.runInProcess([dataPath, 'tasks.txt'], b.master, 1 : 3, 'env');
            for f = 1 : 3
                testCase.verifyEqual(strtrim(fileread(sprintf('%sout_%d', dataPath, f))), sprintf('%sm/|%d', dataPath, f));
            end
        end

        function compiledMasterLeavesTasksToRunningRunners(testCase)
            % compiled mode: a task another job's runner took is not
            % submitted while that job runs; then its batch is submitted
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            b = testCase.recordingBuild([dataPath, 'build/']);
            configFn = [dataPath, 'config.json'];
            testCase.writeLines(configFn, {jsonencode(struct('MCRParam', b.mcr, 'MCCMasterStr', b.master, ...
                'MCRCacheRoot', [dataPath, 'cache'], 'BashLaunchStr', 'echo', 'maxJobNum', 2))});
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : 3, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout_%d.txt', dataPath, x), 1 : 3, 'unif', 0);
            cellfun(@(x) fclose(fopen(x, 'w')), inFns);
            funcStrs = cellfun(@(x, y) sprintf('XR_crop_frame(''%s'',''%s'',[1,1,1,2,2,2])', x, y), inFns, outFns, 'unif', 0);
            slurm.setJobStates('1000_1', [repmat({'R'}, 1, 6), {'done'}], outFns(1));
            slurm.setJobStates('1001_2', {'R', 'done'}, outFns(2));
            slurm.setJobStates('1002_3', {'R'}, outFns(3));
            claimer = timer('StartDelay', 0.5, 'TimerFcn', @(~, ~) testCase.claimInCompiledCall(dataPath, 3, '1000_1'));
            testCase.addTeardown(@() delete(claimer));
            start(claimer);

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''mccMode'', true, ''maxTrialNum'', 2, ''queryInterval'', 1, ', ...
                '''jobLogDir'', dataPath, ''ConfigFile'', configFn)']);

            cmds = slurm.history();
            submit3 = find(startsWith(cmds, 'sbatch') & contains(cmds, '--array=3 '));
            testCase.assertNumElements(submit3, 1, 'Batch 3 should be submitted once.');
            running = find(startsWith(cmds, 'squeue') & contains(cmds, '1000_1'));
            testCase.assertGreaterThanOrEqual(numel(running), 6);
            testCase.verifyGreaterThan(submit3, running(6), 'Batch 3 was submitted while the job that took its task ran.');
            testCase.verifyTrue(all(cellfun(@isfile, outFns)));
        end

        function compiledMasterKeepsItsBatch(testCase)
            % compiled mode: the batch the master runs itself (the last one)
            % has its own key, so other jobs' runners leave its tasks to it;
            % the other batches have the CPUs of their jobs (the compiled
            % wrapper directly: the framework lets only a compiled master
            % compute in compiled mode)
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            b = testCase.recordingBuild([dataPath, 'build/']); %#ok<NASGU>
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : 6, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout_%d.txt', dataPath, x), 1 : 6, 'unif', 0);
            cellfun(@(x) fclose(fopen(x, 'w')), inFns);
            funcStrs = cellfun(@(x, y) sprintf('XR_crop_frame(''%s'',''%s'',[1,1,1,2,2,2])', x, y), inFns, outFns, 'unif', 0); %#ok<NASGU>
            evalc(['mcc_slurm_cluster_generic_computing_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', true, ''taskBatchNum'', 2, ''cpusPerTask'', 3, ''maxCPUNum'', 3, ''maxTrialNum'', 1, ', ...
                '''MCCMasterStr'', b.master, ''MCRParam'', b.mcr, ''MCRCacheRoot'', [dataPath, ''cache''], ', ...
                '''BashLaunchStr'', ''echo'', ''jobLogDir'', dataPath, ''tmpDir'', [dataPath, ''tmp/''])']);
            keyFiles = dir([dataPath, 'tmp/*/keys.txt']);
            testCase.assertNotEmpty(keyFiles, 'The call has no keys file.');
            keys = splitlines(strip(fileread(fullfile(keyFiles(1).folder, keyFiles(1).name))));
            testCase.verifyEqual(keys, {'3'; '3'; 'master'});
            cmds = slurm.calls('sbatch');
            testCase.verifyNotEmpty(regexp([cmds{:}], '--cpus-per-task=3 ', 'once'));
        end

        function compiledLocalCallDoneInItsLastTryRemovesItsTaskFolder(testCase)
            % compiled mode, a local run (as a compiled stitching job runs
            % its cross-correlations): the batch run in its last try ran
            % after its round looked at the outputs, so the call returned it
            % unfinished, kept its task folder in tmp/ (the task lines and the
            % runners' files) and its caller ran it again; the call now looks
            % once more at its end. With task runners (their pass runs every
            % task once) and without
            testCase.assumeLinux();
            for runner = [true, false]
                dataPath = testCase.createTempFolder();
                b = testCase.fakeBuild([dataPath, 'build/']);
                inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : 4, 'unif', 0);
                outFns = arrayfun(@(x) sprintf('%sout_%d.txt', dataPath, x), 1 : 4, 'unif', 0);
                cellfun(@(x) fclose(fopen(x, 'w')), inFns);
                funcStrs = cellfun(@(x) sprintf('pk_task(''%s'',''second'')', x), outFns, 'unif', 0); %#ok<NASGU>
                flags = []; %#ok<NASGU>
                evalc(['flags = mcc_slurm_cluster_generic_computing_wrapper(inFns, outFns, funcStrs, ''parseCluster'', false, ', ...
                    '''masterCompute'', true, ''taskBatchNum'', 2, ''maxTrialNum'', 2, ''unitWaitTime'', 0, ', ...
                    '''taskRunner'', runner, ''MCCMasterStr'', b.master, ''MCRParam'', b.mcr, ', ...
                    '''MCRCacheRoot'', [dataPath, ''cache''], ''BashLaunchStr'', ''echo'', ''jobLogDir'', dataPath)']);
                testCase.verifyEqual(flags(:), true(4, 1), sprintf('runners %d: a batch done in its last try counts as unfinished.', runner));
                testCase.verifyTrue(all(cellfun(@isfile, outFns)));
                left = dir([dataPath, 'tmp/2*']);
                testCase.verifyEmpty(left, sprintf('runners %d: the task folder %s is left.', runner, strjoin({left.name}, ', ')));
            end
        end

        function matlabMasterKeepsItsBatch(testCase)
            % MATLAB jobs: the same
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(3);
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''copy'')', x, y, logFn), ...
                inFns, outFns, 'unif', 0);
            slurm.setJobStates('1000_1', {'R', 'R', 'R', 'done'}, outFns(1));
            slurm.setJobStates('1001_2', {'R', 'R', 'R', 'done'}, outFns(2));
            % (the call's folder is removed at its end)
            copier = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.2, 'TasksToExecute', 50, ...
                'TimerFcn', @(~, ~) testCase.copyCallFile(dataPath, 'keys.txt'));
            testCase.addTeardown(@() delete(copier));
            start(copier);
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', true, ''maxTrialNum'', 2, ''queryInterval'', 1, ''jobLogDir'', dataPath)']);
            stop(copier);
            testCase.assertTrue(isfile([dataPath, 'keys.txt']), 'The call has no keys file.');
            keys = splitlines(strip(fileread([dataPath, 'keys.txt'])));
            testCase.assertNumElements(keys, 3);
            testCase.verifyEqual(keys{3}, 'master');
            testCase.verifyEqual(keys{1}, keys{2});
            testCase.verifyNotEqual(keys{1}, 'master');
        end
    end

    methods (Access = private)
        function assumeLinux(testCase)
            testCase.assumeTrue(isunix && ~ismac, 'Task runners run on Linux only.');
        end

        function [c, b] = makeCall(testCase, nT, BS, modes, kind, excluded, callDir)
            % a call folder (task_runner_run setup) for nT fake tasks in
            % batches of BS; modes: {task, mode; ...} ('error', 'crash',
            % 'slow', 'giveback': a runner gives it back once; else 'ok');
            % kind 'local' for a local run; excluded: tasks without inputs;
            % callDir: the call folder (default: call/ in a new folder)
            if nargin < 5
                kind = '';
            end
            if nargin < 6
                excluded = false(nT, 1);
            end
            dataPath = testCase.createTempFolder();
            if nargin < 7
                callDir = [dataPath, 'call/'];
            end
            b = testCase.fakeBuild([dataPath, 'build/']);
            mkdir([dataPath, 'out']);
            c.outFns = arrayfun(@(x) sprintf('%sout/o_%d.txt', dataPath, x), 1 : nT, 'unif', 0);
            lines = cell(nT, 1);
            for f = 1 : nT
                mode = 'ok';
                if ~isempty(modes)
                    k = find([modes{:, 1}] == f, 1);
                    if ~isempty(k)
                        mode = modes{k, 2};
                    end
                end
                lines{f} = sprintf(['OUT="%s"; if [ ! -f "$OUT" ] && [ ! -d "$OUT" ]; then PETAKIT5D_TASK=%d ', ...
                    'MCR_CACHE_ROOT=/tmp/mcr_cache_${SLURM_JOB_ID:-0}_${PARALLEL_JOBSLOT:-$$} %s "%s" pk_task %s %s ; ', ...
                    'else echo output "$OUT" already exists; fi'], c.outFns{f}, f, b.master, b.mcr, c.outFns{f}, mode);
            end
            tmpInfo = [];
            c.tmpFns = strrep(c.outFns, '.txt', '.tmp');
            if strcmp(kind, 'local')
                tmpInfo = struct('tmpFns', {c.tmpFns}, 'outFns', {c.outFns}, 'unitWaitTime', 3);
            end
            evalc('c.rc = task_runner_run(''setup'', callDir, lines, b.master, BS, tmpInfo, excluded);');
            testCase.assertTrue(c.rc.on, 'The runner files could not be written.');
        end

        function [b] = fakeBuild(testCase, folder)
            % a stand-in compiled master: run_mccMaster.sh as the compiler
            % makes it (its argument handling) and a bash mccMaster that runs
            % 'pk_task <output> <mode>' and, as a runner, the next tasks
            % (mode 'second': the output only in the task's second run)
            mkdir(folder);
            b.master = [folder, 'run_mccMaster.sh'];
            b.mcr = [folder, 'mcr'];
            b.log = [folder, 'runs.log'];
            mkdir([b.mcr, '/bin']);
            testCase.writeLines(b.master, testCase.runMccMasterScript());
            testCase.writeLines([folder, 'mccMaster'], { ...
                '#!/bin/bash'
                'CFG="$PETAKIT5D_RUNNER"'
                'R=0; [ -n "$CFG" ] && R=1'
                'unset PETAKIT5D_RUNNER'
                'run_one() {'
                sprintf('    echo "$2 $3 ${PETAKIT5D_TASK:-} $$ ${CUDA_VISIBLE_DEVICES:-none} $R" >> "%s"', b.log)
                '    case "$3" in'
                '        error) return 1 ;;'
                '        crash) kill -9 $$ ;;'
                '        slow) sleep 1; echo x > "$2" ;;'
                '        second) if [ -f "$2.tried" ]; then echo x > "$2"; else : > "$2.tried"; fi ;;'
                '        *) echo x > "$2" ;;'
                '    esac'
                '}'
                'run_one "$@" || exit 1'
                '[ -z "$CFG" ] && exit 0'
                'IFS="|" read -r D own steal s st L task <<< "$CFG"'
                'printf "d %s\nt 0.1\n" "$task" >> "$st"'
                'while true; do'
                '    bash "$D/run_tasks.sh" claim "$own" "$steal" "$s" 1 "$st" "$L" "$L/next_$s"'
                '    f=$(cat "$L/next_$s")'
                '    [ -z "$f" ] && { echo "x none" >> "$st"; exit 0; }'
                '    echo "c $f" >> "$st"'
                '    mapfile -d "" parts < "$L/dump_$s"'
                '    if [ "${parts[5]}" = giveback ] && [ ! -f "$L/gaveback_$f" ]; then'
                '        touch "$L/gaveback_$f"; rm -f "$D/claims/t$f"; printf "s %s\nx recycle\n" "$f" >> "$st"; exit 0'
                '    fi'
                '    if PETAKIT5D_TASK="${parts[1]}" run_one "${parts[@]:3}"; then'
                '        printf "d %s\nt 0.1\n" "$f" >> "$st"'
                '    else'
                '        printf "e %s\nx failed\n" "$f" >> "$st"'
                '        exit 0'
                '    fi'
                'done'});
            fileattrib(b.master, '+x');
            fileattrib([folder, 'mccMaster'], '+x');
            testCase.writeLines([folder, 'task_runner.txt'], {sprintf('%d', task_runner('version'))});
        end

        function [lines] = runMccMasterScript(~)
            % run_mccMaster.sh as compile_mccMaster leaves it (Linux)
            lines = {
                '#!/bin/sh'
                'exe_name=$0'
                'exe_dir=`dirname "$0"`'
                'if [ "x$1" = "x" ]; then'
                '  echo Usage:'
                'else'
                '  MCRROOT="$1"'
                '  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;'
                '  export LD_LIBRARY_PATH;'
                '  shift 1'
                '  args='
                '  while [ $# -gt 0 ]; do'
                '      token=$1'
                '      args="${args} \"${token}\"" '
                '      shift'
                '  done'
                '  eval "\"${exe_dir}/mccMaster\"" $args'
                'fi'
                'exit'};
        end

        function [b] = argsBuild(testCase, folder)
            % a build whose mccMaster writes the arguments it gets (PK_ARGS_OUT)
            mkdir(folder);
            b.master = [folder, 'run_mccMaster.sh'];
            testCase.writeLines(b.master, testCase.runMccMasterScript());
            testCase.writeLines([folder, 'mccMaster'], {'#!/bin/bash', 'printf ''%s\0'' "$@" > "$PK_ARGS_OUT"'});
            fileattrib(b.master, '+x');
            fileattrib([folder, 'mccMaster'], '+x');
        end

        function [out] = runJob(testCase, c, args)
            % runs the job script (bash run_tasks.sh <args>) as a job would
            [status, out] = system(sprintf('env -u PETAKIT5D_RUNNER bash "%s" %s 2>&1', c.rc.script, args));
            testCase.assertEqual(status, 0, out);
        end

        function [runs] = taskRuns(~, b)
            % the fake tasks that ran: task, output, mode, process, device
            runs = struct('task', {}, 'out', {}, 'mode', {}, 'pid', {}, 'cuda', {}, 'runner', {});
            if ~isfile(b.log)
                return;
            end
            c = textscan(fileread(b.log), '%s %s %f %f %s %f');
            for i = 1 : numel(c{1})
                runs(i) = struct('task', c{3}(i), 'out', c{1}{i}, 'mode', c{2}{i}, 'pid', c{4}(i), 'cuda', c{5}{i}, ...
                    'runner', c{6}(i));
            end
        end

        function [names] = claimNames(~, c)
            dir_info = dir(c.rc.claimsDir);
            names = setdiff({dir_info.name}, {'.', '..'});
        end

        function [lines, taskFn, master] = frameworkTaskLines(testCase, dataPath, kind)
            % task lines the framework writes in compiled mode (fake Slurm;
            % the jobs never run) for tasks with mixed arguments, through a
            % build with task runners; returns them, the call's task file and
            % the build's run_mccMaster.sh
            build = testCase.recordingBuild([dataPath, 'build/']);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture); %#ok<NASGU>
            configFn = [dataPath, 'config.json'];
            testCase.writeLines(configFn, {jsonencode(struct('MCRParam', build.mcr, 'MCCMasterStr', build.master, ...
                'MCRCacheRoot', [dataPath, 'cache'], 'BashLaunchStr', 'echo', 'GNUparallel', true, ...
                'wholeNodeJob', true, 'maxCPUNum', 4))});
            inPath = [dataPath, 'in dir/'];
            outPath = [dataPath, 'in dir/out/'];
            mkdir(outPath);
            nF = 4;
            inFns = arrayfun(@(x) sprintf('%sim_%d.tif', inPath, x), (1 : nF)', 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sim_%d.zarr', outPath, x), (1 : nF)', 'unif', 0);
            for f = 1 : nF
                testCase.writeLines(inFns{f}, {'x'});
            end
            funcStrs = cell(nF, 1);
            for f = 1 : nF
                switch kind
                    case 'mixed'
                        % cell arrays, empty strings and numbers, as the wrappers write
                        % them (and a non-ASCII character: lines are read by byte offset)
                        funcStrs{f} = sprintf(['tiffToZarr(%s,''%s'',[],''BlockSize'',%s,''flipZstack'',%s,', ...
                            '''resampleFactor'',%s,''compressor'',''%s'',''usrFcn'',"%s",''uuid'',''a b', char(233), ''')'], ...
                            sprintf('{''%s''}', inFns{f}), outFns{f}, '[32,32,16]', 'false', '[]', 'zstd', '');
                end
            end
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''mccMode'', true, ''maxTrialNum'', 1, ''jobLogDir'', dataPath, ', ...
                '''tmpDir'', [dataPath, ''tmp/''], ''ConfigFile'', configFn, ''finalOutFullpath'', [dataPath, ''final.txt''])']);
            taskFiles = dir([dataPath, 'tmp/*/tasks.txt']);
            testCase.assertNotEmpty(taskFiles, 'No call folder with task runner files.');
            taskFn = fullfile(taskFiles(1).folder, taskFiles(1).name);
            lines = splitlines(strip(fileread(taskFn)));
            master = build.master;
        end

        function [b] = recordingBuild(testCase, folder)
            % a build with task runners whose master does nothing (the
            % framework only writes the task lines here)
            mkdir([folder, 'mcr/bin']);
            b.master = [folder, 'run_mccMaster.sh'];
            b.mcr = [folder, 'mcr'];
            testCase.writeLines(b.master, testCase.runMccMasterScript());
            testCase.writeLines([folder, 'mccMaster'], {'#!/bin/bash', 'exit 0'});
            fileattrib(b.master, '+x');
            fileattrib([folder, 'mccMaster'], '+x');
            testCase.writeLines([folder, 'task_runner.txt'], {sprintf('%d', task_runner('version'))});
        end

        function [master] = stateTaskLines(testCase, dataPath, outFns, modes)
            % the call folder of task lines that run pktest.runnerStateTask;
            % returns the build's run_mccMaster.sh
            b = testCase.recordingBuild([dataPath, 'build/']);
            n = numel(outFns);
            lines = cell(n, 1);
            for f = 1 : n
                lines{f} = sprintf(['OUT="%s"; if [ ! -f "$OUT" ] && [ ! -d "$OUT" ]; then MCR_CACHE_ROOT=/tmp/x %s "%s" ', ...
                    'pktest.runnerStateTask %s %s ; else echo output "$OUT" already exists; fi'], outFns{f}, b.master, ...
                    b.mcr, outFns{f}, modes{f});
            end
            evalc('rc = task_runner_run(''setup'', dataPath, lines, b.master, n, [], false(n, 1));');
            testCase.assertTrue(rc.on);
            master = b.master;
        end

        function [recorded, status] = runInProcess(testCase, taskFn, master, tasks, kind)
            % the runner in this MATLAB as the job script starts it for the
            % first of tasks (task line run with PETAKIT5D_RUNNER): the first
            % task runs as mccMaster would run it, the next ones through
            % task_runner; mccMaster is replaced by a recorder of its
            % arguments ('state': which also runs the task function)
            if nargin < 5
                kind = 'record';
            end
            if isempty(getenv('PETAKIT5D_RUNNER_MEM_GB'))
                % (this MATLAB's memory varies with the other tests)
                setenv('PETAKIT5D_RUNNER_MEM_GB', '100');
            end
            D = fileparts(taskFn);
            recDir = [tempname, '/'];
            mkdir(recDir);
            testCase.addTeardown(@() rmdir(recDir, 's'));
            body = {'function mccMaster(varargin)', ...
                'n = numel(dir([getenv(''PK_REC_DIR''), ''rec_*.mat'']));', ...
                'save(sprintf(''%srec_%03d.mat'', getenv(''PK_REC_DIR''), n + 1), ''varargin'');'};
            if any(strcmp(kind, {'state', 'env'}))
                body{end + 1} = 'feval(varargin{1}, varargin{2 : end});';
            end
            body{end + 1} = 'end';
            testCase.writeLines([recDir, 'mccMaster.m'], body);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(recDir));
            setenv('PK_REC_DIR', recDir);
            testCase.addTeardown(@setenv, 'PK_REC_DIR', '');
            L = [tempname, '/'];
            mkdir(L);
            testCase.addTeardown(@() rmdir(L, 's'));
            st = [L, 's1_1'];
            f = tasks(1);
            system(sprintf('ln -s "0_0:test:1:1" "%s/claims/t%d"', D, f));
            testCase.writeLines(st, {sprintf('c %d', f)});
            % the first task's arguments from its line, as its process gets them
            line = fileread(taskFn);
            line = splitlines(line);
            line = strrep(line{f}, master, [D, '/argdump/run_mccMaster.sh']);
            testCase.writeLines([L, 'line.sh'], {line});
            [status, out] = system(sprintf('PK_DUMP="%sdump0" bash "%sline.sh"', L, L));
            testCase.assertEqual(status, 0, out);
            args = testCase.readNul([L, 'dump0']);
            % the first task's process has the environment of its line
            names = {'PETAKIT5D_TASK_METRICS', 'PETAKIT5D_TASK', 'PETAKIT5D_TASK_GPU'};
            for i = 1 : 3
                if ~isempty(args{i})
                    setenv(names{i}, args{i});
                end
            end
            % the compiled runner's reset clears the base workspace (a compiled
            % master has nothing there): this MATLAB's comes back after the test
            baseNames = evalin('base', 'who');
            baseValues = cellfun(@(x) evalin('base', x), baseNames, 'unif', 0);
            testCase.addTeardown(@() cellfun(@(x, y) assignin('base', x, y), baseNames, baseValues));
            setenv('PETAKIT5D_RUNNER', sprintf('%s|1|0|1|%s|%s|%d', D, st, L(1 : end - 1), f));
            r = task_runner('start');
            testCase.assertNotEmpty(r);
            mccMaster(args{4 : end});
            evalc('task_runner(''continue'', r, 0.1)');
            recs = dir([recDir, 'rec_*.mat']);
            [~, order] = sort({recs.name});
            recorded = arrayfun(@(x) getfield(load(fullfile(x.folder, x.name), 'varargin'), 'varargin'), ...
                recs(order), 'unif', 0); %#ok<GFLD>
            status = fileread(st);
        end

        function [dataPath, cmds] = compiledCall(testCase, slurm, nT, config, varargin)
            % the framework in compiled mode (whole-node jobs of 4 slots, GNU
            % parallel) on nT fake tasks with a build that has task runners
            % ('oldBuild' true: without); the jobs never run. Returns the
            % sbatch commands
            ip = inputParser;
            ip.addParameter('taskRunner', true);
            ip.addParameter('oldBuild', false);
            ip.parse(varargin{:});
            dataPath = testCase.createTempFolder();
            b = testCase.recordingBuild([dataPath, 'build/']);
            if ip.Results.oldBuild
                delete([dataPath, 'build/task_runner.txt']);
            end
            s = struct('MCRParam', b.mcr, 'MCCMasterStr', b.master, 'MCRCacheRoot', [dataPath, 'cache'], ...
                'BashLaunchStr', 'echo', 'GNUparallel', true, 'wholeNodeJob', true, 'maxCPUNum', 4);
            for i = 1 : 2 : numel(config)
                s.(config{i}) = config{i + 1};
            end
            configFn = [dataPath, 'config.json'];
            testCase.writeLines(configFn, {jsonencode(s)});
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : nT, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout_%d.txt', dataPath, x), 1 : nT, 'unif', 0);
            cellfun(@(x) fclose(fopen(x, 'w')), inFns);
            funcStrs = cellfun(@(x, y) sprintf('XR_crop_frame(''%s'',''%s'',[1,1,1,2,2,2])', x, y), inFns, outFns, 'unif', 0);
            args = {'parseCluster', true, 'masterCompute', false, 'mccMode', true, 'maxTrialNum', 1, ...
                'jobLogDir', dataPath, 'ConfigFile', configFn};
            if ~ip.Results.taskRunner
                args = [args, {'taskRunner', false}];
            end
            evalc('generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, args{:})');
            cmds = slurm.calls('sbatch');
        end

        function claimInCompiledCall(~, dataPath, f, job)
            % claims task f of the (only) compiled call in dataPath for job
            d = dir([dataPath, 'tmp/*/claims']);
            if ~isempty(d)
                system(sprintf('ln -s "%s:n1:5:1" "%s/t%d"', job, d(1).folder, f));
            end
        end

        function removeHeldMemory(~)
            % the memory a 'hold' state task keeps
            if isappdata(0, 'pkRunnerStateHold')
                rmappdata(0, 'pkRunnerStateHold');
            end
        end

        function copyCallFile(~, dataPath, name)
            % copies a file of the (only) MATLAB-job call in dataPath there
            d = dir([dataPath, 'task_runner_*/', name]);
            if ~isempty(d)
                copyfile(fullfile(d(1).folder, d(1).name), [dataPath, name]);
            end
        end

        function claimInCall(~, dataPath, f, job)
            % claims task f of the (only) MATLAB-job call in dataPath for job
            d = dir([dataPath, 'task_runner_*']);
            if ~isempty(d)
                system(sprintf('ln -s "%s:n1:5:1" "%s/%s/claims/t%d"', job, d(1).folder, d(1).name, f));
            end
        end

        function [inFns, outFns, logFn, dataPath] = makeTasks(testCase, n)
            % n input text files and the output paths in a separate folder
            dataPath = testCase.createTempFolder();
            mkdir([dataPath, 'out']);
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : n, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout/out_%d.txt', dataPath, x), 1 : n, 'unif', 0);
            for f = 1 : n
                testCase.writeLines(inFns{f}, {sprintf('task %d', f)});
            end
            logFn = [dataPath, 'calls.log'];
        end

        function [lines] = readLog(~, logFn)
            % output paths of all pktest.frameworkTask calls, in call order
            lines = cell(0, 1);
            if exist(logFn, 'file')
                lines = splitlines(strip(fileread(logFn)));
                lines = lines(~cellfun(@isempty, lines));
            end
        end

        function [args] = readNul(~, fn)
            % NUL separated fields of a file
            fid = fopen(fn, 'r');
            raw = fread(fid, inf, 'uint8=>uint8')';
            fclose(fid);
            args = strsplit(native2unicode(raw(1 : end - 1), 'UTF-8'), char(0), 'CollapseDelimiters', false);
        end

        function writeLines(~, fn, lines)
            fid = fopen(fn, 'w');
            if ~isempty(lines)
                fwrite(fid, [strjoin(lines, newline), newline]);
            end
            fclose(fid);
        end
    end
end


function [txt] = fileread_or_empty(fn)
% the text of a file, '' if it doesn't exist

txt = '';
if isfile(fn)
    txt = fileread(fn);
end

end
