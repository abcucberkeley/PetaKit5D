classdef SingleJobRunnerTest < pktest.TestCase
    % task runners for the single jobs of generic_single_job_submit_wrapper:
    % the single jobs a MATLAB session submits with the same function,
    % resources and compiled master (or MATLAB launch) share a pool (see
    % task_runner_run), so a job, after its task (which runs as before),
    % runs unstarted tasks of the pool's other jobs in a task runner, and a
    % job whose task a runner of another job has ends without starting a
    % process. A pool is made with the second task of a kind (the first
    % runs as before), and a runner takes a task only while Slurm keeps its
    % job in the queue (pending with a reason). The callers' status checks
    % count such a task as running, and the master takes a task back from
    % the pool before it submits it again (as before) or runs it. The tests
    % start a kind with a first task (startKind), so that tasks 1, 2, ...
    % are the pool's tasks 1, 2, ... with jobs 1000, 1001, ...
    % The tests give the wrapper their job log folder (jobLogDir, as the
    % callers do), so that pools and logs go there and never to the default
    % one, ../job_logs beside the PetaKit5D folder, where it is there.
    % Compiled jobs use a stand-in build (pktest.fakeRunnerBuild); the jobs'
    % commands from a fake Slurm (pktest.FakeSlurmFixture) run in bash as
    % Slurm would run them, with Slurm's job variables.

    properties (Constant)
        EnvNames = {'PETAKIT5D_RUNNER', 'PETAKIT5D_TASK_METRICS', 'PETAKIT5D_TASK', 'PETAKIT5D_TASK_GPU', ...
            'PETAKIT5D_RUNNER_MEM_GB', 'PETAKIT5D_RUNNER_MINUTES', 'SLURM_JOB_END_TIME', 'SLURM_JOB_START_TIME'}
    end

    methods (TestMethodSetup)
        function clearSession(testCase)
            % the session's pools and the config files are kept in
            % persistent variables
            clear task_runner_run generic_single_job_submit_wrapper single_job_pools
            testCase.addTeardown(@clear, 'task_runner_run', 'generic_single_job_submit_wrapper', 'single_job_pools');
        end

        function restoreEnvironment(testCase)
            for i = 1 : numel(testCase.EnvNames)
                testCase.addTeardown(@setenv, testCase.EnvNames{i}, getenv(testCase.EnvNames{i}));
                setenv(testCase.EnvNames{i}, '');
            end
        end
    end

    methods (Test, TestTags = {'Unit'})
        function compiledJobsRunTheUnstartedTasksOfTheirPool(testCase)
            % four single jobs of the same kind: the first one to start runs
            % its task as before (a process of its own), then the tasks of
            % the others in their order in one more process (a task
            % runner), and cancels their jobs while they are pending; if
            % they start, they end without starting a process. The jobs
            % have the resources of before
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            for k = 1 : 4
                [~, ~, submitted] = testCase.submit(c, k, 'ok', -1);
                testCase.assertEqual(submitted, 1);
            end
            for k = 2 : 4
                slurm.setJobStates(sprintf('%d_%d', 999 + k, k), {'PD'});
            end
            cmds = slurm.calls('sbatch');
            testCase.assertNumElements(cmds, 4);
            for k = 1 : 4
                testCase.verifyNotEmpty(regexp(cmds{k}, sprintf(['^--array=%d -o \\S+ -e \\S+ --cpus-per-task=1 -p abc ', ...
                    '.*-t 24:0:00 +--wrap=date; echo Matlab command: +"pk_task\\(.*o_%d\\.txt.*\\)"; ', ...
                    'bash \\S+/task_runner_single_\\S+/run_tasks\\.sh job %d 1 0 3 1$'], k, k, k), 'once'), cmds{k});
            end
            testCase.verifyNumElements(dir([c.logDir, 'task_runner_single_*']), 1, 'The jobs are not in one pool.');

            out = testCase.runJob(c, slurm, 1);
            runs = testCase.taskRuns(c);
            testCase.verifyEqual([runs.task], [1, 2, 3, 4], out);
            testCase.verifyEqual([runs.runner], [0, 1, 1, 1], 'The own task did not run as before.');
            testCase.verifyNumElements(unique([runs(2 : end).pid]), 1, 'The other tasks did not run in one process.');
            testCase.verifyNotEqual(runs(1).pid, runs(end).pid);
            testCase.verifyEqual(slurm.calls('scancel'), {'--state=PENDING 1001_2'; '--state=PENDING 1002_3'; ...
                '--state=PENDING 1003_4'});
            for k = 2 : 4
                out = testCase.runJob(c, slurm, k);
                testCase.verifySubstring(out, 'no task left');
            end
            testCase.verifyNumElements(testCase.taskRuns(c), 4, 'A task ran more than once.');
            testCase.verifyTrue(all(cellfun(@isfile, c.outFns(1 : 4))));
        end

        function firstTaskOfAKindRunsAsBefore(testCase)
            % a session's first task of a kind runs as before: a session with
            % one task pays nothing (no pool; neither its submission nor its
            % status checks load the runner code); the second task makes the
            % pool, with itself as the pool's task 1
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            clear task_runner_run
            testCase.submit(c, 1, 'ok', -1);
            testCase.verifyEmpty(dir([c.logDir, 'task_runner_single_*']), 'The first task made a pool.');
            testCase.verifyEqual(check_single_job_status(1000, 1), -1);
            testCase.verifyFalse(ismember('task_runner_run', inmem), 'The first task loaded the runner code.');
            testCase.submit(c, 2, 'ok', -1);
            % (the commands of before end with a line break)
            cmds = slurm.calls('sbatch');
            cmds = cmds(~cellfun(@isempty, strtrim(cmds)));
            testCase.assertNumElements(cmds, 2);
            testCase.verifyEmpty(strfind(cmds{1}, 'run_tasks.sh'), cmds{1});
            testCase.verifyNotEmpty(regexp(cmds{2}, 'bash \S+/task_runner_single_\S+/run_tasks\.sh job 1 1 0 3 1$', 'once'), cmds{2});
            pool = testCase.poolDir(c);
            testCase.verifyEqual(strip(fileread([pool, '/runnable.txt'])), '1');
            testCase.verifyEqual(strip(fileread([pool, '/jobs.txt'])), '1 1001_2');
        end

        function runnersTakeOnlyTasksWhoseJobsWaitInTheQueue(testCase)
            % a job's runner takes another job's task only while Slurm keeps
            % that job in the queue (pending with a reason: resources,
            % priority, a limit), not a job its scheduler hasn't looked at
            % yet (reason None: it may start at once; in idle streaming,
            % runners took the tasks of tiles submitted a moment before), a
            % job that runs, or a task whose job is not recorded yet
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            for k = 1 : 6
                testCase.submit(c, k, 'ok', -1);
            end
            % 2: just submitted; 3: pending (resources); 4: its job runs;
            % 5: pending, its job not recorded; 6: pending (priority)
            states = {'PD:None', 'PD', 'R', 'PD:Priority', 'PD:Priority'};
            for k = 2 : 6
                slurm.setJobStates(sprintf('%d_%d', 999 + k, k), states(k - 1));
            end
            pool = testCase.poolDir(c);
            jobs = splitlines(strip(fileread([pool, '/jobs.txt'])));
            testCase.writeLines([pool, '/jobs.txt'], jobs(~startsWith(jobs, '5 ')));
            out = testCase.runJob(c, slurm, 1);
            runs = testCase.taskRuns(c);
            testCase.verifyEqual([runs.task], [1, 3, 6], out);
            testCase.verifyEqual([runs.runner], [0, 1, 1]);
            testCase.verifyEqual(slurm.calls('scancel'), {'--state=PENDING 1002_3'; '--state=PENDING 1005_6'});
        end

        function jobWhoseTaskARunnerHasEndsAtOnce(testCase)
            % a job whose own task another job's runner took ends without
            % starting a process, and takes no other task (those have jobs
            % of their own)
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            for k = 1 : 3
                testCase.submit(c, k, 'ok', -1);
            end
            system(sprintf('ln -s "99_1:othernode:1:1" "%s/claims/t2"', testCase.poolDir(c)));
            out = testCase.runJob(c, slurm, 2);
            testCase.verifySubstring(out, 'no task left');
            testCase.verifyEmpty(testCase.taskRuns(c));
        end

        function slotStartsWithItsOwnTask(testCase)
            % the first task of a single job's slot is its own: if another
            % process took it in the meantime, the slot takes no other task
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            for k = 1 : 3
                testCase.submit(c, k, 'ok', -1);
                slurm.setJobStates(sprintf('%d_%d', 999 + k, k), {'PD'});
            end
            pool = testCase.poolDir(c);
            system(sprintf('ln -s "99_1:othernode:1:1" "%s/claims/t2"', pool));
            L = testCase.createTempFolder();
            [status, out] = system(sprintf('env -u PETAKIT5D_RUNNER bash "%s/run_tasks.sh" slot 2 0 3 1 "%s" 1 1 2>&1', ...
                pool, L(1 : end - 1)));
            testCase.assertEqual(status, 0, out);
            testCase.verifyEmpty(testCase.taskRuns(c), 'The slot took a task of another job first.');
        end

        function matlabJobGoesOnAfterAProcessEnded(testCase)
            % a MATLAB job whose process ended early (its memory, time or
            % open file guard) starts another one for the pool's next tasks
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [funcStrs, ~, logDir] = testCase.matlabTasks(3);
            % a stand-in MATLAB: records the task and ends as a recycled runner
            fakeFn = [logDir, 'fake_matlab.sh'];
            recFn = [logDir, 'started.log'];
            testCase.writeLines(fakeFn, {
                'args=$(echo "$1" | sed -n "s/.*task_runner(''matlab'',''[^'']*'',''[^'']*'',[0-9]*,''\([^'']*\)'',''[^'']*'',\([0-9]*\)).*/\1 \2/p")'
                'read -r st f <<< "$args"'
                sprintf('echo "$f" >> "%s"', recFn)
                'printf "d %s\nt 0.1\nx recycle\n" "$f" >> "$st"'});
            testCase.startMatlabKind(logDir, slurm, 'MatlabLaunchStr', sprintf('bash %s', fakeFn));
            for k = 1 : 3
                testCase.submitMatlab(funcStrs{k}, logDir, k, 'MatlabLaunchStr', sprintf('bash %s', fakeFn));
                slurm.setJobStates(sprintf('%d_%d', 999 + k, k), {'PD'});
            end
            [status, out] = system(sprintf('env -u PETAKIT5D_RUNNER TMPDIR="%s" bash "%s/run_tasks.sh" matlab 1 3 2>&1', ...
                tempdir, testCase.poolDir(logDir)));
            testCase.assertEqual(status, 0, out);
            testCase.verifyEqual(sscanf(fileread(recFn), '%d')', [1, 2, 3], out);
        end

        function jobWithNothingToTakeRunsItsTaskAsBefore(testCase)
            % with no unstarted task in its pool, a job runs its task line
            % as before (a process for the task, no task runner)
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            testCase.submit(c, 1, 'ok', -1);
            testCase.runJob(c, slurm, 1);
            runs = testCase.taskRuns(c);
            testCase.assertNumElements(runs, 1);
            testCase.verifyEqual(runs.runner, 0, 'The process was a task runner.');
        end

        function tasksSubmittedWhileAJobRunsAreTaken(testCase)
            % (streaming) a task added while a job runs is taken by it (by
            % its task runner, after the tasks submitted before it)
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            testCase.submit(c, 1, 'slow', -1);
            testCase.submit(c, 2, 'ok', -1);
            slurm.setJobStates('1001_2', {'PD'});
            slurm.setJobStates('1002_3', {'PD'});
            doneFn = [c.root, 'job_1.done'];
            testCase.runJob(c, slurm, 1, doneFn);
            t0 = tic;
            while isempty(testCase.taskRuns(c)) && toc(t0) < 20
                pause(0.05);
            end
            testCase.submit(c, 3, 'ok', -1);
            while ~isfile(doneFn) && toc(t0) < 30
                pause(0.05);
            end
            testCase.assertTrue(isfile(doneFn), 'The job did not end.');
            runs = testCase.taskRuns(c);
            testCase.verifyEqual([runs.task], [1, 2, 3]);
            testCase.verifyEqual([runs.runner], [0, 1, 1]);
            testCase.verifyNumElements(unique([runs(2 : end).pid]), 1);
        end

        function poolsKeepKindsApart(testCase)
            % a pool has the tasks of one function, resources and master;
            % another CPU count, function or job log folder makes another
            % (each with its second task)
            testCase.assumeLinux();
            testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            kinds = {{}, {'cpusPerTask', 2}, {'funcName', 'pk_other'}, {}, {'cpusPerTask', 2}, {}, ...
                {'cpusPerTask', 2}, {'funcName', 'pk_other'}};
            for k = 1 : 8
                testCase.submit(c, k, 'ok', -1, kinds{k}{:});
            end
            pools = dir([c.logDir, 'task_runner_single_*']);
            testCase.assertNumElements(pools, 3);
            nTasks = arrayfun(@(x) numel(splitlines(strip(fileread(fullfile(x.folder, x.name, 'runnable.txt'))))), pools);
            testCase.verifyEqual(sort(nTasks(:))', [1, 2, 2]);
        end

        function offSwitchOldBuildsAndRelativeLogFoldersSubmitAsBefore(testCase)
            % taskRunner false (option or config key), a compiled master
            % built without task runners, or a job log folder that is not an
            % absolute path: the job command of before also for the second
            % task of a kind (which would join a pool), and no pool
            testCase.assumeLinux();
            % setup arguments, a build with runners, submission arguments
            % (the caller's false wins over the config's true; a first
            % submission of a task with a job of before is submitted again;
            % relative log names stay so with a job log folder that isn't
            % there, as with the default ../job_logs where it isn't there)
            cases = {{}, true, {'taskRunner', false}; {}, false, {}; {'configTaskRunner', false}, true, {}; ...
                {}, true, {'jobLogFname', 'job_%A_%a.out', 'jobLogDir', 'no_job_log_folder'}; ...
                {'configTaskRunner', true}, true, {'taskRunner', false}; {}, true, {'jobId', 555}};
            for i = 1 : size(cases, 1)
                slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
                c = testCase.compiledSetup(cases{i, 2}, cases{i, 1}{:});
                jobId = -1;
                args = cases{i, 3};
                if ~isempty(args) && strcmp(args{1}, 'jobId')
                    jobId = args{2};
                    args = {};
                end
                testCase.submit(c, 1, 'ok', -1, args{:});
                testCase.submit(c, 2, 'ok', jobId, args{:});
                % (the commands of before end with a line break)
                cmds = slurm.calls('sbatch');
                cmds = cmds(~cellfun(@isempty, strtrim(cmds)));
                testCase.assertNumElements(cmds, 2);
                testCase.verifyNotEmpty(regexp(cmds{2}, sprintf(['--wrap=date; echo Matlab command: +"pk_task\\(.*\\)"; ', ...
                    'MCR_CACHE_ROOT=\\S+/mcr_cache_\\$\\{SLURM_JOB_ID:-0\\} %s %s pk_task ''\\S+/o_2\\.txt'' ''ok''\\s*$'], ...
                    c.b.master, c.b.mcr), 'once'), sprintf('Case %d: %s', i, cmds{2}));
                testCase.verifyEmpty(dir([c.logDir, 'task_runner_single_*']), sprintf('Case %d made a pool.', i));
                clear task_runner_run generic_single_job_submit_wrapper single_job_pools
            end
        end

        function taskTakenBetweenSubmissionTriesIsSubmittedAsBefore(testCase)
            % sbatch fails once, and meanwhile a runner (whose job ended
            % since) claimed the new task: the next try submits it as before
            % (its own job could not run it)
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            % an sbatch that, the first time, has the task claimed and fails
            stubDir = [c.root, 'stub/'];
            mkdir(stubDir);
            testCase.writeLines([stubDir, 'sbatch'], {'#!/bin/sh'
                sprintf('if [ ! -f "%sfailed" ]; then', stubDir)
                sprintf('    touch "%sfailed"', stubDir)
                sprintf('    for d in "%s"task_runner_single_*/claims; do ln -s "999_9:othernode:5:1" "$d/t1"; done', c.logDir)
                '    exit 1'
                'fi'
                sprintf('exec "%s/sbatch" "$@"', slurm.Folder)});
            fileattrib([stubDir, 'sbatch'], '+x');
            testCase.addTeardown(@setenv, 'PATH', getenv('PATH'));
            setenv('PATH', [stubDir, pathsep, getenv('PATH')]);
            [~, ~, submitted] = testCase.submit(c, 1, 'ok', -1, 'maxTrialNum', 2);
            testCase.verifyEqual(submitted, 1);
            cmds = slurm.calls('sbatch');
            testCase.assertNumElements(cmds, 1);
            testCase.verifyEmpty(strfind(cmds{1}, 'run_tasks.sh'), cmds{1});
            owner = testCase.claimOf(testCase.poolDir(c), 1);
            testCase.verifyEqual(owner(end - 1 : end), ':0', 'The master did not take the task.');
        end

        function poolThatCannotBeMadeLeavesJobsAsBefore(testCase)
            % a job log folder where no pool can be made (read-only): a
            % warning, and the jobs are submitted as before
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            fileattrib(c.logDir, '-w', 'a');
            testCase.addTeardown(@fileattrib, c.logDir, '+w', 'u');
            testCase.verifyWarning(@() testCase.submit(c, 1, 'ok', -1), 'task_runner_run:single');
            testCase.verifyWarningFree(@() testCase.submit(c, 2, 'ok', -1));
            % (the commands of before end with a line break)
            cmds = slurm.calls('sbatch');
            cmds = cmds(~cellfun(@isempty, strtrim(cmds)));
            testCase.assertNumElements(cmds, 2);
            testCase.verifyEmpty(cell2mat(strfind(cmds, 'run_tasks.sh')), 'A job was submitted with task runners.');
        end

        function claimsFolderThatFailsLeavesTasksToTheMaster(testCase)
            % a claims folder that takes no new claims (here read-only):
            % the master can't claim the task of an ended job either, but as
            % no process has it, it is submitted again as before
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            testCase.submit(c, 1, 'ok', -1);
            claimsDir = [testCase.poolDir(c), '/claims'];
            fileattrib(claimsDir, '-w', 'a');
            testCase.addTeardown(@fileattrib, claimsDir, '+w', 'u');
            [jobId, ~, submitted] = testCase.submit(c, 1, 'ok', 1000);
            testCase.verifyEqual([jobId, submitted], [1001, 1]);
            cmds = slurm.calls('sbatch');
            testCase.verifyEmpty(strfind(cmds{end}, 'run_tasks.sh'), cmds{end});
        end

        function taskWithALineBreakIsSubmittedAsBefore(testCase)
            % a pool keeps a task's line in a line of its task file, so a
            % task string with a line break is submitted as before
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            testCase.submit(c, 1, sprintf('o\nk'), -1);
            testCase.verifyEqual(slurm.callCount('sbatch') > 0, true);
            log = fileread([slurm.Folder, '/sbatch.log']);
            testCase.verifyEmpty(strfind(log, 'run_tasks.sh'), log);
            testCase.verifySubstring(log, 'run_mccMaster.sh');
            testCase.verifyEmpty(dir([c.logDir, 'task_runner_single_*']), 'The task joined a pool.');
        end

        function matlabJobsRunTheJobScript(testCase)
            % MATLAB jobs: the job script starts MATLAB (with the launch
            % string of before) only for a task it claimed; off: the
            % command of before
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [funcStrs, ~, logDir] = testCase.matlabTasks(3);
            testCase.startMatlabKind(logDir, slurm, 'MatlabLaunchStr', 'echo STARTED');
            for k = 1 : 3
                testCase.submitMatlab(funcStrs{k}, logDir, k, 'MatlabLaunchStr', 'echo STARTED');
            end
            cmds = slurm.calls('sbatch');
            testCase.assertNumElements(cmds, 3);
            for k = 1 : 3
                testCase.verifyNotEmpty(regexp(cmds{k}, sprintf('; bash \\S+/run_tasks\\.sh matlab %d 3$', k), 'once'), cmds{k});
            end
            pool = dir([logDir, 'task_runner_single_*']);
            testCase.assertNumElements(pool, 1);
            system(sprintf('ln -s "99_1:othernode:1:1" "%s/%s/claims/t2"', logDir, pool.name));
            out = testCase.runWrap(cmds{2}, 1001, 2, logDir);
            testCase.verifyFalse(contains(out, 'STARTED'), out);
            testCase.verifySubstring(out, 'no task left');
            out = testCase.runWrap(cmds{1}, 1000, 1, logDir);
            testCase.verifySubstring(out, 'STARTED setup([]);t0_=tic;task_runner(''matlab''');

            clear task_runner_run generic_single_job_submit_wrapper single_job_pools
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.submitMatlab(funcStrs{1}, logDir, 1, 'taskRunner', false);
            testCase.submitMatlab(funcStrs{2}, logDir, 2, 'taskRunner', false);
            cmds = slurm.calls('sbatch');
            testCase.verifySubstring(cmds{2}, sprintf('-nojvm -batch "setup([]);t0_=tic;%s;toc(t0_)"', funcStrs{2}));
        end

        function matlabRunnerTakesTaskStringsFromThePool(testCase)
            % a MATLAB job's runner runs its task, then the pool's other
            % tasks in their order, each from the string its job got
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [funcStrs, outFns, logDir, logFn, inFns] = testCase.matlabTasks(3);
            testCase.startMatlabKind(logDir, slurm);
            for k = 1 : 3
                testCase.submitMatlab(funcStrs{k}, logDir, k);
                slurm.setJobStates(sprintf('%d_%d', 999 + k, k), {'PD'});
            end
            testCase.runMatlabRunner(logDir, 1);
            testCase.verifyEqual(testCase.readLog(logFn), outFns([1, 2, 3])');
            testCase.verifyEqual(cellfun(@fileread, outFns, 'unif', 0), cellfun(@fileread, inFns, 'unif', 0));
        end

        function matlabPoolTasksRecordTheirMetrics(testCase)
            % with task metrics, a task a MATLAB runner runs records itself
            % as the task of its single job
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [funcStrs, ~, logDir] = testCase.matlabTasks(2);
            testCase.startMatlabKind(logDir, slurm, 'taskMetrics', true, 'cpusPerTask', 2, 'memAllocate', 30);
            testCase.submitMatlab(funcStrs{1}, logDir, 7, 'taskMetrics', true, 'cpusPerTask', 2, 'memAllocate', 30);
            testCase.submitMatlab(funcStrs{2}, logDir, 8, 'taskMetrics', true, 'cpusPerTask', 2, 'memAllocate', 30);
            testCase.assertNumElements(slurm.calls('sbatch'), 2);
            slurm.setJobStates('1001_8', {'PD'});
            setenv('SLURM_ARRAY_JOB_ID', '1000');
            setenv('SLURM_ARRAY_TASK_ID', '7');
            testCase.addTeardown(@setenv, 'SLURM_ARRAY_JOB_ID', getenv('SLURM_ARRAY_JOB_ID'));
            testCase.addTeardown(@setenv, 'SLURM_ARRAY_TASK_ID', getenv('SLURM_ARRAY_TASK_ID'));
            testCase.runMatlabRunner(logDir, 1);
            t = sortrows(read_task_metrics([logDir, 'task_metrics']), 'task');
            testCase.assertEqual(height(t), 2);
            testCase.verifyEqual(t.task, [7; 8]);
            testCase.verifyEqual([t.mode, t.status], ["slurm_matlab", "ok"; "slurm_matlab", "ok"]);
            testCase.verifyEqual([t.cpusPerTask, t.memAllocate], [2, 30; 2, 30], 'The submissions were not joined.');
        end

        function liveRunnerKeepsTheTaskOfAnEndedJob(testCase)
            % after the own job of a task ended, a runner of another job
            % that took the task may still run it: the task counts as
            % running (no new job) until that job ends; then it is
            % submitted again as before (no task runner), and the master
            % has it in the pool (no runner takes it)
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            testCase.submit(c, 1, 'ok', -1);
            testCase.submit(c, 2, 'ok', -1);
            pool = testCase.poolDir(c);
            system(sprintf('ln -s "1000_1:othernode:5:1" "%s/claims/t2"', pool));
            slurm.setJobStates('1000_1', {'R', 'R', 'done'});
            for i = 1 : 2
                [jobId, status, submitted] = testCase.submit(c, 2, 'ok', 1001);
                testCase.verifyEqual([jobId, status, submitted], [1001, 1, 0], sprintf('Check %d.', i));
            end
            [jobId, ~, submitted] = testCase.submit(c, 2, 'ok', 1001);
            testCase.verifyEqual([jobId, submitted], [1002, 1]);
            cmds = slurm.calls('sbatch');
            testCase.assertNumElements(cmds, 3);
            testCase.verifyEmpty(strfind(cmds{3}, 'run_tasks.sh'), 'The task was submitted again with task runners.');
            testCase.verifySubstring(cmds{3}, 'pk_task ''');
            [owner, marked] = testCase.claimOf(pool, 2);
            testCase.verifyEqual(owner(end - 1 : end), ':0', 'The master did not take the task.');
            testCase.verifyTrue(marked, 'The task is not marked submitted again.');
        end

        function runnerKeepsItsTaskWhileSqueueFails(testCase)
            % a failed squeue (e.g., a timeout of a busy controller) read as
            % the job of a runner that took a pool task gone, so the master
            % took the task back and ran or submitted it again; the runner
            % keeps it now
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            testCase.submit(c, 1, 'ok', -1);
            testCase.submit(c, 2, 'ok', -1);
            pool = testCase.poolDir(c);
            system(sprintf('ln -s "1000_1:othernode:5:1" "%s/claims/t2"', pool));
            slurm.setJobStates('1000_1', {'R'});
            funcStr = sprintf('pk_task(''%s'',''ok'')', c.outFns{2});
            slurm.failQueries(1);
            [live, claimed] = task_runner_run('single_live', funcStr);
            testCase.verifyEqual([live, claimed], [true, true]);
            slurm.failQueries(1);
            testCase.verifyFalse(task_runner_run('single_take', funcStr), 'The master took the task of a runner squeue failed to answer about.');
            [owner, marked] = testCase.claimOf(pool, 2);
            testCase.verifyEqual(owner, '1000_1:othernode:5:1');
            testCase.verifyFalse(marked);
        end

        function pendingJobGetsBackTheTaskOfARunnerThatEnded(testCase)
            % a runner of another job took the task of a pending job and
            % ended without finishing it (its job was killed): the master
            % removes the claim, so that the pending job runs the task when
            % it starts (as without runners); nothing is submitted
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            testCase.submit(c, 1, 'ok', -1);
            testCase.submit(c, 2, 'ok', -1);
            pool = testCase.poolDir(c);
            % (task 1 runs elsewhere)
            system(sprintf('ln -s "77_1:othernode:5:1" "%s/claims/t1"', pool));
            system(sprintf('ln -s "1000_1:othernode:5:1" "%s/claims/t2"', pool));
            slurm.setJobStates('1001_2', {'PD'});
            [jobId, status, submitted] = testCase.submit(c, 2, 'ok', 1001);
            testCase.verifyEqual([jobId, status, submitted], [1001, 0, 0]);
            testCase.verifyEmpty(testCase.claimOf(pool, 2), 'The claim of the runner that ended is left.');
            testCase.runJob(c, slurm, 2);
            testCase.verifyEqual([testCase.taskRuns(c).task], 2);
            testCase.verifyEqual(slurm.callCount('sbatch'), 2);
            testCase.verifyEqual(slurm.callCount('scancel'), 0);
        end

        function pendingJobOfATaskARunnerFailedIsReplaced(testCase)
            % a runner of another job failed the task of a pending job: the
            % master takes the task back, cancels the pending job and
            % submits the task again as before (as a failed job of its own)
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            for k = 1 : 3
                testCase.submit(c, k, 'ok', -1);
            end
            pool = testCase.poolDir(c);
            % (the runner's job still runs: it ends after the failure)
            slurm.setJobStates('1000_1', {'R'});
            for k = 2 : 3
                system(sprintf('ln -s "1000_1:othernode:5:1" "%s/claims/t%d"', pool, k));
                system(sprintf('ln -s "1000_1:othernode:5:1" "%s/claims/f%d"', pool, k));
                slurm.setJobStates(sprintf('%d_%d', 999 + k, k), {'PD'});
            end
            [jobId, ~, submitted] = testCase.submit(c, 2, 'ok', 1001);
            testCase.verifyEqual([jobId, submitted], [1003, 1]);
            testCase.verifyEqual(slurm.calls('scancel'), {'1001_2'});
            cmds = slurm.calls('sbatch');
            testCase.verifyEmpty(strfind(cmds{end}, 'run_tasks.sh'), cmds{end});
            % the callers' status check: the task has ended (its last try)
            testCase.verifyEqual(check_single_job_status(1002, 3), -1);
            testCase.verifyEqual(slurm.calls('scancel'), {'1001_2'; '1002_3'});
            owner = testCase.claimOf(pool, 3);
            testCase.verifyEqual(owner(end - 1 : end), ':0', 'The master did not take the task.');
        end

        function failedTaskIsSubmittedAgainAsBefore(testCase)
            % a task that failed in its job is submitted again as before
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            testCase.submit(c, 1, 'error', -1);
            testCase.runJob(c, slurm, 1);
            testCase.verifyTrue(isfolder(testCase.poolDir(c)));
            [owner, ~, failed] = testCase.claimOf(testCase.poolDir(c), 1);
            testCase.verifyTrue(startsWith(owner, '1000_1:') && failed, 'The failed task is not marked.');
            [jobId, ~, submitted] = testCase.submit(c, 1, 'error', 1000);
            testCase.verifyEqual([jobId, submitted], [1001, 1]);
            cmds = slurm.calls('sbatch');
            testCase.verifyEmpty(strfind(cmds{2}, 'run_tasks.sh'), cmds{2});
            % (the job of before runs it again; the task fails again)
            testCase.runJob(c, slurm, 2, '', true);
            testCase.verifyEqual([testCase.taskRuns(c).task], [1, 1]);
        end

        function singleJobStatusCountsTaskRunners(testCase)
            % check_single_job_status: a task that a runner of another job
            % has counts as running; once none has, the master takes it
            % from the pool (-1), so that no job's runner takes it after
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            for k = 1 : 3
                testCase.submit(c, k, 'ok', -1);
            end
            pool = testCase.poolDir(c);
            system(sprintf('ln -s "1000_1:othernode:5:1" "%s/claims/t2"', pool));
            slurm.setJobStates('1000_1', {'R', 'done'});
            slurm.setJobStates('1002_3', {'PD'});
            testCase.verifyEqual(check_single_job_status(1001, 2), 1);
            testCase.verifyEqual(check_single_job_status(1001, 2), -1);
            testCase.verifyEqual(check_single_job_status(1002, 3), 0);
            % a job's task the master took: its job (or a runner) leaves it
            testCase.verifyEqual(check_single_job_status(1000, 1), -1);
            out = testCase.runJob(c, slurm, 1);
            testCase.verifySubstring(out, 'no task left');
            testCase.verifyEmpty(testCase.taskRuns(c));
            % a job of no pool: as check_slurm_job_status
            slurm.setJobStates('77_4', {'R', 'done'});
            testCase.verifyEqual([check_single_job_status(77, 4), check_single_job_status(77, 4)], [1, -1]);
        end

        function masterTakesTheLastFileFromThePool(testCase)
            % masterCompute with the last file: the master cancels the
            % file's pending job and runs it itself, after it took the task
            % from the pool; a task a runner of another job has stays with it
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            testCase.submit(c, 1, 'ok', -1);
            testCase.submit(c, 2, 'ok', -1);
            pool = testCase.poolDir(c);
            slurm.setJobStates('1001_2', {'PD'});
            slurm.setJobStates('1000_1', {'R'});
            system(sprintf('ln -s "1000_1:othernode:5:1" "%s/claims/t2"', pool));
            [~, status] = testCase.submit(c, 2, 'ok', 1001, 'masterCompute', true, 'lastFile', true);
            testCase.verifyEqual(status, 1);
            testCase.verifyEqual(slurm.callCount('scancel'), 0);
            system(sprintf('rm -f "%s/claims/t2"', pool));
            [~, status] = testCase.submit(c, 2, 'ok', 1001, 'masterCompute', true, 'lastFile', true);
            testCase.verifyEqual(status, 0);
            testCase.verifyEqual(slurm.calls('scancel'), {'1001_2'});
            owner = testCase.claimOf(pool, 2);
            testCase.verifyEqual(owner(end - 1 : end), ':0', 'The master did not take the task.');
            testCase.verifyEqual(slurm.callCount('sbatch'), 2);
        end

        function failedSubmissionLeavesTheTaskToTheMaster(testCase)
            % a task whose job couldn't be submitted (sbatch failed each
            % try) is taken back from the pool; the caller's next try
            % submits it as before
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            slurm.rejectSubmissions(2);
            [jobId, ~, submitted] = testCase.submit(c, 1, 'ok', -1, 'maxTrialNum', 2);
            testCase.verifyEqual([jobId, submitted], [-1, 0]);
            owner = testCase.claimOf(testCase.poolDir(c), 1);
            testCase.verifyEqual(owner(end - 1 : end), ':0', 'The master did not take the task.');
            [~, ~, submitted] = testCase.submit(c, 1, 'ok', -1);
            testCase.verifyEqual(submitted, 1);
            cmds = slurm.calls('sbatch');
            testCase.verifyEmpty(strfind(cmds{end}, 'run_tasks.sh'), cmds{end});
        end

        function endedPoolGoesWithItsLastJob(testCase)
            % at the end of a call (single_job_pools('end', mark)) the
            % master leaves the pools the call made, which take no more
            % tasks: a job of a pool still in the queue keeps it, and the
            % last of its jobs to end removes it; a later task of the kind
            % makes a new pool
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            mark = single_job_pools('mark');
            testCase.startKind(c, slurm);
            for k = 1 : 2
                testCase.submit(c, k, 'ok', -1);
            end
            D = testCase.poolDir(c);
            slurm.setJobStates('1001_2', {'R'});
            out = testCase.runJob(c, slurm, 1);
            single_job_pools('end', mark);
            testCase.verifyTrue(isfolder(D), 'The pool went while a job of it ran.');
            out = [out, testCase.runJob(c, slurm, 2)];
            testCase.verifyFalse(isfolder(D), out);
            testCase.verifyEqual([testCase.taskRuns(c).task], [1, 2], out);
            testCase.submit(c, 3, 'ok', -1);
            testCase.verifyNotEqual(testCase.poolDir(c), D, 'The ended pool took a task.');
        end

        function endedPoolWithoutJobsGoesAtOnce(testCase)
            % a pool whose jobs are no longer in the queue goes when the
            % master leaves it, also with the marks of jobs that never ran
            % (a runner took their tasks and cancelled them)
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            mark = single_job_pools('mark');
            testCase.startKind(c, slurm);
            for k = 1 : 3
                testCase.submit(c, k, 'ok', -1);
            end
            for k = 2 : 3
                slurm.setJobStates(sprintf('%d_%d', 999 + k, k), {'PD'});
            end
            D = testCase.poolDir(c);
            out = testCase.runJob(c, slurm, 1);
            testCase.verifyEqual([testCase.taskRuns(c).task], [1, 2, 3], out);
            % (cancelled)
            for k = 2 : 3
                slurm.setJobStates(sprintf('%d_%d', 999 + k, k), {'done'});
            end
            testCase.verifyTrue(isfolder(D), 'The pool went before the master left it.');
            single_job_pools('end', mark);
            testCase.verifyFalse(isfolder(D), 'The pool is still there.');
        end

        function pendingJobKeepsAnEndedPool(testCase)
            % a job whose task a runner took, but whose cancel was missed,
            % is still pending when the master leaves the pool: the pool
            % stays until the job starts, ends at once (its task is taken)
            % and removes it
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            mark = single_job_pools('mark');
            testCase.startKind(c, slurm);
            for k = 1 : 2
                testCase.submit(c, k, 'ok', -1);
            end
            slurm.setJobStates('1001_2', {'PD'});
            D = testCase.poolDir(c);
            out = testCase.runJob(c, slurm, 1);
            testCase.verifyEqual([testCase.taskRuns(c).task], [1, 2], out);
            single_job_pools('end', mark);
            testCase.verifyTrue(isfolder(D), 'The pool went while a job of it was pending.');
            out = testCase.runJob(c, slurm, 2);
            testCase.verifySubstring(out, 'no task left');
            testCase.verifyFalse(isfolder(D), out);
            testCase.verifyNumElements(testCase.taskRuns(c), 2, 'A task ran more than once.');
        end

        function matlabPoolGoesWithItsLastJob(testCase)
            % the same for a pool of MATLAB jobs
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [funcStrs, ~, logDir] = testCase.matlabTasks(2);
            mark = single_job_pools('mark');
            testCase.startMatlabKind(logDir, slurm, 'MatlabLaunchStr', 'echo STARTED');
            for k = 1 : 2
                testCase.submitMatlab(funcStrs{k}, logDir, k, 'MatlabLaunchStr', 'echo STARTED');
            end
            D = testCase.poolDir(logDir);
            slurm.setJobStates('1000_1', {'R'});
            single_job_pools('end', mark);
            testCase.verifyTrue(isfolder(D), 'The pool went while a job of it ran.');
            cmds = slurm.calls('sbatch');
            out = testCase.runWrap(cmds{1}, 1000, 1, logDir);
            testCase.verifySubstring(out, 'STARTED');
            testCase.verifyFalse(isfolder(D), out);
        end

        function poolOfAnEnclosingCallStays(testCase)
            % a call that ends (an inner one, as a stitching wrapper the
            % pipeline's master runs itself) leaves only the pools it made
            testCase.assumeLinux();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            c = testCase.compiledSetup(true);
            testCase.startKind(c, slurm);
            testCase.submit(c, 1, 'ok', -1);
            D = testCase.poolDir(c);
            mark = single_job_pools('mark');
            single_job_pools('end', mark);
            testCase.verifyTrue(isfolder([D, '/live']));
            testCase.submit(c, 2, 'ok', -1);
            testCase.verifyEqual(strip(fileread([D, '/runnable.txt'])), sprintf('1\n2'), 'The pool took no more tasks.');
        end
    end

    methods (Access = private)
        function assumeLinux(testCase)
            testCase.assumeTrue(isunix && ~ismac, 'Task runners run on Linux only.');
        end

        function [c] = compiledSetup(testCase, withRunners, varargin)
            % a stand-in build, its config file (configTaskRunner: a
            % taskRunner key) and a job log folder
            ip = inputParser;
            ip.addParameter('configTaskRunner', []);
            ip.parse(varargin{:});
            c.root = testCase.createTempFolder();
            c.b = pktest.fakeRunnerBuild([c.root, 'build'], withRunners);
            c.logDir = [c.root, 'logs/'];
            mkdir(c.logDir);
            mkdir([c.root, 'out']);
            c.outFns = arrayfun(@(x) sprintf('%sout/o_%d.txt', c.root, x), 1 : 10, 'unif', 0);
            s = struct('MCRParam', c.b.mcr, 'MCCMasterStr', c.b.master, 'MCRCacheRoot', [c.root, 'cache']);
            if ~isempty(ip.Results.configTaskRunner)
                s.taskRunner = ip.Results.configTaskRunner;
            end
            c.config = [c.root, 'config.json'];
            testCase.writeLines(c.config, {jsonencode(s)});
            % (the stand-in build has no parsers manifest)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture('mcc:unverifiedBuild'));
        end

        function [jobId, status, submitted] = submit(~, c, k, mode, jobId, varargin)
            % generic_single_job_submit_wrapper for compiled task k
            % (pk_task('<out k>','<mode>'), array task k) and its job jobId
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.addParameter('funcName', 'pk_task');
            ip.addParameter('maxTrialNum', 1);
            ip.addParameter('jobLogDir', c.logDir);
            ip.addParameter('jobLogFname', [c.logDir, 'job_%A_%a.out']);
            ip.parse(varargin{:});
            args = [fieldnames(ip.Unmatched), struct2cell(ip.Unmatched)]';
            funcStr = sprintf('%s(''%s'',''%s'')', ip.Results.funcName, c.outFns{k}, mode); %#ok<NASGU>
            [~, jobId, status, submitted] = evalc(['generic_single_job_submit_wrapper(funcStr, jobId, k, ', ...
                '''jobLogDir'', ip.Results.jobLogDir, ''jobLogFname'', ip.Results.jobLogFname, ', ...
                '''jobErrorFname'', [c.logDir, ''job_%A_%a.err''], ', ...
                '''maxTrialNum'', ip.Results.maxTrialNum, ''mccMode'', true, ''ConfigFile'', c.config, args{:})']);
        end

        function [] = startKind(testCase, c, slurm, varargin)
            % the first task of the kind (task 10, with the submission
            % options varargin): it runs as before, without a pool; then the
            % fake Slurm forgets it, so that the test's tasks 1, 2, ... are
            % the pool's tasks 1, 2, ..., with jobs 1000, 1001, ...
            [~, ~, submitted] = testCase.submit(c, 10, 'ok', -1, varargin{:});
            testCase.assertEqual(submitted, 1);
            testCase.assertEmpty(dir([c.logDir, 'task_runner_single_*']), 'The first task of a kind made a pool.');
            slurm.clearCalls();
        end

        function [] = startMatlabKind(testCase, logDir, slurm, varargin)
            % the same for MATLAB tasks (a pktest.frameworkTask string)
            funcStr = sprintf('pktest.frameworkTask(''%sfirst_in.txt'',''%sfirst_out.txt'',''%sfirst.log'',''copy'')', ...
                logDir, logDir, logDir);
            testCase.submitMatlab(funcStr, logDir, 99, varargin{:});
            testCase.assertNumElements(slurm.calls('sbatch'), 1);
            testCase.assertEmpty(dir([logDir, 'task_runner_single_*']), 'The first task of a kind made a pool.');
            slurm.clearCalls();
        end

        function [] = submitMatlab(~, funcStr, logDir, k, varargin) %#ok<INUSL>
            % generic_single_job_submit_wrapper for MATLAB task string funcStr
            evalc(['generic_single_job_submit_wrapper(funcStr, -1, k, ''jobLogDir'', logDir, ', ...
                '''jobLogFname'', [logDir, ''job_%A_%a.out''], ''jobErrorFname'', [logDir, ''job_%A_%a.err''], ', ...
                '''maxTrialNum'', 1, varargin{:})']);
        end

        function [out] = runJob(testCase, c, slurm, k, doneFn, mayFail)
            % runs the k-th submitted job's command (array job 999 + k) in
            % bash, as Slurm would; with doneFn in the background, which it
            % makes when the job ends; mayFail: the command may fail
            if nargin < 5
                doneFn = '';
            end
            if nargin < 6
                mayFail = false;
            end
            cmds = slurm.calls('sbatch');
            arr = regexp(cmds{k}, '--array=(\d+)', 'tokens', 'once');
            out = testCase.runWrap(cmds{k}, 999 + k, str2double(arr{1}), c.root, doneFn, mayFail);
        end

        function [out] = runWrap(testCase, cmd, jobId, arrayId, folder, doneFn, mayFail)
            % runs a job's --wrap command in bash with Slurm's job variables
            if nargin < 6
                doneFn = '';
            end
            if nargin < 7
                mayFail = false;
            end
            jobFn = sprintf('%sjob_%d_%d.sh', folder, jobId, arrayId);
            body = {extractAfter(cmd, '--wrap=')};
            if ~isempty(doneFn)
                body{end + 1} = sprintf('touch "%s"', doneFn);
            end
            testCase.writeLines(jobFn, body);
            envStr = sprintf(['env -u PETAKIT5D_RUNNER SLURM_JOB_ID=%d SLURM_ARRAY_JOB_ID=%d SLURM_ARRAY_TASK_ID=%d ', ...
                'SLURMD_NODENAME=testnode SLURM_JOB_END_TIME= TMPDIR="%s"'], jobId, jobId, arrayId, tempdir);
            if ~isempty(doneFn)
                system(sprintf('%s bash "%s" > "%s.log" 2>&1 &', envStr, jobFn, jobFn));
                out = '';
                return;
            end
            [status, out] = system(sprintf('%s bash "%s" 2>&1', envStr, jobFn));
            if ~mayFail
                testCase.assertEqual(status, 0, out);
            end
        end

        function [] = runMatlabRunner(testCase, logDir, f)
            % the MATLAB part of job f of the (only) MATLAB pool in this
            % MATLAB: the job script claims task f (writing its task string
            % for the runner), then task_runner('matlab') runs it and the
            % next ones
            pool = dir([logDir, 'task_runner_single_*']);
            testCase.assertNumElements(pool, 1);
            D = fullfile(pool.folder, pool.name);
            L = [testCase.createTempFolder(), 'local'];
            mkdir(L);
            st = [L, '/s1_1'];
            [status, out] = system(sprintf('bash "%s/run_tasks.sh" claim %d 0 1 0 "%s" "%s" "%s/next_1"', D, f, st, L, L));
            testCase.assertEqual(status, 0, out);
            testCase.assertEqual(strtrim(fileread([L, '/next_1'])), sprintf('%d', f));
            testCase.writeLines(st, {sprintf('c %d', f)});
            setenv('PETAKIT5D_RUNNER_MEM_GB', '100');
            evalc('task_runner(''matlab'', D, sprintf(''%d'', f), 3, st, L, f)');
            testCase.verifySubstring(fileread(st), 'x none');
        end

        function [funcStrs, outFns, logDir, logFn, inFns] = matlabTasks(testCase, n)
            % n pktest.frameworkTask strings (copy an input) and a job log folder
            dataPath = testCase.createTempFolder();
            mkdir([dataPath, 'out']);
            logDir = [dataPath, 'logs/'];
            mkdir(logDir);
            logFn = [dataPath, 'calls.log'];
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : n, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout/out_%d.txt', dataPath, x), 1 : n, 'unif', 0);
            for f = 1 : n
                testCase.writeLines(inFns{f}, {sprintf('task %d', f)});
            end
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''copy'')', x, y, logFn), ...
                inFns, outFns, 'unif', 0);
        end

        function [D] = poolDir(testCase, c)
            % the (only) pool in the job log folder of setup c (or folder c)
            logDir = c;
            if isstruct(c)
                logDir = c.logDir;
            end
            pool = dir([logDir, 'task_runner_single_*']);
            testCase.assertNumElements(pool, 1);
            D = fullfile(pool.folder, pool.name);
        end

        function [owner, marked, failed] = claimOf(~, D, i)
            % the owner of task i's claim, and its resubmitted and failed marks
            [~, owner] = system(sprintf('readlink "%s/claims/t%d"', D, i));
            owner = strtrim(owner);
            marked = isfile(sprintf('%s/claims/r%d', D, i));
            [status, ~] = system(sprintf('[ -L "%s/claims/f%d" ]', D, i));
            failed = status == 0;
        end

        function [runs] = taskRuns(~, c)
            % the stand-in tasks that ran: task (from its output name),
            % mode, process, runner
            runs = struct('task', {}, 'mode', {}, 'pid', {}, 'runner', {});
            if ~isfile(c.b.log)
                return;
            end
            lines = splitlines(strip(fileread(c.b.log)));
            for i = 1 : numel(lines)
                p = split(lines{i});
                k = regexp(p{1}, 'o_(\d+)\.txt$', 'tokens', 'once');
                runs(i) = struct('task', str2double(k{1}), 'mode', p{2}, 'pid', str2double(p{4}), ...
                    'runner', str2double(p{6}));
            end
        end

        function [lines] = readLog(~, logFn)
            % output paths of all pktest.frameworkTask calls, in call order
            lines = cell(0, 1);
            if exist(logFn, 'file')
                lines = splitlines(strip(fileread(logFn)));
                lines = lines(~cellfun(@isempty, lines));
            end
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
