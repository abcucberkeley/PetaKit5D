classdef TaskMetricsTest < pktest.TestCase
    % task metrics of the job framework (task_metrics, read_task_metrics):
    % with taskMetrics, every task run leaves a record of its resources in
    % <jobLogDir>/task_metrics/<run>/ (local loop, MATLAB jobs, compiled
    % tasks, local pool, single jobs), and without it nothing changes. Slurm
    % is faked (pktest.FakeSlurmFixture): the MATLAB command of a job is run
    % in this MATLAB, a compiled task line with bash and a fake compiled master
    % that logs the metrics environment it gets.

    properties (Constant)
        EnvNames = {'PETAKIT5D_TASK_METRICS', 'PETAKIT5D_TASK', 'PETAKIT5D_TASK_GPU', ...
            'SLURM_ARRAY_JOB_ID', 'SLURM_ARRAY_TASK_ID'}
    end

    methods (TestMethodSetup)
        function clearCaches(testCase)
            % config files are cached in persistent variables; task_metrics
            % keeps this process's file name and the current task,
            % task_metrics_run the single-job folders
            clear generic_computing_frameworks_wrapper generic_single_job_submit_wrapper task_metrics task_metrics_run
        end

        function restoreMetricsEnvironment(testCase)
            for i = 1 : numel(testCase.EnvNames)
                testCase.addTeardown(@setenv, testCase.EnvNames{i}, getenv(testCase.EnvNames{i}));
                setenv(testCase.EnvNames{i}, '');
            end
        end
    end

    methods (Test, TestTags = {'Unit'})
        function localTasksAreRecorded(testCase)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(3);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');
            t0 = posixtime(datetime('now', 'TimeZone', 'local'));

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, 'parseCluster', false, ...
                'masterCompute', true, 'taskMetrics', true);

            t1 = posixtime(datetime('now', 'TimeZone', 'local'));
            testCase.verifyEqual(flags, true(3, 1));
            testCase.verifyEqual(cellfun(@fileread, outFns, 'unif', 0), cellfun(@fileread, inFns, 'unif', 0));
            % without a job log folder the metrics go to <data>/job_logs
            [t, runs] = read_task_metrics([dataPath, 'job_logs']);
            testCase.assertEqual(height(t), 3);
            t = sortrows(t, 'task');
            testCase.verifyEqual(t.task', 1 : 3);
            testCase.verifyEqual(t.func, repmat("pktest.frameworkTask", 3, 1));
            testCase.verifyEqual(t.mode, repmat("local", 3, 1));
            testCase.verifyEqual(t.status, repmat("ok", 3, 1));
            testCase.verifyEqual(t.exitCode, zeros(3, 1));
            testCase.verifyEqual(t.output, string(outFns(:)));
            testCase.verifyEqual(t.input, string(inFns(:)));
            testCase.verifyEqual(t.pid, repmat(feature('getpid'), 3, 1));
            testCase.verifyTrue(all(t.node ~= ""), 'No node name.');
            testCase.verifyTrue(all(t.start >= t0 - 1 & t.end <= t1 + 1 & t.end >= t.start), ...
                'The task times are not seconds since 1970 within the call.');
            testCase.verifyTrue(all(abs(t.end - t.start - t.wallSec) < 1e-3));
            testCase.verifyTrue(all(t.cpuSec >= 0 & t.threads == maxNumCompThreads));
            testCase.verifyTrue(all(t.metricsSec < 0.5), 'The metrics took longer than expected.');
            testCase.assertEqual(height(runs), 1);
            testCase.verifyEqual([runs.nTasks, runs.nDone, runs.nRecords], [3, 3, 3]);
            testCase.verifyEqual(runs.backend, "matlab");
            testCase.verifyEqual(runs.func, "pktest.frameworkTask");
            testCase.verifyEqual(runs.parseCluster, 0);
            % the master's Slurm job (if any) tells the tasks it runs itself
            rec = jsondecode(strtok(fileread([char(runs.folder), 'run.ndjson']), newline));
            testCase.verifyEqual(string(rec.jobId), string(getenv('SLURM_JOB_ID')));
            % strict JSON for other readers: null, not NaN, for missing values
            taskFiles = dir([char(runs.folder), 'task_*.ndjson']);
            testCase.verifyFalse(contains(fileread([char(runs.folder), taskFiles(1).name]), 'NaN'));
        end

        function metricsAreOffByDefault(testCase)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');

            flags = testCase.verifyWarningFree(@() generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true));

            testCase.verifyEqual(flags, true(2, 1));
            testCase.verifyFalse(exist([dataPath, 'job_logs'], 'dir') == 7, 'A metrics folder was made.');
        end

        function configFileTurnsMetricsOn(testCase)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');
            configFn = [dataPath, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(struct('taskMetrics', true)));
            fclose(fid);

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', false, ', ...
                '''masterCompute'', true, ''ConfigFile'', configFn)']);

            t = read_task_metrics([dataPath, 'job_logs']);
            testCase.verifyEqual(sort(t.task)', 1 : 2);
        end

        function failedLocalTaskIsRecordedAndStillFails(testCase)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(1);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'error');

            testCase.verifyError(@() generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'masterCompute', true, 'taskMetrics', true), 'pktest:frameworkTask:failed');

            t = read_task_metrics([dataPath, 'job_logs']);
            testCase.assertEqual(height(t), 1);
            testCase.verifyEqual(t.status, "error");
            testCase.verifyEqual(t.exitCode, 1);
            testCase.verifySubstring(char(t.error), 'pktest:frameworkTask:failed');
            % the task ran in try and again in catch, in one record
            testCase.verifyNumElements(testCase.readLog(logFn), 2);
        end

        function bashTaskExitStatusIsRecorded(testCase)
            [inFns, outFns, ~, dataPath] = testCase.makeTasks(1);

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, {''exit 3''}, ''parseCluster'', false, ', ...
                '''masterCompute'', true, ''language'', ''bash'', ''maxTrialNum'', 1, ''unitWaitTime'', 0, ''taskMetrics'', true)']);

            t = read_task_metrics([dataPath, 'job_logs']);
            testCase.assertEqual(height(t), 1);
            testCase.verifyEqual([t.status, t.func], ["failed", "exit 3"]);
            testCase.verifyEqual(t.exitCode, 3);
        end

        function memoryAndCpuOnLinux(testCase)
            % the peak resident memory is the task's own: it is reset at the
            % start of each task (a task that allocates 512 MB, then one that
            % doesn't); user and system CPU add up to the CPU time (a task
            % that computes, mostly user time)
            testCase.assumeTrue(isunix && ~ismac, 'Resource use beyond time is recorded on Linux only.');
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(3);
            funcStrs = [{'ones(2^26,1)'}, testCase.taskStrings(inFns(2), outFns(2), logFn, 'copy'), ...
                {'sum(sqrt(rand(1,4e7)))'}];

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', false, ', ...
                '''masterCompute'', true, ''maxTrialNum'', 1, ''unitWaitTime'', 0, ''taskMetrics'', true)']);

            t = sortrows(read_task_metrics([dataPath, 'job_logs']), 'task');
            testCase.assertEqual(height(t), 3);
            testCase.verifyGreaterThan(t.rssPeakGB(1) - t.rssStartGB(1), 0.45, 'The 512 MB array is not in the peak.');
            testCase.verifyLessThan(t.rssPeakGB(2) - t.rssStartGB(2), 0.2, 'The peak was not reset for the second task.');
            testCase.verifyLessThan(abs(t.userSec + t.sysSec - t.cpuSec), 0.03 + 0.05 * t.cpuSec);
            testCase.verifyGreaterThan(t.userSec(3), 0.5 * t.cpuSec(3), 'The computing task should be mostly user time.');
            % (transparent huge pages fault 2 MB at a time)
            testCase.verifyGreaterThan(t.minorFaults(1), 200, 'The page faults of filling 512 MB are missing.');
            testCase.verifyTrue(all(t.cpus >= 1 & t.osThreads >= 1 & t.memTotalGB > 0));
            testCase.verifyTrue(all(t.processStart < t.start));
        end

        function unwritableMetricsDoNotFailTasks(testCase)
            % a file where the metrics folder would be: no metrics, but the
            % tasks run as usual
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');
            fclose(fopen([dataPath, 'job_logs'], 'w'));

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, 'parseCluster', false, ...
                'masterCompute', true, 'taskMetrics', true);

            testCase.verifyEqual(flags, true(2, 1));
            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort(outFns(:)));
        end

        function failedRecordWriteOnlyWarns(testCase)
            dataPath = testCase.createTempFolder();
            fclose(fopen([dataPath, 'notAFolder'], 'w'));
            tm = task_metrics('start', 'f(1)', [dataPath, 'notAFolder/run/'], 1, 'local', false);

            testCase.verifyWarning(@() task_metrics('stop', tm), 'task_metrics:failed');
        end

        function matlabJobsRunEachTaskWithMetrics(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''taskBatchNum'', 2, ''jobLogDir'', dataPath, ', ...
                '''taskMetrics'', true)']);

            sbatch = slurm.calls('sbatch');
            testCase.assertNumElements(sbatch, 1);
            % the job's MATLAB command, run here
            jobCmd = regexp(sbatch{1}, '-batch "setup\(\[\]\);(t0_=tic;.*?;toc\(t0_\))"', 'tokens', 'once');
            testCase.assertNotEmpty(jobCmd, 'No MATLAB command in the job.');
            evalc(jobCmd{1});
            testCase.verifyEqual(cellfun(@fileread, outFns, 'unif', 0), cellfun(@fileread, inFns, 'unif', 0));
            [t, runs] = read_task_metrics([dataPath, 'task_metrics']);
            testCase.assertEqual(height(t), 2);
            testCase.verifyEqual(sort(t.task)', 1 : 2);
            testCase.verifyEqual(t.mode, repmat("slurm_matlab", 2, 1));
            testCase.verifyEqual(t.status, repmat("ok", 2, 1));
            testCase.verifyEqual(sort(t.output), sort(string(outFns(:))));
            testCase.verifyEqual(runs.backend, "matlab");
            testCase.verifyEqual(runs.parseCluster, 1);
        end

        function masterComputedTasksAreRecordedAsMaster(testCase)
            % in cluster mode with masterCompute the master runs the last
            % batch itself (here task 2; the job of task 1 writes its output)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');
            slurm.setJobStates('1000_1', {'R', 'done'}, outFns(1));

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', true, ''maxTrialNum'', 1, ''queryInterval'', 0, ''jobLogDir'', dataPath, ', ...
                '''taskMetrics'', true)']);

            t = read_task_metrics([dataPath, 'task_metrics']);
            testCase.assertEqual(height(t), 1);
            testCase.verifyEqual([t.task, t.batchSize], [2, 1]);
            testCase.verifyEqual(t.mode, "master");
        end

        function matlabJobStatementRunsWithoutTaskMetrics(testCase)
            % a job whose MATLAB can't find task_metrics still runs its task
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(1);
            funcStr = testCase.taskStrings(inFns, outFns, logFn, 'copy');
            str = task_metrics_run('matlab_task', [dataPath, 'metrics/'], 1, funcStr{1}, false);

            eval(strrep(str, 'task_metrics(', 'no_such_task_metrics_function('));

            testCase.verifyEqual(fileread(outFns{1}), fileread(inFns{1}));
        end

        function matlabJobStatementRethrowsTaskErrors(testCase)
            % a failed task fails its job as before, after its record
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(1);
            funcStr = testCase.taskStrings(inFns, outFns, logFn, 'error');
            runDir = [dataPath, 'metrics/'];
            mkdir(runDir);
            str = task_metrics_run('matlab_task', runDir, 4, funcStr{1}, false);

            % (in a function workspace, as the statement runs in a job's base workspace)
            testCase.verifyError(@() eval_statement(str), 'pktest:frameworkTask:failed');

            t = read_task_metrics(runDir);
            testCase.verifyEqual([t.status, t.mode], ["error", "slurm_matlab"]);
            testCase.verifyEqual(t.task, 4);
        end

        function matlabJobCommandIsUnchangedWithoutMetrics(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'copy');

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''taskBatchNum'', 2, ''jobLogDir'', dataPath)']);

            sbatch = slurm.calls('sbatch');
            testCase.assertNumElements(sbatch, 1);
            testCase.verifySubstring(sbatch{1}, sprintf('-batch "setup([]);t0_=tic;%s;%s;toc(t0_)"', funcStrs{:}));
            testCase.verifyFalse(contains(sbatch{1}, 'task_metrics'));
            testCase.verifyFalse(exist([dataPath, 'task_metrics'], 'dir') == 7);
        end

        function compiledTaskLinesCarryMetricsEnvironment(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();

            [lines, logFn] = testCase.compiledTaskLines(dataPath, 2, {'taskMetrics', true});

            testCase.assertNumElements(lines, 2);
            runDir = testCase.runFolder(dataPath);
            for f = 1 : 2
                testCase.runTaskLine(lines{f}, 'SLURM_JOB_ID=42 PARALLEL_JOBSLOT=1');
            end
            envs = splitlines(strip(fileread(logFn)));
            testCase.verifyEqual(envs, {sprintf('%s|1||XR_crop_frame', runDir); sprintf('%s|2||XR_crop_frame', runDir)});
        end

        function gpuTaskLinesAskForTheGPU(testCase)
            % a GPU task records its GPU (the flag is only in GPU task lines)
            testCase.assumeFalse(ispc, 'The task lines of Windows are batch files.');
            testCase.verifyEqual(task_metrics_run('env', '/m/run/', 3, true), ...
                'PETAKIT5D_TASK_METRICS="/m/run/" PETAKIT5D_TASK=3 PETAKIT5D_TASK_GPU=1 ');
            testCase.verifyEqual(task_metrics_run('env', '/m/run/', 3, false), 'PETAKIT5D_TASK_METRICS="/m/run/" PETAKIT5D_TASK=3 ');
        end

        function compiledLocalRunKeepsAJoblog(testCase)
            % without a cluster the master runs the task lines with GNU
            % parallel (here a stand-in that logs its arguments) and keeps a
            % joblog in the metrics folder
            testCase.assumeFalse(ispc, 'Needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            binDir = [dataPath, 'bin/'];
            mkdir(binDir);
            fid = fopen([binDir, 'parallel'], 'w');
            fprintf(fid, '#!/bin/sh\necho "$*" >> "%sparallel.log"\n', dataPath);
            fclose(fid);
            fileattrib([binDir, 'parallel'], '+x');
            path0 = getenv('PATH');
            testCase.addTeardown(@setenv, 'PATH', path0);
            setenv('PATH', [binDir, pathsep, path0]);
            [~, ~, configFn] = testCase.fakeCompiledMaster(dataPath, {'wholeNodeJob', true, 'GNUparallel', true, ...
                'maxCPUNum', 4, 'taskMetrics', true});
            inFns = {[dataPath, 'in_1.txt'], [dataPath, 'in_2.txt']};
            cellfun(@(x) fclose(fopen(x, 'w')), inFns);
            outFns = {[dataPath, 'out_1.txt'], [dataPath, 'out_2.txt']};
            funcStrs = cellfun(@(x, y) sprintf('XR_crop_frame(''%s'',''%s'',[1,1,1,2,2,2])', x, y), inFns, outFns, 'unif', 0);

            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', false, ', ...
                '''mccMode'', true, ''maxTrialNum'', 1, ''jobLogDir'', dataPath, ''ConfigFile'', configFn)']);

            args = strip(fileread([dataPath, 'parallel.log']));
            testCase.verifySubstring(args, sprintf('--joblog +%sparallel_master.joblog --retries 2 --jobs 4 ', ...
                testCase.runFolder(dataPath)));
        end

        function compiledTaskLinesAreUnchangedWithoutMetrics(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();

            [lines, logFn] = testCase.compiledTaskLines(dataPath, 1, {});

            testCase.assertNumElements(lines, 1);
            testCase.verifyTrue(startsWith(lines{1}, sprintf(['OUT="%sout_1.txt"; if [ ! -f "$OUT" ] && [ ! -d "$OUT" ]; ', ...
                'then MCR_CACHE_ROOT=%scache/mcr_cache_'], dataPath, dataPath)), lines{1});
            testCase.verifyFalse(contains(lines{1}, 'PETAKIT5D'));
            testCase.runTaskLine(lines{1}, 'SLURM_JOB_ID=42');
            testCase.verifyEqual(strip(fileread(logFn)), '|||XR_crop_frame');
        end

        function compiledJobsKeepAJoblog(testCase)
            % GNU parallel records each task line's start, run time and exit
            % status in the run's metrics folder
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            parallelConfig = {'wholeNodeJob', true, 'GNUparallel', true, 'maxCPUNum', 4};

            testCase.compiledTaskLines(dataPath, 2, [parallelConfig, {'taskMetrics', true}]);

            sbatch = slurm.calls('sbatch');
            testCase.assertNumElements(sbatch, 1);
            testCase.verifySubstring(sbatch{1}, sprintf(['parallel --joblog "%sparallel_${SLURM_JOB_ID}_', ...
                '${SLURM_ARRAY_TASK_ID}.joblog" --ungroup --retries 2 --jobs 4 '], testCase.runFolder(dataPath)));

            dataPath = testCase.createTempFolder();
            testCase.compiledTaskLines(dataPath, 2, parallelConfig);
            sbatch = slurm.calls('sbatch');
            testCase.verifySubstring(sbatch{end}, 'parallel --ungroup --retries 2 --jobs 4 ');
            testCase.verifyFalse(contains(sbatch{end}, 'joblog'));
        end

        function mccMasterRecordsTaskFromEnvironment(testCase)
            % a compiled task gets its metrics settings from its task line's
            % environment
            dataPath = testCase.createTempFolder();
            runDir = [dataPath, 'metrics/'];
            mkdir(runDir);
            inFn = [dataPath, 'x.tif'];
            writetiff(uint16(reshape(1 : 512, 8, 8, 8)), inFn);
            setenv('PETAKIT5D_TASK_METRICS', runDir);
            setenv('PETAKIT5D_TASK', '7');

            evalc('mccMaster(''XR_crop_frame'', inFn, [dataPath, ''y.tif''], ''[1,1,1,4,4,4]'')');

            testCase.verifyEqual(size(readtiff([dataPath, 'y.tif'])), [4, 4, 4]);
            % a task that fails
            setenv('PETAKIT5D_TASK', '8');
            testCase.verifyError(@() evalc('mccMaster(''no_such_task_function'', ''a'')'), ?MException);

            t = sortrows(read_task_metrics(runDir), 'task');
            testCase.assertEqual(height(t), 2);
            testCase.verifyEqual(t.task', [7, 8]);
            testCase.verifyEqual(t.func, ["XR_crop_frame"; "no_such_task_function"]);
            mode = "mcc";
            if ~isempty(getenv('SLURM_JOB_ID'))
                mode = "slurm_mcc";
            end
            testCase.verifyEqual(t.mode, [mode; mode]);
            testCase.verifyEqual(t.status, ["ok"; "error"]);
            testCase.verifySubstring(char(t.error(2)), 'no_such_task_function');
        end

        function callWithoutMetricsClearsTaskEnvironment(testCase)
            % a compiled task with metrics runs a framework call without
            % metrics: the tasks of that call (which get this environment,
            % e.g., through sbatch) must not record into the task's run
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            setenv('PETAKIT5D_TASK_METRICS', [dataPath, 'outer_run/']);
            setenv('PETAKIT5D_TASK', '5');

            [lines, logFn] = testCase.compiledTaskLines(dataPath, 1, {});

            testCase.verifyEqual(getenv('PETAKIT5D_TASK_METRICS'), '');
            testCase.verifyEqual(getenv('PETAKIT5D_TASK'), '');
            testCase.runTaskLine(lines{1}, 'SLURM_JOB_ID=42');
            testCase.verifyEqual(strip(fileread(logFn)), '|||XR_crop_frame');
        end

        function singleJobWithoutMetricsClearsTaskEnvironment(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            setenv('PETAKIT5D_TASK_METRICS', [dataPath, 'outer_run/']);
            setenv('PETAKIT5D_TASK', '5');

            evalc('generic_single_job_submit_wrapper(''disp(1)'', -1, 7, ''maxTrialNum'', 1, ''jobLogDir'', dataPath)');

            testCase.verifyEqual(getenv('PETAKIT5D_TASK_METRICS'), '');
            testCase.verifyFalse(exist([dataPath, 'task_metrics'], 'dir') == 7);
        end

        function readerLabelsTasksTheMasterRan(testCase)
            % compiled tasks in the master's own Slurm job ran in the master
            runDir = [testCase.createTempFolder(), 'run_2/'];
            mkdir(runDir);
            testCase.writeLines([runDir, 'run.ndjson'], {jsonencode(struct('kind', 'run_start', 'runId', 'run_2', ...
                'time', 100, 'jobId', '77', 'parseCluster', true, 'outputFullpaths', {{'/d/a', '/d/b'}}))});
            testCase.writeLines([runDir, 'task_n1_1_1.ndjson'], { ...
                jsonencode(struct('kind', 'task', 'runId', 'run_2', 'task', 1, 'mode', 'slurm_mcc', 'jobId', '77')), ...
                jsonencode(struct('kind', 'task', 'runId', 'run_2', 'task', 2, 'mode', 'slurm_mcc', 'jobId', '78'))});

            t = sortrows(read_task_metrics(runDir), 'task');

            testCase.verifyEqual(t.mode, ["master"; "slurm_mcc"]);
        end

        function singleJobsRecordSubmissionAndTask(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(1);
            funcStr = testCase.taskStrings(inFns, outFns, logFn, 'copy');

            evalc(['generic_single_job_submit_wrapper(funcStr{1}, -1, 7, ''maxTrialNum'', 1, ''jobLogDir'', dataPath, ', ...
                '''cpusPerTask'', 2, ''memAllocate'', 30, ''taskMetrics'', true)']);

            sbatch = slurm.calls('sbatch');
            testCase.assertNumElements(sbatch, 1);
            % the job runs its command with Slurm's array job and task ids
            jobCmd = regexp(sbatch{1}, '-batch "setup\(\[\]\);(t0_=tic;.*?;toc\(t0_\))"', 'tokens', 'once');
            testCase.assertNotEmpty(jobCmd, 'No MATLAB command in the job.');
            setenv('SLURM_ARRAY_JOB_ID', '1000');
            setenv('SLURM_ARRAY_TASK_ID', '7');
            evalc(jobCmd{1});
            [t, runs] = read_task_metrics([dataPath, 'task_metrics']);
            testCase.assertEqual(height(t), 1);
            testCase.verifyEqual([t.mode, t.status, t.arrayJobId], ["slurm_matlab", "ok", "1000"]);
            testCase.verifyEqual(t.task, 7);
            testCase.verifyEqual([t.cpusPerTask, t.memAllocate], [2, 30], 'The submission was not joined.');
            testCase.assertEqual(height(runs), 1);
            testCase.verifyEqual([runs.kind, runs.func], ["single", "pktest.frameworkTask"]);
        end

        function singleCompiledJobCarriesMetricsEnvironment(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            dataPath = testCase.createTempFolder();
            [masterFn, logFn, configFn] = testCase.fakeCompiledMaster(dataPath, {});

            evalc(['generic_single_job_submit_wrapper(''XR_crop_frame(''''in'''',''''out'''',[1,1,1,2,2,2],''''pad'''',false)'', ', ...
                '-1, 7, ''maxTrialNum'', 1, ''mccMode'', true, ''jobLogDir'', dataPath, ''ConfigFile'', configFn, ''taskMetrics'', true)']);

            testCase.verifyTrue(isfile(masterFn));
            calls = slurm.calls('sbatch');
            testCase.assertNumElements(calls, 1);
            % run the job's command as Slurm would, in the job's shell
            jobFn = [dataPath, 'job.sh'];
            fid = fopen(jobFn, 'w');
            fprintf(fid, '%s\n', extractAfter(calls{1}, '--wrap='));
            fclose(fid);
            [status, out] = system(sprintf('env -u PETAKIT5D_TASK_METRICS SLURM_JOB_ID=42 bash %s', jobFn));
            testCase.assertEqual(status, 0, out);
            testCase.verifyEqual(strip(fileread(logFn)), sprintf('%s|7||XR_crop_frame', testCase.runFolder(dataPath)));
        end

        function nestedRunRecordsItsParent(testCase)
            % a task that runs framework calls of its own (large-file tasks
            % do, one after another): the inner calls record the outer task
            % as their parent
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(1);
            innerOut = [dataPath, 'out/inner.txt'];
            funcStr = sprintf('pktest.nestedFrameworkTask(''%s'',''%s'',''%s'')', inFns{1}, innerOut, logFn);

            evalc(['generic_computing_frameworks_wrapper(inFns, {innerOut}, {funcStr}, ''parseCluster'', false, ', ...
                '''masterCompute'', true, ''taskMetrics'', true)']);

            [t, runs] = read_task_metrics([dataPath, 'job_logs']);
            testCase.assertEqual(height(runs), 3);
            outer = runs(runs.func == "pktest.nestedFrameworkTask", :);
            inner = runs(runs.func == "pktest.frameworkTask", :);
            testCase.assertEqual([height(outer), height(inner)], [1, 2]);
            testCase.verifyEqual(outer.parent, "");
            testCase.verifyEqual(inner.parent, repmat(outer.runId + "/1", 2, 1));
            testCase.verifyEqual(sort(t.func), ["pktest.frameworkTask"; "pktest.frameworkTask"; "pktest.nestedFrameworkTask"]);
            testCase.verifyEqual(t.parent(t.func == "pktest.frameworkTask"), repmat(outer.runId + "/1", 2, 1));
            testCase.verifyEqual(t.parent(t.func == "pktest.nestedFrameworkTask"), "");
        end

        function readerJoinsJoblogAndRunSettings(testCase)
            % a compiled run: GNU parallel's joblog gives the time from the
            % end of a task to the end of its process, and lines without a
            % record (the process was killed, or the output existed)
            runDir = [testCase.createTempFolder(), 'run_1/'];
            mkdir(runDir);
            testCase.writeLines([runDir, 'run.ndjson'], {jsonencode(struct('kind', 'run_start', 'runId', 'run_1', ...
                'time', 100, 'elapsedSec', 1, 'backend', 'mcc', 'cpusPerTask', 2, 'memAllocate', 5, 'mccMode', true, ...
                'inputFullpaths', {{'/d/in_1.tif', '/d/in_2.tif', '/d/in_3.tif'}}, ...
                'outputFullpaths', {{'/d/out_1.tif', '/d/out_2.tif', '/d/out_3.tif'}}, ...
                'functionStrs', {{'f(1)', 'f(2)', 'f(3)'}})), ...
                jsonencode(struct('kind', 'run_end', 'time', 140, 'elapsedSec', 41, 'nTasks', 3, 'nDone', 2))});
            testCase.writeLines([runDir, 'task_n1_11_5.ndjson'], {jsonencode(struct('kind', 'task', 'runId', 'run_1', ...
                'task', 1, 'func', 'f', 'mode', 'slurm_mcc', 'status', 'ok', 'node', 'n1', 'pid', 11, 'jobId', '55', ...
                'start', 110, 'end', 112, 'wallSec', 2, 'processStart', 106.5))});
            testCase.writeLines([runDir, 'parallel_55_1.joblog'], { ...
                sprintf('Seq\tHost\tStarttime\tJobRuntime\tSend\tReceive\tExitval\tSignal\tCommand'), ...
                sprintf('1\t:\t106.400\t6.100\t0\t0\t0\t0\tOUT="/d/out_1.tif"; PETAKIT5D_TASK=1 run'), ...
                sprintf('2\t:\t106.600\t20.000\t0\t0\t137\t9\tOUT="/d/out_2.tif"; PETAKIT5D_TASK=2 run'), ...
                sprintf('3\t:\t106.800\t0.010\t0\t0\t0\t0\tOUT="/d/out_3.tif"; PETAKIT5D_TASK=3 run')});

            [t, runs] = read_task_metrics(runDir);

            t = sortrows(t, 'task');
            testCase.assertEqual(height(t), 3);
            testCase.verifyEqual(t.status, ["ok"; "killed"; "skipped"]);
            testCase.verifyEqual(t.startupSec(1), 3.5, 'AbsTol', 1e-9);
            testCase.verifyEqual(t.lineSec(1), 6.1, 'AbsTol', 1e-9);
            testCase.verifyEqual(t.exitSec(1), 0.5, 'AbsTol', 1e-9);
            testCase.verifyEqual(t.exitCode(2), 137);
            testCase.verifyEqual(t.output, ["/d/out_1.tif"; "/d/out_2.tif"; "/d/out_3.tif"]);
            testCase.verifyEqual(t.input(1), "/d/in_1.tif");
            testCase.verifyEqual([t.cpusPerTask(1), t.memAllocate(1)], [2, 5]);
            testCase.verifyEqual([runs.start, runs.end, runs.elapsedSec, runs.nDone], [99, 140, 41, 2]);
        end
    end

    methods (Test, TestTags = {'Integration', 'Slow'})
        function poolTasksAreRecordedOnTheWorkers(testCase)
            testCase.assumeTrue(license('test', 'Distrib_Computing_Toolbox') && ~isempty(which('gcp')), ...
                'Needs the Parallel Computing Toolbox.');
            p = gcp('nocreate');
            if ~isempty(p) && ~isa(p, 'parallel.ProcessPool')
                delete(p);
            end
            if isempty(gcp('nocreate'))
                % the framework starts it, with job files of its own (MATLAB
                % sessions starting pools in a shared home can clash)
                testCase.addTeardown(@() delete(gcp('nocreate')));
            end
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(4);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'where');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, 'parseCluster', false, ...
                'localWorkerNum', 2, 'taskMetrics', true);

            testCase.verifyEqual(flags, true(4, 1));
            [t, runs] = read_task_metrics([dataPath, 'job_logs']);
            testCase.assertEqual(height(t), 4);
            testCase.verifyEqual(sort(t.task)', 1 : 4);
            testCase.verifyEqual(t.mode, repmat("pool", 4, 1));
            testCase.verifyFalse(any(t.pid == feature('getpid')), 'The tasks did not run on the workers.');
            testCase.verifyEqual(sort(t.output), sort(string(outFns(:))));
            testCase.verifyEqual(runs.backend, "matlab");
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

        function [funcStrs] = taskStrings(~, inFns, outFns, logFn, mode)
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''%s'')', ...
                x, y, logFn, mode), inFns, outFns, 'unif', 0);
        end

        function [lines] = readLog(~, logFn)
            % output paths of all task calls, in call order
            lines = cell(0, 1);
            if exist(logFn, 'file')
                lines = splitlines(strip(fileread(logFn)));
                lines = lines(~cellfun(@isempty, lines));
            end
        end

        function [masterFn, logFn, configFn] = fakeCompiledMaster(~, dataPath, config)
            % a fake MATLAB Runtime and compiled master that logs the metrics
            % environment and the function of each task, and a config file
            mkdir([dataPath, 'mcr/bin']);
            logFn = [dataPath, 'envs.log'];
            masterFn = [dataPath, 'run_mccMaster.sh'];
            fid = fopen(masterFn, 'w');
            fprintf(fid, '#!/bin/bash\necho "$PETAKIT5D_TASK_METRICS|$PETAKIT5D_TASK|$PETAKIT5D_TASK_GPU|$2" >> "%s"\n', logFn);
            fclose(fid);
            fileattrib(masterFn, '+x');
            s = struct('MCRParam', [dataPath, 'mcr/'], 'MCCMasterStr', masterFn, ...
                'MCRCacheRoot', [dataPath, 'cache'], 'BashLaunchStr', 'echo');
            for i = 1 : 2 : numel(config)
                s.(config{i}) = config{i + 1};
            end
            configFn = [dataPath, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(s));
            fclose(fid);
        end

        function [lines, logFn] = compiledTaskLines(testCase, dataPath, nTasks, config)
            % the task lines of the framework in mcc mode on nTasks tasks (the
            % fake Slurm never runs the jobs, so the framework gives up after
            % one try and leaves the task file behind)
            [~, logFn, configFn] = testCase.fakeCompiledMaster(dataPath, config);
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : nTasks, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout_%d.txt', dataPath, x), 1 : nTasks, 'unif', 0);
            cellfun(@(x) fclose(fopen(x, 'w')), inFns);
            funcStrs = cellfun(@(x, y) sprintf('XR_crop_frame(''%s'',''%s'',[1,1,1,2,2,2],''pad'',false)', x, y), ...
                inFns, outFns, 'unif', 0);
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''mccMode'', true, ''maxTrialNum'', 1, ''jobLogDir'', dataPath, ', ...
                '''ConfigFile'', configFn)']);
            taskFiles = dir([dataPath, 'tmp/*/input_*.txt']);
            lines = {};
            for i = 1 : numel(taskFiles)
                lines = [lines; splitlines(strip(fileread(fullfile(taskFiles(i).folder, taskFiles(i).name))))]; %#ok<AGROW>
            end
        end

        function [runDir] = runFolder(testCase, dataPath)
            % the only run in <dataPath>/task_metrics
            dir_info = dir([dataPath, 'task_metrics']);
            dir_info = dir_info([dir_info.isdir] & ~startsWith({dir_info.name}, '.'));
            testCase.assertNumElements(dir_info, 1, 'There should be one metrics folder.');
            runDir = [dataPath, 'task_metrics/', dir_info.name, '/'];
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

        function writeLines(~, fn, lines)
            fid = fopen(fn, 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
        end
    end
end


function [] = eval_statement(str)
% runs a MATLAB job's statement in a workspace of its own

eval(str);

end
