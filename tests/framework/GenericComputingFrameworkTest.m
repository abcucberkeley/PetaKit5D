classdef GenericComputingFrameworkTest < pktest.TestCase
    % generic_computing_frameworks_wrapper in local mode (parseCluster =
    % false): tasks run sequentially in this MATLAB session. Tasks call
    % pktest.frameworkTask, which logs every call so runs can be counted.

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
            if ~exist(logFn, 'file')
                lines = cell(0, 1);
                return;
            end
            lines = splitlines(strip(fileread(logFn)));
            lines = lines(~cellfun(@isempty, lines));
        end
    end
end
