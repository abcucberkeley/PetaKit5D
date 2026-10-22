classdef LocalPoolFrameworkTest < pktest.TestCase
    % generic_computing_frameworks_wrapper with localWorkerNum in local mode
    % (parseCluster = false): the tasks run on a local MATLAB process pool
    % (matlab_pool_generic_computing_wrapper), and the framework retries the
    % unfinished ones in this MATLAB. Needs the Parallel Computing Toolbox;
    % a pool of 2 workers is started once for the class (~15 s).

    methods (TestClassSetup)
        function needsParallelComputingToolbox(testCase)
            testCase.assumeTrue(license('test', 'Distrib_Computing_Toolbox') && ~isempty(which('gcp')), ...
                'Needs the Parallel Computing Toolbox.');
            testCase.addTeardown(@() delete(gcp('nocreate')));
        end
    end

    methods (TestMethodSetup)
        function startProcessPool(testCase)
            p = gcp('nocreate');
            if ~isempty(p) && ~isa(p, 'parallel.ProcessPool')
                delete(p);
                p = [];
            end
            if isempty(p)
                parpool('Processes', 2);
            end
        end
    end

    methods (Test, TestTags = {'Integration', 'Slow'})
        function tasksRunOnTheLocalPool(testCase)
            [inFns, outFns, logFn] = testCase.makeTasks(4);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'where');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'localWorkerNum', 2);

            testCase.verifyEqual(flags, true(4, 1));
            testCase.verifyEqual(cellfun(@fileread, outFns, 'unif', 0), repmat({'worker'}, 1, 4), ...
                'The tasks should run on the pool workers.');
            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort(outFns(:)), 'Each task should run once.');
        end

        function failedTasksAreRetriedInThisMatlab(testCase)
            % the tasks fail on the workers; the framework runs them again
            [inFns, outFns, logFn] = testCase.makeTasks(3);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'clientonly');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'localWorkerNum', 2);

            testCase.verifyEqual(flags, true(3, 1));
            testCase.verifyEqual(cellfun(@fileread, outFns, 'unif', 0), repmat({'client'}, 1, 3));
            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort([outFns(:); outFns(:)]), ...
                'Each task should run once on a worker and once in this MATLAB.');
        end

        function failedTasksAreReportedUnfinished(testCase)
            [inFns, outFns, logFn] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'clientonly');

            flags = matlab_pool_generic_computing_wrapper(inFns, outFns, funcStrs, 'workerNum', 2);

            testCase.verifyEqual(flags, false(2, 1));
            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort(outFns(:)));
        end

        function tasksThatDoNotFitInMemoryRunInThisMatlab(testCase)
            % 8 GB of memory: half of it holds no worker with a 10 GB task
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 8));
            [inFns, outFns, logFn] = testCase.makeTasks(3);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'where');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'localWorkerNum', 2, 'memAllocate', 10);

            testCase.verifyEqual(flags, true(3, 1));
            testCase.verifyEqual(cellfun(@fileread, outFns, 'unif', 0), repmat({'client'}, 1, 3));
        end

        function existingFinalOutputSkipsAllTasks(testCase)
            [inFns, outFns, logFn, dataPath] = testCase.makeTasks(2);
            finalFn = [dataPath, 'final.txt'];
            fclose(fopen(finalFn, 'w'));
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'where');

            flags = matlab_pool_generic_computing_wrapper(inFns, outFns, funcStrs, 'workerNum', 2, ...
                'finalOutFullpath', finalFn);

            testCase.verifyEqual(flags, true(2, 1));
            testCase.verifyEmpty(testCase.readLog(logFn));
        end

        function tasksWithMissingInputsAreSkipped(testCase)
            % as in the framework, which reports them as done
            [inFns, outFns, logFn] = testCase.makeTasks(3);
            delete(inFns{2});
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'where');

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'localWorkerNum', 2);

            testCase.verifyEqual(flags, true(3, 1));
            testCase.verifyEqual(sort(testCase.readLog(logFn)), sort(outFns([1, 3])'));
            testCase.verifyFalse(exist(outFns{2}, 'file') == 2);
        end

        function threadPoolIsNotUsed(testCase)
            % thread-based workers cannot run mex functions: no task is sent
            % to them, the framework runs the tasks in this MATLAB
            delete(gcp('nocreate'));
            parpool('Threads', 2);
            [inFns, outFns, logFn] = testCase.makeTasks(2);
            funcStrs = testCase.taskStrings(inFns, outFns, logFn, 'where');

            out = evalc('flags = matlab_pool_generic_computing_wrapper(inFns, outFns, funcStrs, ''workerNum'', 2);');
            testCase.verifyEqual(flags, false(2, 1));
            testCase.verifyFalse(contains(out, 'local pool'), 'Tasks were sent to the thread pool.');
            testCase.verifyEmpty(testCase.readLog(logFn));

            flags = generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ...
                'parseCluster', false, 'localWorkerNum', 2);

            testCase.verifyEqual(flags, true(2, 1));
            testCase.verifyEqual(cellfun(@fileread, outFns, 'unif', 0), repmat({'client'}, 1, 2));
            testCase.verifyClass(gcp('nocreate'), 'parallel.ThreadPool', 'The thread pool was replaced.');
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
                lines = {};
                return;
            end
            lines = splitlines(strtrim(fileread(logFn)));
            lines = lines(~cellfun(@isempty, lines));
        end
    end
end
