classdef MccClusterWrapperTest < pktest.TestCase
    % mcc_slurm_cluster_generic_computing_wrapper, the job framework's
    % compiled mode, with a fake Slurm (pktest.FakeSlurmFixture) and a fake
    % compiled master that records its environment. The jobs never run;
    % the tests run the task lines the framework writes for them.

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
            % the framework in mcc mode on nTasks tasks, with a fake MATLAB
            % Runtime and a fake compiled master that logs MCR_CACHE_ROOT
            dataPath = testCase.createTempFolder();
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
