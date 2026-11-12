classdef CoordinatorNiceTest < pktest.TestCase
    % jobs that wait for jobs they submit (coordinating jobs) go after
    % those: coordinator_nice is the part of a job's priority that its age
    % and size can add (Slurm's multifactor priority), and the job
    % framework's nice option lowers the priority of its jobs by it (sbatch
    % --nice). Pending coordinating jobs were older than the jobs of the
    % running ones and outranked them.

    methods (Test, TestTags = {'Unit'})
        function niceIsTheWeightsOfAgeAndSize(testCase)
            testCase.assumeFalse(ispc, 'The scontrol stub needs a POSIX shell.');
            % priority type and weights in `scontrol show config`, nice
            cases = {'priority/multifactor', {'PriorityWeightAge       = 1000', 'PriorityWeightJobSize   = 0'}, 1000; ...
                'priority/multifactor', {'PriorityWeightAge = 1000', 'PriorityWeightJobSize = 500'}, 1500; ...
                'priority/basic', {'PriorityWeightAge = 1000'}, 0};
            for i = 1 : size(cases, 1)
                testCase.stubScontrol([{'PriorityDecayHalfLife   = 1-00:00:00', ...
                    sprintf('PriorityType            = %s', cases{i, 1})}, cases{i, 2}, ...
                    {'PriorityWeightFairShare = 100000'}], 0);
                testCase.verifyEqual(coordinator_nice(), cases{i, 3}, sprintf('Case %d.', i));
            end
            % (scontrol fails: no Slurm)
            testCase.stubScontrol({}, 1);
            testCase.verifyEqual(coordinator_nice(), 0);
        end

        function singleJobsGetTheirNice(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            logDir = testCase.createTempFolder(); %#ok<NASGU>
            evalc('generic_single_job_submit_wrapper(''disp(1)'', -1, 7, ''maxTrialNum'', 1, ''jobLogDir'', logDir, ''nice'', 1000)');
            evalc('generic_single_job_submit_wrapper(''disp(2)'', -1, 8, ''maxTrialNum'', 1, ''jobLogDir'', logDir)');
            calls = slurm.calls('sbatch');
            testCase.assertNumElements(calls, 2);
            testCase.verifySubstring(calls{1}, ' --nice=1000');
            testCase.verifyFalse(contains(calls{2}, '--nice'), 'A job without nice got one.');
        end

        function frameworkJobsGetTheirNice(testCase)
            % (the jobs don't run: the framework gives up after its try)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 100));
            testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', [], 100));
            folder = testCase.createTempFolder();
            inputs = {[folder, 'in_1.txt'], [folder, 'in_2.txt']};
            outputs = {[folder, 'out_1.txt'], [folder, 'out_2.txt']}; %#ok<NASGU>
            cellfun(@(x) fclose(fopen(x, 'w')), inputs);
            evalc(['generic_computing_frameworks_wrapper(inputs, outputs, {''disp(1)'', ''disp(2)''}, ', ...
                '''parseCluster'', true, ''masterCompute'', false, ''maxTrialNum'', 1, ''unitWaitTime'', 0, ', ...
                '''jobLogDir'', [folder, ''logs''], ''nice'', 1000)']);
            calls = slurm.calls('sbatch');
            testCase.assertNotEmpty(calls);
            testCase.verifyTrue(all(contains(calls, ' --nice=1000')), 'A job of the call went without its nice.');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function jobsAtAnotherPlaceKeepTheirNice(testCase)
            % with SlurmAlternatives a job goes to the first place with room
            % (here the alternative cpu_queue, as gpu8 is full): it keeps its nice
            % there, as a single job and as a framework job (a moved
            % pending job keeps it, as Slurm keeps its nice)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            clear slurm_alternatives generic_computing_frameworks_wrapper generic_single_job_submit_wrapper ...
                task_runner_run single_job_pools
            testCase.addTeardown(@clear, 'slurm_alternatives', 'generic_computing_frameworks_wrapper', ...
                'generic_single_job_submit_wrapper', 'task_runner_run', 'single_job_pools');
            slurm_alternatives('set', 'viewAge', 0);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({pktest.FakeSlurmFixture.nodeLine('n1', 'cpu_queue', 24, 0, 514048, 0, 0, 0, 'MIXED', 'c24'), ...
                pktest.FakeSlurmFixture.nodeLine('g1', 'gpu8', 128, 128, 1031600, 128 * 8000, 8, 0, 'MIXED', 'gpu8')});
            folder = testCase.createTempFolder();
            cpuParam = '-p cpu_queue --qos qos_normal -n1 --mem-per-cpu=21418M';
            conf = struct('SlurmParam', '-p gpu8 --qos qos_normal -n1 --mem-per-cpu=7800M', 'memPerCPU', 7.6, ...
                'maxCPUNum', 32, 'gpuHeadroom', 0, 'SlurmAlternatives', {{struct('SlurmParam', cpuParam, ...
                'memPerCPU', 20.9, 'maxCPUNum', 24)}});
            configFn = [folder, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(conf));
            fclose(fid);
            evalc(['generic_single_job_submit_wrapper(''disp(1)'', -1, 3, ''jobLogDir'', folder, ''memAllocate'', 70, ', ...
                '''ConfigFile'', configFn, ''taskRunner'', false, ''nice'', 1000)']);
            inputs = {[folder, 'in_1.txt']};
            fclose(fopen(inputs{1}, 'w'));
            outputs = {[folder, 'out_1.txt']}; %#ok<NASGU>
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 100));
            evalc(['generic_computing_frameworks_wrapper(inputs, outputs, {''disp(1)''}, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''memAllocate'', 70, ''jobLogDir'', folder, ', ...
                '''taskRunner'', false, ''queryInterval'', 0, ''ConfigFile'', configFn, ''nice'', 1000)']);
            jobs = slurm.calls('sbatch');
            jobs = jobs(~contains(jobs, '--test-only'));
            testCase.assertNumElements(jobs, 2);
            testCase.verifyTrue(all(contains(jobs, cpuParam)), 'The jobs did not go to the place with room.');
            testCase.verifyTrue(all(contains(jobs, ' --nice=1000')), strjoin(jobs, newline));
        end
    end

    methods (Access = private)
        function [] = stubScontrol(testCase, lines, status)
            % an scontrol first on the PATH that prints lines and exits with status
            folder = testCase.createTempFolder();
            fn = [folder, 'scontrol'];
            fid = fopen(fn, 'w');
            fprintf(fid, '#!/bin/bash\n');
            for i = 1 : numel(lines)
                fprintf(fid, 'echo "%s"\n', lines{i});
            end
            fprintf(fid, 'exit %d\n', status);
            fclose(fid);
            fileattrib(fn, '+x');
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('PATH', ...
                [folder, pathsep, getenv('PATH')]));
        end
    end
end
