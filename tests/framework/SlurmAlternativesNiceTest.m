classdef SlurmAlternativesNiceTest < pktest.TestCase
    % the nice of a call (its nice option, sbatch --nice) with
    % SlurmAlternatives (see slurm_alternatives): it is the call's, not a
    % place's. The jobs get it at every place, and pending jobs move between
    % places as without it (the main place's options had it and the
    % alternatives' didn't, and places that differ in options don't trade
    % jobs, so no job of a call with nice moved). The cluster is
    % pktest.FakeSlurmFixture's, with a place of CPU nodes (cpu: 24 CPUs,
    % 20.9 GB per CPU) and one of GPU nodes that takes CPU jobs (gpu: 128
    % CPUs, 7.6 GB per CPU); a job of 70 GB is 4 CPUs at cpu.

    properties (Constant)
        CpuParam = '-p cpu --qos normal -n1 --mem-per-cpu=21418M'
        GpuParam = '-p gpu --qos normal -n1 --mem-per-cpu=7800M'
    end

    methods (TestMethodSetup)
        function clearSession(testCase)
            % the view of the cluster, the recorded jobs, the timing
            % settings and the config cache are kept in persistent variables
            clear slurm_alternatives generic_computing_frameworks_wrapper generic_single_job_submit_wrapper ...
                task_runner_run single_job_pools
            testCase.addTeardown(@clear, 'slurm_alternatives', 'generic_computing_frameworks_wrapper', ...
                'generic_single_job_submit_wrapper', 'task_runner_run', 'single_job_pools');
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function pendingJobOfACallWithNiceMoves(testCase)
            % the job of batch 1 waits at cpu (full, gpu too); then gpu has
            % room (the second look at the nodes): the wait loop moves it
            % there, and it keeps its nice
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.cpuNode(24), testCase.gpuNode(128)}, {testCase.cpuNode(24), testCase.gpuNode(0)});
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 40));
            slurm.setJobStates('1000_1', [repmat({'PD:Priority'}, 1, 6), {'done'}]);
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0, 'viewAge', 0);
            testCase.runCall(testCase.config(testCase.GpuParam), 1000);
            jobs = testCase.submissions(slurm);
            testCase.assertNumElements(jobs, 1);
            testCase.verifySubstring(jobs{1}, [testCase.CpuParam, ' --nice=1000']);
            updates = slurm.calls('scontrol');
            updates = updates(startsWith(updates, 'update'));
            testCase.verifyEqual(updates, {['update JobId=1000_1 Partition=gpu QOS=normal NumCPUs=4 CPUsPerTask=4 ', ...
                'MinMemoryCPU=21418']}, 'The pending job of a call with nice did not move.');
        end

        function jobsGetTheCallsNiceAtEveryPlace(testCase)
            % a place's own --nice gave way to the call's only where it had
            % none: here the jobs go to gpu (cpu is full), whose SlurmParam
            % has --nice=5, and get the call's 1000 there, as framework jobs
            % and as single jobs
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.cpuNode(24), testCase.gpuNode(0)});
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 40));
            slurm_alternatives('set', 'viewAge', 0);
            configFn = testCase.config([testCase.GpuParam, ' --nice=5']);
            testCase.runCall(configFn, 1000);
            folder = fileparts(configFn); %#ok<NASGU>
            evalc(['generic_single_job_submit_wrapper(''disp(1)'', -1, 3, ''jobLogDir'', folder, ''memAllocate'', 70, ', ...
                '''ConfigFile'', configFn, ''taskRunner'', false, ''nice'', 1000)']);
            jobs = testCase.submissions(slurm);
            testCase.assertNumElements(jobs, 2);
            for i = 1 : 2
                testCase.verifySubstring(jobs{i}, '-p gpu ', jobs{i});
                nices = regexp(jobs{i}, '--nice=(\d+)', 'tokens');
                testCase.verifyEqual(cellfun(@(x) x{1}, nices, 'unif', 0), {'1000'}, jobs{i});
            end
        end
    end

    methods (Access = private)
        function [configFn] = config(testCase, altParam)
            % a cluster config of cpu with gpu (altParam) as its alternative
            folder = testCase.createTempFolder();
            conf = struct('SlurmParam', testCase.CpuParam, 'memPerCPU', 20.9, 'maxCPUNum', 24, ...
                'SlurmAlternatives', {{struct('SlurmParam', altParam, 'memPerCPU', 7.6, 'maxCPUNum', 32, ...
                'gpuHeadroom', 0, 'cpuJobs', true)}});
            configFn = [folder, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '%s', jsonencode(conf));
            fclose(fid);
        end

        function [] = runCall(~, configFn, nice)
            % a framework call of one task of 70 GB with nice (its job
            % never runs: the call gives up after its try)
            folder = [fileparts(configFn), '/'];
            inFns = {[folder, 'in_1.txt']};
            fclose(fopen(inFns{1}, 'w'));
            outFns = {[folder, 'out_1.txt']}; %#ok<NASGU>
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, {''disp(1)''}, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''memAllocate'', 70, ''jobLogDir'', folder, ', ...
                '''taskRunner'', false, ''queryInterval'', 0, ''ConfigFile'', configFn, ''nice'', nice)']);
        end

        function [jobs] = submissions(~, slurm)
            jobs = slurm.calls('sbatch');
            jobs = jobs(~contains(jobs, '--test-only'));
        end

        function [line] = cpuNode(~, cpuAlloc)
            line = pktest.FakeSlurmFixture.nodeLine('n1', 'cpu', 24, cpuAlloc, 514048, cpuAlloc * 21418, 0, 0, ...
                'MIXED', 'cpu24');
        end

        function [line] = gpuNode(~, cpuAlloc)
            line = pktest.FakeSlurmFixture.nodeLine('g1', 'gpu', 128, cpuAlloc, 1031600, cpuAlloc * 8000, 8, 0, ...
                'MIXED', 'gpu8');
        end
    end
end
