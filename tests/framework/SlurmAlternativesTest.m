classdef SlurmAlternativesTest < pktest.TestCase
    % jobs of a call on more than one place of a Slurm cluster (the
    % SlurmAlternatives of a cluster config, see slurm_alternatives): each
    % job goes to the first place with room for it now, and pending jobs
    % move to a place where they can start now. A job is shaped for the main
    % place; at a place with more memory per CPU it gets that place's own
    % shape (fewer CPUs), at one with less it keeps its CPUs and memory.
    % The cluster is pktest.FakeSlurmFixture's: nodes as scontrol show node
    % gives them, the user's pending jobs (squeue -O) and sbatch
    % --test-only answers. Places: abc (24-CPU nodes, 20.9 GB per CPU) and
    % gpu8 (128-CPU nodes with 8 GPUs, 7.6 GB per CPU). A job of 70 GB is 4
    % CPUs at abc and 10 at gpu8 (16 in the batch framework, rounded to a
    % factor of 32).

    properties (Constant)
        AbcParam = '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M'
        GpuParam = '-p gpu8 --qos abc_normal -n1 --mem-per-cpu=7800M'
    end

    methods (TestMethodSetup)
        function clearSession(testCase)
            % the view of the cluster, the recorded jobs and the timing
            % settings, the config cache and the single jobs' session state
            % are kept in persistent variables
            clear slurm_alternatives generic_computing_frameworks_wrapper generic_single_job_submit_wrapper ...
                task_runner_run single_job_pools
            testCase.addTeardown(@clear, 'slurm_alternatives', 'generic_computing_frameworks_wrapper', ...
                'generic_single_job_submit_wrapper', 'task_runner_run', 'single_job_pools');
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            % (each call reads the cluster again, unless a test keeps a view)
            slurm_alternatives('set', 'viewAge', 0);
        end
    end

    methods (Test, TestTags = {'Unit'})
        function placesTakeTheMainSettingsForKeysTheyLeaveOut(testCase)
            main = testCase.mainSettings();
            places = slurm_alternatives('places', main, []);
            testCase.verifyNumElements(places, 1);
            places = slurm_alternatives('places', main, {});
            testCase.verifyNumElements(places, 1);
            % a JSON list of objects with other keys is a cell array
            alts = jsondecode(sprintf(['[{"SlurmParam": "%s", "memPerCPU": 7.6, "maxCPUNum": 32}, ', ...
                '{"slurmparam": "-p x --qos q --mem-per-cpu=32G"}]'], testCase.GpuParam));
            places = slurm_alternatives('places', main, alts);
            testCase.assertNumElements(places, 3);
            testCase.verifyEqual({places.SlurmParam}, {testCase.AbcParam, testCase.GpuParam, '-p x --qos q --mem-per-cpu=32G'});
            % (memPerCPU: an entry's, else its --mem-per-cpu)
            testCase.verifyEqual([places.memPerCPU], [20.9, 7.6, 32]);
            testCase.verifyEqual([places.maxCPUNum], [24, 32, 24]);
            testCase.verifyEqual([places.minCPUNum], [1, 1, 1]);
            testCase.verifyEqual(places(2).opts.partition, 'gpu8');
            testCase.verifyEqual(places(2).opts.qos, 'abc_normal');
            testCase.verifyEqual(places(2).opts.memPerCPU, 7800);
            % a struct array (same keys) is the same; an unknown key is ignored with a warning
            places = testCase.verifyWarning(@() slurm_alternatives('places', main, ...
                struct('SlurmParam', testCase.GpuParam, 'maxCPU', 32)), 'slurm_alternatives:key');
            testCase.verifyEqual([places.maxCPUNum], [24, 24]);
        end

        function memoryPerCPUOfAPlaceWithoutOneComesFromItsNodes(testCase)
            % (the smallest RealMemory / CPUs of its nodes: 1031600 / 128 MB)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 0, 0), ...
                pktest.FakeSlurmFixture.nodeLine('g2', 'gpu8', 64, 0, 1031600, 0, 4, 0)});
            places = slurm_alternatives('places', testCase.mainSettings(), struct('SlurmParam', ...
                '-p gpu8 --qos abc_normal -n1 --mem=100G'));
            testCase.verifyEqual(places(2).memPerCPU, 1031600 / 128 / 1024, 'AbsTol', 1e-12);
        end

        function optionsOfASlurmParamAreRead(testCase)
            places = slurm_alternatives('places', testCase.mainSettings(), struct('SlurmParam', ...
                '-A lab --partition=gpu4 -q high -n1 --mem=78G --gres=gpu:a100:2 -C c24 -x g3.cluster -J pk -t 1:00:00 --exclusive', ...
                'memPerCPU', 30));
            o = places(2).opts;
            testCase.verifyEqual({o.account, o.partition, o.qos, o.gres, o.constraint, o.exclude, o.jobName, o.time, o.memToken}, ...
                {'lab', 'gpu4', 'high', 'gpu:a100:2', 'c24', 'g3.cluster', 'pk', '1:00:00', '--mem=78G'});
            testCase.verifyEqual([o.mem, o.gpus, o.ntasks], [78 * 1024, 2, 1]);
            testCase.verifyEqual(o.other, {'--exclusive'});
        end

        function jobGoesToTheFirstPlaceWithRoom(testCase)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            places = testCase.abcFirst(0);
            % abc has room: the main place
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 0, 0)});
            testCase.verifyEqual(testCase.choose(places), 1);
            % abc is full: gpu8
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            testCase.verifyEqual(testCase.choose(places), 2);
            % Slurm checked the place it chose
            testCase.verifyTrue(any(contains(slurm.calls('sbatch'), '--test-only') & contains(slurm.calls('sbatch'), '-p gpu8')));
        end

        function jobKeepsItsShapeAtAPlaceWithLessMemoryPerCPU(testCase)
            % (shaped for abc: 4 CPUs of 21418 MB; at gpu8 it keeps them,
            % not the 10 CPUs of 7.8 GB a gpu8 config gives it)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            [k, cpus, param] = slurm_alternatives('choose', testCase.abcFirst(0), [4, 10], {'', ''});
            testCase.verifyEqual({k, cpus, param}, {2, 4, '-p gpu8 --qos abc_normal -n1 --mem-per-cpu=21418M'});
            testOnly = slurm.calls('sbatch');
            testCase.verifySubstring(testOnly{end}, '--cpus-per-task=4 -p gpu8 --qos abc_normal -n1 --mem-per-cpu=21418M');
        end

        function jobIsReshapedAtAPlaceWithMoreMemoryPerCPU(testCase)
            % (shaped for gpu8: 10 CPUs of 7.8 GB; at abc 4 CPUs of 21418 MB)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 128, 0)});
            [k, cpus, param] = slurm_alternatives('choose', testCase.gpuFirst(), [10, 4], {'', ''});
            testCase.verifyEqual({k, cpus, param}, {2, 4, testCase.AbcParam});
        end

        function plannedNodesHaveNoRoom(testCase)
            % a node the backfill scheduler planned for a pending job keeps
            % its free CPUs for it
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 16, 'MIXED+PLANNED'), testCase.gpuNode('g1', 0, 0)});
            testCase.verifyEqual(testCase.choose(testCase.abcFirst(0)), 2);
        end

        function pendingJobsOfTheUserStartFirst(testCase)
            % abc has 8 free CPUs, which a pending job of the user takes
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 16), testCase.gpuNode('g1', 0, 0)});
            slurm.setJobStates('900_1', {'PD:Priority'});
            slurm.setJobAttributes('900_1', 'abc', 8, 8 * 21418);
            testCase.verifyEqual(testCase.choose(testCase.abcFirst(0)), 2);
            % one that fits no node waits, and a new job would wait behind
            % it (Slurm starts them by priority): no room at abc either
            slurm.setJobAttributes('900_1', 'abc', 32, 1000);
            testCase.verifyEqual(testCase.choose(testCase.abcFirst(0)), 2);
            % (only waiting jobs: a held one doesn't count)
            slurm.setJobStates('900_1', {'PD:JobHeldUser'});
            testCase.verifyEqual(testCase.choose(testCase.abcFirst(0)), 1);
        end

        function slurmsOwnCheckDecides(testCase)
            % the nodes have room, Slurm says the job would start later
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 0, 0)});
            slurm.setTestOnly('abc', 'later');
            testCase.verifyEqual(testCase.choose(testCase.abcFirst(0)), 2);
        end

        function placeSlurmRefusesIsNotUsed(testCase)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            slurm.setTestOnly('gpu8', 'error');
            places = testCase.abcFirst(0);
            k = testCase.verifyWarning(@() slurm_alternatives('choose', places, [4, 10], {'', ''}), ...
                'slurm_alternatives:testOnly');
            testCase.verifyEqual(k, 1);
            n = slurm.callCount('sbatch');
            testCase.verifyEqual(testCase.choose(places), 1);
            testCase.verifyEqual(slurm.callCount('sbatch'), n, 'Slurm was asked again about a place it refused.');
        end

        function mainPlaceIsNotLeftOutWhenSlurmRefusesItOnce(testCase)
            % (no room there then, but it is the user's place: asked again)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 0, 0)});
            slurm.setTestOnly('abc', 'error');
            places = testCase.abcFirst(0);
            testCase.verifyEqual(testCase.choose(places), 2);
            slurm.setTestOnly('abc', 'now');
            testCase.verifyEqual(testCase.choose(places), 1);
        end

        function burstOfJobsFillsTheRoomInOrder(testCase)
            % (one view of the cluster for all jobs of a burst) abc has room
            % for 2 jobs (4 CPUs of 21418 MB), gpu8 for 3 of them (its free
            % memory); the others wait at abc (not in the GPU nodes' queue)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 16), testCase.gpuNode('g1', 96, 0)});
            places = testCase.abcFirst(0);
            slurm_alternatives('set', 'viewAge', 60);
            ks = arrayfun(@(x) testCase.choose(places), 1 : 7);
            testCase.verifyEqual(ks, [1, 1, 2, 2, 2, 1, 1]);
            testCase.verifyEqual(slurm.callCount('scontrol'), 1, 'The nodes were read more than once in a view.');
            testOnly = contains(slurm.calls('sbatch'), '--test-only');
            testCase.verifyEqual(nnz(testOnly), 2, 'Slurm was asked more than once per place and shape in a view.');
        end

        function cpuJobsLeaveIdleGPUsTheirShare(testCase)
            % a 8-GPU node with 6 idle GPUs keeps 16 CPUs and 1/8 of its
            % memory per idle GPU (the node's share), here all its free CPUs
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 32, 2)});
            testCase.verifyEqual(testCase.choose(testCase.abcFirst([])), 1);
            % with the share of the run's GPU jobs (4 CPUs, 32 GB): gpu8
            testCase.verifyEqual(testCase.choose(testCase.abcFirst([4, 32])), 2);
            % 6 x (4 CPUs, 32 GB) kept: 578992 MB free for 6 jobs of 85672 MB
            slurm_alternatives('set', 'viewAge', 60);
            places = testCase.abcFirst([4, 32]);
            ks = arrayfun(@(x) testCase.choose(places), 1 : 8);
            testCase.verifyEqual(ks, [2, 2, 2, 2, 2, 2, 1, 1]);
            % (CPUs kept) 32 free CPUs and all memory, 4 idle GPUs: 16 CPUs for 4 jobs
            slurm.setNodes({testCase.abcNode('n1', 24), pktest.FakeSlurmFixture.nodeLine('g1', 'gpu8', 128, 96, ...
                1031600, 0, 8, 4)});
            slurm_alternatives('set', 'viewAge', 60);
            ks = arrayfun(@(x) testCase.choose(places), 1 : 5);
            testCase.verifyEqual(ks, [2, 2, 2, 2, 1]);
        end

        function featuresAndNodeListsOfAPlaceCount(testCase)
            % a place's -C and -x: here no node of gpu8 is one of its nodes
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0), ...
                testCase.gpuNode('g2', 0, 0)});
            for alt = {[testCase.GpuParam, ' -C gpu1'], [testCase.GpuParam, ' -x g[1-2]'], ...
                    [testCase.GpuParam, ' -x g1,g2.cluster']}
                places = slurm_alternatives('places', testCase.mainSettings(), struct('SlurmParam', alt{1}, ...
                    'memPerCPU', 7.6, 'gpuHeadroom', 0));
                testCase.verifyEqual(testCase.choose(places), 1, alt{1});
            end
            places = slurm_alternatives('places', testCase.mainSettings(), struct('SlurmParam', ...
                [testCase.GpuParam, ' -C gpu8 -x g1'], 'memPerCPU', 7.6, 'gpuHeadroom', 0));
            testCase.verifyEqual(testCase.choose(places), 2);
        end

        function gpuJobsNeedAFreeGPU(testCase)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.gpuNode('g1', 8, 8), testCase.gpuNode('g2', 8, 7, 'MIXED', 'gpu9')});
            main = testCase.mainSettings();
            main.SlurmParam = '-p gpu8 --qos abc_normal -n1 --mem-per-cpu=8000M --gres=gpu:1';
            main.memPerCPU = 7.8;
            places = slurm_alternatives('places', main, struct('SlurmParam', ...
                '-p gpu9 --qos abc_normal -n1 --mem-per-cpu=8000M --gres=gpu:1')); %#ok<NASGU>
            [~, k] = evalc('slurm_alternatives(''choose'', places, [4, 4], {'''', ''''})');
            testCase.verifyEqual(k, 2);
        end

        function limitThatPendingJobsWaitForClosesAPlace(testCase)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 0, 0)});
            slurm.setJobStates('900_1', {'PD:QOSMaxCpuPerUserLimit'});
            slurm.setJobAttributes('900_1', 'abc', 4, 4 * 21418, 0, 'abc_normal');
            testCase.verifyEqual(testCase.choose(testCase.abcFirst(0)), 2);
        end

        function cpuJobsWaitOutsideGPUNodes(testCase)
            % main place gpu8 (GPU nodes), nothing has room:
            % the job waits at abc (which has GPU nodes too, not only)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 128, 0), ...
                pktest.FakeSlurmFixture.nodeLine('t1', 'abc', 20, 20, 515638, 0, 4, 4)});
            [~, k] = evalc('slurm_alternatives(''choose'', testCase.gpuFirst(), [10, 4], {'''', ''''})');
            testCase.verifyEqual(k, 2);
        end

        function pendingJobMovesToAPlaceWithRoom(testCase)
            % (at gpu8 it keeps its abc shape)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            places = testCase.abcFirst(0);
            testCase.submitPending(slurm, places, 1000, 1, 1);
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
            evalc('slurm_alternatives(''move'', 1000, 1)');
            testCase.verifyEqual(testCase.updates(slurm), {['update JobId=1000_1 Partition=gpu8 QOS=abc_normal ', ...
                'NumCPUs=4 CPUsPerTask=4 MinMemoryCPU=21418']});
            attrs = slurm.jobAttributes('1000_1');
            testCase.verifyEqual({attrs.partition, attrs.cpus, attrs.mem}, {'gpu8', 4, 4 * 21418});
            % not again within 2 minutes of the move
            slurm_alternatives('set', 'minWait', 120);
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 128, 0)});
            n = slurm.callCount('scontrol');
            evalc('slurm_alternatives(''move'', 1000, 1)');
            testCase.verifyEqual(slurm.callCount('scontrol'), n, 'A job moved again within 2 minutes.');
        end

        function onlyAsManyJobsMoveAsFit(testCase)
            % gpu8 has room for 2 jobs of 4 CPUs of 21418 MB
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 106, 0)});
            places = testCase.abcFirst(0);
            for a = 1 : 3
                testCase.submitPending(slurm, places, 1000 + a, a, 1);
            end
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
            evalc('slurm_alternatives(''move'', [1001; 1002; 1003], [1; 2; 3])');
            testCase.verifyEqual(regexprep(testCase.updates(slurm), '^update JobId=(\S+) .*$', '$1'), {'1001_1'; '1002_2'});
        end

        function jobsThatDoNotWaitForResourcesStay(testCase)
            % held, waiting for a start time or dependency, not seen by the
            % scheduler yet, running, or with room where they are (abc has
            % room for 2: the one not seen yet and the last one)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 16), testCase.gpuNode('g1', 0, 0)});
            places = testCase.abcFirst(0);
            states = {'PD:JobHeldUser', 'PD:BeginTime', 'PD:Dependency', 'PD:None', 'R', 'PD:Priority'};
            for a = 1 : numel(states)
                testCase.submitPending(slurm, places, 1000 + a, a, 1);
                slurm.setJobStates(sprintf('%d_%d', 1000 + a, a), states(a));
            end
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
            evalc('slurm_alternatives(''move'', 1000 + (1 : 6)'', (1 : 6)'')');
            testCase.verifyEmpty(testCase.updates(slurm));
        end

        function jobWithRoomOnlyOutsideItsPlaceMoves(testCase)
            % (its place is abc without n1; the room at abc is on n1)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.abcNode('n2', 24), testCase.gpuNode('g1', 0, 0)});
            main = testCase.mainSettings();
            main.SlurmParam = [testCase.AbcParam, ' -x n1'];
            places = slurm_alternatives('places', main, struct('SlurmParam', [testCase.GpuParam, ' -x n1'], ...
                'memPerCPU', 7.6, 'maxCPUNum', 32, 'gpuHeadroom', 0));
            testCase.submitPending(slurm, places, 1000, 1, 1);
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
            evalc('slurm_alternatives(''move'', 1000, 1)');
            testCase.verifyNumElements(testCase.updates(slurm), 1);
        end

        function pendingJobsAreFoundByTheirIds(testCase)
            % (the first pending job in squeue's list is held; the second
            % one waits for resources and moves)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            places = testCase.abcFirst(0);
            testCase.submitPending(slurm, places, 1001, 1, 1);
            slurm.setJobStates('1001_1', {'PD:JobHeldUser'});
            testCase.submitPending(slurm, places, 1002, 2, 1);
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
            evalc('slurm_alternatives(''move'', [1001; 1002], [1; 2])');
            testCase.verifyEqual(regexprep(testCase.updates(slurm), '^update JobId=(\S+) .*$', '$1'), {'1002_2'});
        end

        function noSlurmCallsBeforeAJobHasWaited(testCase)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            places = testCase.abcFirst(0);
            testCase.submitPending(slurm, places, 1000, 1, 1);
            slurm.clearCalls();
            evalc('slurm_alternatives(''move'', 1000, 1)');
            testCase.verifyEmpty(slurm.history(), 'Slurm was asked about a job that just came.');
            % jobs not recorded (of other calls or users) are never looked at
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
            evalc('slurm_alternatives(''move'', 2000, 1)');
            testCase.verifyEmpty(slurm.history());
        end

        function cpuJobLeavesTheQueueOfGPUNodes(testCase)
            % a CPU-only job waiting at gpu8 (main place) without room there
            % or anywhere: it waits at abc instead, shaped for abc
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 128, 0)});
            places = testCase.gpuFirst();
            slurm_alternatives('submitted', 1000, 1, places, [10, 4], 1, {'', ''});
            slurm.setJobStates('1000_1', {'PD:Resources'});
            slurm.setJobAttributes('1000_1', 'gpu8', 10, 78000, 0, 'abc_normal');
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
            evalc('slurm_alternatives(''move'', 1000, 1)');
            testCase.verifyEqual(testCase.updates(slurm), {['update JobId=1000_1 Partition=abc ', ...
                'QOS=abc_normal NumCPUs=4 CPUsPerTask=4 MinMemoryCPU=21418']});
        end

        function placesWithOtherOptionsDoNotTrade(testCase)
            % options a pending job can't change (here --exclusive), and a
            % place without a QOS for jobs of another partition (Slurm keeps
            % the job's QOS, also one the new partition doesn't allow)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            for alt = {[testCase.GpuParam, ' --exclusive'], '-p gpu8 -n1 --mem-per-cpu=7800M'}
                places = slurm_alternatives('places', testCase.mainSettings(), struct('SlurmParam', alt{1}, ...
                    'memPerCPU', 7.6, 'maxCPUNum', 32, 'gpuHeadroom', 0));
                testCase.submitPending(slurm, places, 1000, 1, 1);
                slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
                evalc('slurm_alternatives(''move'', 1000, 1)');
                testCase.verifyEmpty(testCase.updates(slurm), alt{1});
            end
        end

        function jobMovedBySomethingElseIsLeftToIt(testCase)
            % (recorded at abc, found at gpu8: e.g. a user's script moved
            % it; gpu9 has room)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 128, 0), ...
                testCase.gpuNode('g2', 0, 0, 'MIXED', 'gpu9')});
            places = slurm_alternatives('places', testCase.mainSettings(), struct('SlurmParam', ...
                {testCase.GpuParam, strrep(testCase.GpuParam, 'gpu8', 'gpu9')}, 'memPerCPU', 7.6, 'gpuHeadroom', 0));
            testCase.submitPending(slurm, places, 1000, 1, 1);
            slurm.setJobAttributes('1000_1', 'gpu8', 6, 6 * 21418, 0, 'abc_normal');
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
            evalc('slurm_alternatives(''move'', 1000, 1)');
            % (then abc has room: still not moved)
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 128, 0), ...
                testCase.gpuNode('g2', 0, 0, 'MIXED', 'gpu9')});
            slurm_alternatives('set', 'minWait', 0);
            evalc('slurm_alternatives(''move'', 1000, 1)');
            testCase.verifyEmpty(testCase.updates(slurm));
        end

        function failedMoveLeavesTheJob(testCase)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            places = testCase.abcFirst(0);
            testCase.submitPending(slurm, places, 1000, 1, 1);
            slurm.failUpdates();
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
            out = evalc('slurm_alternatives(''move'', 1000, 1)');
            testCase.verifySubstring(out, 'was not moved');
            testCase.verifyEqual(slurm.jobAttributes('1000_1').partition, 'abc');
        end

        function frameworkSubmitsToThePlaceWithRoom(testCase)
            % a main place of GPU nodes (for CPU jobs), whose idle
            % GPUs keep their share; abc has room: the jobs go to abc with
            % the CPUs of an abc config (4, not gpu8's 16), and the time
            % limit of its SlurmParam (not the framework's)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 0, 0)});
            [inFns, outFns, funcStrs, dataPath] = testCase.makeTasks(2); %#ok<ASGLU>
            conf = testCase.gpuFirstConfig();
            conf.SlurmAlternatives{1}.SlurmParam = [testCase.AbcParam, ' -t 2:00:00'];
            configFn = testCase.writeConfig(dataPath, conf); %#ok<NASGU>
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''memAllocate'', 70, ''jobLogDir'', dataPath, ', ...
                '''taskRunner'', false, ''ConfigFile'', configFn)']);
            jobs = testCase.submissions(slurm);
            testCase.assertNumElements(jobs, 2);
            testCase.verifyTrue(all(contains(jobs, ['--cpus-per-task=4 ', testCase.AbcParam, ' -t 2:00:00 '])), strjoin(jobs, newline));
            testCase.verifyFalse(any(contains(jobs, ' -t 24:00:00 ')), strjoin(jobs, newline));
        end

        function compiledFrameworkSubmitsToThePlaceWithRoom(testCase)
            % (compiled mode: a fake MATLAB Runtime and compiled master)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 0, 0)});
            dataPath = testCase.createTempFolder();
            mkdir([dataPath, 'mcr/bin']);
            masterFn = [dataPath, 'run_mccMaster.sh'];
            fid = fopen(masterFn, 'w');
            fprintf(fid, '#!/bin/bash\necho "$MCR_CACHE_ROOT" >> "%scaches.log"\n', dataPath);
            fclose(fid);
            fileattrib(masterFn, '+x');
            conf = testCase.gpuFirstConfig();
            conf.MCRParam = [dataPath, 'mcr/'];
            conf.MCCMasterStr = masterFn;
            conf.MCRCacheRoot = [dataPath, 'cache'];
            conf.BashLaunchStr = 'echo';
            configFn = testCase.writeConfig(dataPath, conf); %#ok<NASGU>
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : 2, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout_%d.txt', dataPath, x), 1 : 2, 'unif', 0);
            cellfun(@(x) fclose(fopen(x, 'w')), inFns);
            funcStrs = cellfun(@(x, y) sprintf('XR_crop_frame(''%s'',''%s'',[1,1,1,2,2,2],''pad'',false)', x, y), ...
                inFns, outFns, 'unif', 0); %#ok<NASGU>
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''mccMode'', true, ''maxTrialNum'', 1, ''memAllocate'', 70, ', ...
                '''jobLogDir'', dataPath, ''ConfigFile'', configFn)']);
            jobs = testCase.submissions(slurm);
            testCase.assertNumElements(jobs, 2);
            testCase.verifyTrue(all(contains(jobs, ['--cpus-per-task=4 ', testCase.AbcParam, ' '])), strjoin(jobs, newline));
            testCase.verifyTrue(all(contains(jobs, ' -t 24:00:00 ')), strjoin(jobs, newline));
        end

        function frameworkMovesPendingJobs(testCase)
            % the job of batch 1 waits at abc (full, gpu8 too); then gpu8 has
            % room (the second look at the nodes): the wait loop moves it
            % there, with its abc shape, and not again once it fits there
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 128, 0)}, ...
                {testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 40));
            slurm.setJobStates('1000_1', [repmat({'PD:Priority'}, 1, 6), {'done'}]);
            [inFns, outFns, funcStrs, dataPath] = testCase.makeTasks(1); %#ok<ASGLU>
            configFn = testCase.writeConfig(dataPath, testCase.abcFirstConfig(0)); %#ok<NASGU>
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0, 'viewAge', 0);
            evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                '''masterCompute'', false, ''maxTrialNum'', 1, ''memAllocate'', 70, ''jobLogDir'', dataPath, ', ...
                '''taskRunner'', false, ''queryInterval'', 0, ''ConfigFile'', configFn)']);
            jobs = testCase.submissions(slurm);
            testCase.assertNumElements(jobs, 1);
            testCase.verifySubstring(jobs{1}, ['--cpus-per-task=4 ', testCase.AbcParam, ' ']);
            testCase.verifyEqual(testCase.updates(slurm), {['update JobId=1000_1 Partition=gpu8 ', ...
                'QOS=abc_normal NumCPUs=4 CPUsPerTask=4 MinMemoryCPU=21418']});
        end

        function singleJobGoesToThePlaceWithRoomAndMoves(testCase)
            % the pipeline's and decon's single jobs (main place gpu8 without
            % GPU headroom): submitted to abc (gpu8 full) with an abc shape,
            % then moved back to gpu8 once it has room and abc not (keeping
            % the abc shape), by the wrapper's and check_single_job_status's
            % looks at pending jobs
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 128, 0)});
            dataPath = testCase.createTempFolder();
            conf = testCase.gpuFirstConfig();
            conf.gpuHeadroom = 0;
            configFn = testCase.writeConfig(dataPath, conf);
            args = {'jobLogDir', dataPath, 'memAllocate', 70, 'ConfigFile', configFn, 'taskRunner', false}; %#ok<NASGU>
            [~, jobId] = evalc('generic_single_job_submit_wrapper(''disp(1)'', -1, 3, args{:})');
            jobs = testCase.submissions(slurm);
            testCase.assertNumElements(jobs, 1);
            testCase.verifySubstring(jobs{1}, ['--cpus-per-task=4 ', testCase.AbcParam, ' ']);
            slurm.setJobStates(sprintf('%d_3', jobId), {'PD:Priority'});
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            slurm_alternatives('set', 'minWait', 0, 'lookInterval', 0);
            [~, ~, status] = evalc('generic_single_job_submit_wrapper(''disp(1)'', jobId, 3, args{:})');
            testCase.verifyEqual(status, 0);
            testCase.verifyEqual(testCase.updates(slurm), {sprintf(['update JobId=%d_3 Partition=gpu8 ', ...
                'QOS=abc_normal NumCPUs=4 CPUsPerTask=4 MinMemoryCPU=21418'], jobId)});
            % and check_single_job_status looks at it as well
            slurm.setNodes({testCase.abcNode('n1', 0), testCase.gpuNode('g1', 128, 0)});
            slurm_alternatives('set', 'minWait', 0);
            evalc('check_single_job_status(jobId, 3)');
            testCase.verifyNumElements(testCase.updates(slurm), 2);
        end

        function configWithoutAlternativesIsAsBefore(testCase)
            % an empty list: the same jobs as without the key, and no Slurm
            % command but sinfo, sbatch and squeue -j
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setNodes({testCase.abcNode('n1', 24), testCase.gpuNode('g1', 0, 0)});
            [inFns, outFns, funcStrs, dataPath] = testCase.makeTasks(2); %#ok<ASGLU>
            conf = testCase.abcFirstConfig(0);
            conf.SlurmAlternatives = {};
            configs = {testCase.writeConfig(dataPath, rmfield(conf, 'SlurmAlternatives'), 'a.json'), ...
                testCase.writeConfig(dataPath, conf, 'b.json')};
            jobs = cell(1, 2);
            for i = 1 : 2
                slurm.clearCalls();
                configFn = configs{i}; %#ok<NASGU>
                evalc(['generic_computing_frameworks_wrapper(inFns, outFns, funcStrs, ''parseCluster'', true, ', ...
                    '''masterCompute'', false, ''maxTrialNum'', 1, ''memAllocate'', 70, ''jobLogDir'', dataPath, ', ...
                    '''taskRunner'', false, ''uuid'', ''u1'', ''ConfigFile'', configFn)']);
                jobs{i} = slurm.calls('sbatch');
                cmds = slurm.history();
                testCase.verifyTrue(all(startsWith(cmds, {'sinfo', 'sbatch', 'squeue -j', 'scancel'})), strjoin(cmds, newline));
                testCase.verifyFalse(any(contains(cmds, '--test-only')));
            end
            testCase.verifyEqual(jobs{2}, jobs{1});
        end
    end

    methods (Access = private)
        function [main] = mainSettings(testCase)
            main = struct('SlurmParam', testCase.AbcParam, 'SlurmConstraint', '', 'memPerCPU', 20.9, ...
                'maxCPUNum', 24, 'minCPUNum', 1, 'gpuHeadroom', []);
        end

        function [places] = abcFirst(testCase, headroom)
            % abc first, then gpu8 with the given GPU headroom
            places = slurm_alternatives('places', testCase.mainSettings(), struct('SlurmParam', ...
                testCase.GpuParam, 'memPerCPU', 7.6, 'maxCPUNum', 32, 'gpuHeadroom', headroom));
        end

        function [places] = gpuFirst(testCase)
            main = struct('SlurmParam', testCase.GpuParam, 'SlurmConstraint', '', 'memPerCPU', 7.6, ...
                'maxCPUNum', 32, 'minCPUNum', 1, 'gpuHeadroom', []);
            places = slurm_alternatives('places', main, struct('SlurmParam', testCase.AbcParam, ...
                'memPerCPU', 20.9, 'maxCPUNum', 24));
        end

        function [k] = choose(~, places) %#ok<INUSD>
            % a job of 70 GB: 4 CPUs at abc, 10 at gpu8 (single-job shape)
            [~, k] = evalc('slurm_alternatives(''choose'', places, [4, 10], {'''', ''''})');
        end

        function [] = submitPending(~, slurm, places, jobId, arrayId, k)
            % a job submitted to place k (shaped for abc: 4 CPUs of 21418
            % MB) that waits for resources
            slurm_alternatives('submitted', jobId, arrayId, places, [4, 10], k, {'', ''});
            id = sprintf('%d_%d', jobId, arrayId);
            slurm.setJobStates(id, {'PD:Priority'});
            slurm.setJobAttributes(id, places(k).opts.partition, 4, 4 * 21418, 0, places(k).opts.qos);
        end

        function [updates] = updates(~, slurm)
            updates = slurm.calls('scontrol');
            updates = updates(startsWith(updates, 'update'));
        end

        function [jobs] = submissions(~, slurm)
            jobs = slurm.calls('sbatch');
            jobs = jobs(~contains(jobs, '--test-only'));
        end

        function [line] = abcNode(~, name, cpuAlloc, state)
            if nargin < 4
                state = 'MIXED';
            end
            line = pktest.FakeSlurmFixture.nodeLine(name, 'abc', 24, cpuAlloc, 514048, cpuAlloc * 21418, ...
                0, 0, state, 'c24');
        end

        function [line] = gpuNode(~, name, cpuAlloc, gpuAlloc, state, partition)
            if nargin < 5
                state = 'MIXED';
            end
            if nargin < 6
                partition = 'gpu8';
            end
            line = pktest.FakeSlurmFixture.nodeLine(name, partition, 128, cpuAlloc, 1031600, ...
                cpuAlloc * 8000, 8, gpuAlloc, state, 'gpu8');
        end

        function [conf] = abcFirstConfig(testCase, headroom)
            conf = struct('SlurmParam', testCase.AbcParam, 'memPerCPU', 20.9, 'maxCPUNum', 24, ...
                'SlurmAlternatives', {{struct('SlurmParam', testCase.GpuParam, 'memPerCPU', 7.6, ...
                'maxCPUNum', 32, 'gpuHeadroom', headroom)}});
        end

        function [conf] = gpuFirstConfig(testCase)
            % (a config for CPU jobs on GPU nodes, with abc as an alternative)
            conf = struct('SlurmParam', testCase.GpuParam, 'memPerCPU', 7.6, 'maxCPUNum', 32, ...
                'SlurmAlternatives', {{struct('SlurmParam', testCase.AbcParam, 'memPerCPU', 20.9, 'maxCPUNum', 24)}});
        end

        function [fn] = writeConfig(~, dataPath, conf, name)
            if nargin < 4
                name = 'config.json';
            end
            fn = [dataPath, name];
            fid = fopen(fn, 'w');
            fprintf(fid, '%s', jsonencode(conf));
            fclose(fid);
        end

        function [inFns, outFns, funcStrs, dataPath] = makeTasks(testCase, n)
            dataPath = testCase.createTempFolder();
            mkdir([dataPath, 'out']);
            inFns = arrayfun(@(x) sprintf('%sin_%d.txt', dataPath, x), 1 : n, 'unif', 0);
            outFns = arrayfun(@(x) sprintf('%sout/out_%d.txt', dataPath, x), 1 : n, 'unif', 0);
            for f = 1 : n
                fclose(fopen(inFns{f}, 'w'));
            end
            funcStrs = cellfun(@(x, y) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''noop'')', ...
                x, y, [dataPath, 'calls.log']), inFns, outFns, 'unif', 0);
        end
    end
end
