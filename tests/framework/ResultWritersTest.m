classdef ResultWritersTest < pktest.TestCase
    % result_writers: a task that makes a result in parts with jobs of its
    % own records its process, and another task of the same result (another
    % uuid) waits while the first may still work on it. Owners outside Slurm
    % are live while their process runs on this node (Linux); owners in
    % Slurm while their job is in the queue (a fake Slurm); none once their
    % record says they have ended.

    methods (TestMethodSetup)
        function forgetSqueueAnswers(~)
            % (squeue's answers are kept per job for a minute in the
            % process, and the tests' fake job ids repeat)
            clear('result_writers');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function liveWritersAreThoseOfOtherUuids(testCase)
            testCase.assumeTrue(isunix && ~ismac, 'Writers outside Slurm are checked on Linux only.');
            [resultFn, flagRoot] = testCase.makeResultFolder();
            pid = testCase.startProcess();
            node = testCase.nodeName();
            testCase.writeRecord(flagRoot, resultFn, 'aaa', sprintf('0_:%s:%d', node, pid));
            testCase.writeRecord(flagRoot, resultFn, 'bbb', sprintf('0_:%s:%d', node, 999999999));
            testCase.writeRecord(flagRoot, resultFn, 'ccc', sprintf('0_:othernode:%d', pid));
            testCase.writeRecord(flagRoot, [resultFn(1 : end - 5), '_x.zarr'], 'ddd', sprintf('0_:%s:%d', node, pid));
            testCase.writeRecord(flagRoot, resultFn, 'eee', sprintf('0_:%s:%d:1', node, pid));
            fff = sprintf('0_:%s:%d:%s', node, pid, testCase.processStart(pid));
            testCase.writeRecord(flagRoot, resultFn, 'fff', fff);
            testCase.writeRecord(flagRoot, resultFn, 'ggg', fff, 0, 'ended=1');
            % live: aaa, and fff (its process, started at the time
            % recorded); bbb's process ended, ccc runs elsewhere (not
            % checkable), ddd is of another result, eee's process id is used
            % by another process now (started at another time), ggg's
            % writer has ended (its process stays, e.g. an interactive
            % session)
            testCase.verifyEqual(sort(result_writers('live', resultFn, 'zzz')), sort({sprintf('0_:%s:%d', node, pid), fff}));
            % (and not one's own)
            testCase.verifyEqual(result_writers('live', resultFn, 'aaa'), {fff});
        end

        function ownerJobsAreLiveWhileInTheQueue(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [resultFn, flagRoot] = testCase.makeResultFolder();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('777_3', {'PD', 'done'});
            slurm.setJobStates('778', {'R'});
            testCase.writeRecord(flagRoot, resultFn, 'aaa', '777_3:node1:11');
            testCase.writeRecord(flagRoot, resultFn, 'bbb', '778_:node2:12');
            testCase.verifyEqual(sort(result_writers('live', resultFn, 'zzz')), {'777_3:node1:11', '778_:node2:12'});
            testCase.verifyEqual(result_writers('live', resultFn, 'zzz'), {'778_:node2:12'}, ...
                'A job that ended still counted.');
            calls = slurm.calls('squeue');
            testCase.verifyTrue(any(contains(calls, '-j 778')) && ~any(contains(calls, '-j 778_')), ...
                'A job that is not an array job was checked as one.');
        end

        function writersOfThisProcesssJobAreItsProcesses(testCase)
            % a record of another uuid owned by this process's own Slurm job
            % (a call this process made before and stopped, or a job started
            % again with the same id): squeue lists the job, so this process
            % waited for itself. Now the process counts, as outside Slurm:
            % this one is not waited for, another of the job on this node is
            % while it runs; squeue is not asked about the job.
            testCase.assumeTrue(isunix && ~ismac, 'Writers outside Slurm are checked on Linux only.');
            [resultFn, flagRoot] = testCase.makeResultFolder();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('4242', {'R'});
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_JOB_ID', '4242'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_ARRAY_JOB_ID', ''));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_ARRAY_TASK_ID', ''));
            node = testCase.nodeName();
            me = feature('getpid');
            testCase.writeRecord(flagRoot, resultFn, 'aaa', sprintf('4242_:%s:%d:%s', node, me, testCase.processStart(me)), -60);
            mkdir([flagRoot, 'result_bbb/']);
            result_writers('register', resultFn, 'bbb', [flagRoot, 'result_bbb/']);
            testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', [], 0));
            testCase.verifyFalse(result_writers('wait', resultFn, 'bbb'), 'This process waited for itself.');
            pid = testCase.startProcess();
            other = sprintf('4242_:%s:%d:%s', node, pid, testCase.processStart(pid));
            testCase.writeRecord(flagRoot, resultFn, 'ccc', other, -60);
            testCase.verifyEqual(result_writers('live', resultFn, 'bbb'), {other});
            testCase.verifyFalse(any(contains(slurm.calls('squeue'), '4242')), 'squeue was asked about this process''s job.');
        end

        function waitEndsWhenTheResultComes(testCase)
            % a writer recorded before waits; the wait looks again each time
            % wait_for_files returns (a spy that makes the result at its
            % second call)
            testCase.assumeTrue(isunix && ~ismac, 'Writers outside Slurm are checked on Linux only.');
            [resultFn, flagRoot] = testCase.makeResultFolder();
            pid = testCase.startProcess();
            testCase.writeRecord(flagRoot, resultFn, 'aaa', sprintf('0_:%s:%d', testCase.nodeName(), pid), -60);
            flagPath = [flagRoot, 'result_bbb/'];
            mkdir(flagPath);
            result_writers('register', resultFn, 'bbb', flagPath);
            waits = testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', @() pktest.atCall( ...
                [flagRoot, 'calls.txt'], 2, @() mkdir(resultFn)), 10));
            done = [];
            out = evalc('done = result_writers(''wait'', resultFn, ''bbb'');');
            testCase.verifyTrue(done);
            testCase.verifyNumElements(waits.getCalls(), 2);
            testCase.verifySubstring(out, 'wait for it');
        end

        function writerRecordedLaterDoesNotWaitForAnEarlierOne(testCase)
            % two writers: the later one waits for the earlier one, not the
            % other way round (else both would wait), and of two recorded at
            % the same time, the one with the larger uuid
            testCase.assumeTrue(isunix && ~ismac, 'Writers outside Slurm are checked on Linux only.');
            [resultFn, flagRoot] = testCase.makeResultFolder();
            pid = testCase.startProcess();
            node = testCase.nodeName();
            mkdir([flagRoot, 'result_aaa/']);
            result_writers('register', resultFn, 'aaa', [flagRoot, 'result_aaa/']);
            testCase.writeRecord(flagRoot, resultFn, 'bbb', sprintf('0_:%s:%d', node, pid), 60);
            % aaa came first: it doesn't wait for bbb
            testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', [], 0));
            testCase.verifyFalse(result_writers('wait', resultFn, 'aaa'));
            % the same time: aaa < bbb, so aaa goes first
            aaaFn = [flagRoot, 'result_aaa/writer.txt'];
            fileObj = java.io.File([flagRoot, 'result_bbb/writer.txt']);
            fileObj.setLastModified(java.io.File(aaaFn).lastModified());
            testCase.verifyFalse(result_writers('wait', resultFn, 'aaa'));
            testCase.verifyEqual(numel(result_writers('live', resultFn, 'aaa')), 1);
        end

        function releaseMarksOnlyThisProcesssRecordEnded(testCase)
            % a writer that ends but by publishing (an error, Ctrl+C) marks
            % its record ended, once; not a record of another process, and
            % no record again once the flags are gone (the result published)
            [resultFn, flagRoot] = testCase.makeResultFolder();
            flagPath = [flagRoot, 'result_aaa/'];
            mkdir(flagPath);
            result_writers('register', resultFn, 'aaa', flagPath);
            % (of another uuid: not this one's)
            result_writers('release', resultFn, 'zzz', flagPath);
            testCase.verifyEmpty(strfind(fileread([flagPath, 'writer.txt']), 'ended='), 'A record of another uuid was marked.');
            result_writers('release', resultFn, 'aaa', flagPath);
            result_writers('release', resultFn, 'aaa', flagPath);
            txt = fileread([flagPath, 'writer.txt']);
            testCase.verifyTrue(endsWith(txt, sprintf('\nended=1\n')), txt);
            testCase.verifyNumElements(strfind(txt, 'ended='), 1);
            testCase.writeRecord(flagRoot, resultFn, 'bbb', '777_3:node1:11');
            before = fileread([flagRoot, 'result_bbb/writer.txt']);
            result_writers('release', resultFn, 'bbb', [flagRoot, 'result_bbb/']);
            testCase.verifyEqual(fileread([flagRoot, 'result_bbb/writer.txt']), before, 'Another process''s record was changed.');
            rmdir(flagPath, 's');
            result_writers('release', resultFn, 'aaa', flagPath);
            testCase.verifyFalse(isfolder(flagPath), 'The flags of a published result were made again.');
        end

        function blocksOfAWriterWhoseJobEndedAreOrphans(testCase)
            % a block's writer whose Slurm job squeue no longer lists (e.g.,
            % cancelled), or whose record says it has ended, is gone:
            % nothing will use the block. A writer still in the queue, one
            % outside Slurm, one that is this process's own job, no record,
            % and a failing squeue: not orphans.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [resultFn, flagRoot] = testCase.makeResultFolder();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('778_2', {'R'});
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_JOB_ID', '4242'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_ARRAY_JOB_ID', ''));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_ARRAY_TASK_ID', ''));
            cases = {'777_3:node1:11', '', true; '778_2:node1:12', '', false; '0_:node1:13', '', false; ...
                '4242_:node1:14', '', false; '778_2:node1:12', 'ended=1', true; '0_:node1:13', 'ended=1', true};
            for i = 1 : size(cases, 1)
                uuid = sprintf('u%d', i);
                testCase.writeRecord(flagRoot, resultFn, uuid, cases{i, 1}, 0, cases{i, 2});
                flagFn = sprintf('%sresult_%s/blocks_1_5.mat', flagRoot, uuid);
                orphaned = [];
                evalc('orphaned = result_writers(''orphaned'', flagFn);');
                testCase.verifyEqual(orphaned, cases{i, 3}, sprintf('%s %s', cases{i, 1 : 2}));
            end
            testCase.verifyFalse(result_writers('orphaned', [flagRoot, 'result_none/blocks_1_5.mat']));
        end

        function failingSqueueDoesNotMakeOrphans(testCase)
            % (the test Slurm stubs fail, as squeue does when the controller is busy)
            testCase.assumeFalse(ispc, 'Slurm commands are stubbed on Linux and macOS only.');
            [resultFn, flagRoot] = testCase.makeResultFolder();
            testCase.writeRecord(flagRoot, resultFn, 'u1', '777_3:node1:11');
            testCase.verifyFalse(result_writers('orphaned', [flagRoot, 'result_u1/blocks_1_5.mat']));
        end

        function squeueAnswersAreKeptForAMinute(testCase)
            % a job that runs many blocks of one writer (task runners) asked
            % squeue at each block; squeue's answer about a job is now kept
            % in the process for a minute: 2 checks make 1 call. An older
            % answer (maxAge) is asked again, and a record that says its
            % writer has ended counts at once (the record is read each time).
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [resultFn, flagRoot] = testCase.makeResultFolder();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('778_2', {'R'});
            testCase.writeRecord(flagRoot, resultFn, 'u1', '778_2:node1:12');
            flagFn = [flagRoot, 'result_u1/blocks_1_5.mat'];
            testCase.verifyFalse(result_writers('orphaned', flagFn));
            testCase.verifyFalse(result_writers('orphaned', [flagRoot, 'result_u1/blocks_6_10.mat']));
            testCase.verifyEqual(slurm.callCount('squeue'), 1, 'squeue was asked again within the minute.');
            testCase.verifyFalse(result_writers('orphaned', flagFn, 0));
            testCase.verifyEqual(slurm.callCount('squeue'), 2, 'An older answer was used.');
            fid = fopen([flagRoot, 'result_u1/writer.txt'], 'a');
            fprintf(fid, 'ended=1\n');
            fclose(fid);
            orphaned = [];
            evalc('orphaned = result_writers(''orphaned'', flagFn);');
            testCase.verifyTrue(orphaned, 'A writer whose record says it has ended counted as running.');
            testCase.verifyEqual(slurm.callCount('squeue'), 2);
        end

        function failedSqueueIsKeptAsRunning(testCase)
            % a failed squeue (the test Slurm stubs fail) counts as running,
            % and is kept as that for a minute; then squeue is asked again
            testCase.assumeFalse(ispc, 'Slurm commands are stubbed on Linux and macOS only.');
            [resultFn, flagRoot] = testCase.makeResultFolder();
            testCase.writeRecord(flagRoot, resultFn, 'u1', '777_3:node1:11');
            flagFn = [flagRoot, 'result_u1/blocks_1_5.mat'];
            testCase.verifyFalse(result_writers('orphaned', flagFn));
            % (now squeue answers, without the job)
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.verifyFalse(result_writers('orphaned', flagFn), 'A failed squeue was not kept as running.');
            testCase.verifyEqual(slurm.callCount('squeue'), 0);
            orphaned = [];
            evalc('orphaned = result_writers(''orphaned'', flagFn, 0);');
            testCase.verifyTrue(orphaned);
        end

        function noWaitWithoutAnotherWriter(testCase)
            [resultFn, flagRoot] = testCase.makeResultFolder();
            flagPath = [flagRoot, 'result_aaa/'];
            mkdir(flagPath);
            result_writers('register', resultFn, 'aaa', flagPath);
            testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', [], 0));
            testCase.verifyFalse(result_writers('wait', resultFn, 'aaa'));
            txt = fileread([flagPath, 'writer.txt']);
            testCase.verifySubstring(txt, sprintf('result=%s\nuuid=aaa\nowner=', resultFn));
        end
    end

    methods (Access = private)
        function [resultFn, flagRoot] = makeResultFolder(testCase)
            folder = testCase.createTempFolder();
            resultFn = [folder, 'result.zarr'];
            flagRoot = [folder, 'zarr_flag/'];
            mkdir(flagRoot);
        end

        function [] = writeRecord(~, flagRoot, resultFn, uuid, owner, age, extra)
            % the record of uuid's writer, as result_writers writes it,
            % dated age seconds from now (optional), with extra lines
            % (optional, e.g. 'ended=1')
            [~, fsn] = fileparts(resultFn);
            flagPath = sprintf('%s%s_%s/', flagRoot, fsn, uuid);
            mkdir(flagPath);
            fid = fopen([flagPath, 'writer.txt'], 'w');
            fprintf(fid, 'result=%s\nuuid=%s\nowner=%s\n', resultFn, uuid, owner);
            if nargin > 6 && ~isempty(extra)
                fprintf(fid, '%s\n', extra);
            end
            fclose(fid);
            if nargin > 5 && age ~= 0
                fileObj = java.io.File([flagPath, 'writer.txt']);
                fileObj.setLastModified(fileObj.lastModified() + age * 1000);
            end
        end

        function [pid] = startProcess(testCase)
            % a process that runs during the test (killed at the end)
            [~, out] = system('sleep 120 > /dev/null 2>&1 & echo $!');
            pid = str2double(strtrim(out));
            testCase.addTeardown(@() system(sprintf('kill %d 2> /dev/null', pid)));
        end

        function [start] = processStart(~, pid)
            % a process's start time (the 22nd field of /proc/<pid>/stat)
            txt = fileread(sprintf('/proc/%d/stat', pid));
            fields = strsplit(strtrim(txt(find(txt == ')', 1, 'last') + 1 : end)));
            start = fields{20};
        end

        function [node] = nodeName(~)
            node = getenv('SLURMD_NODENAME');
            if isempty(node)
                [~, node] = system('uname -n');
                node = strtrim(node);
            end
        end
    end
end
