classdef IsolationHarnessTest < pktest.TestCase
    % the isolation harness of tests/isolation (run_isolation_harness.sh,
    % speed/reports/isolation.md): its snapshot sees what a task changes in
    % the process, its resets undo it, it finds tasks whose outputs depend on
    % what ran before them in the process, and its comparison finds any
    % difference in the outputs

    methods (TestClassSetup)
        function addHarnessToPath(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(mfilename('fullpath'))));
        end
    end

    methods (Test, TestTags = {'Unit'})
        function snapshotSeesWhatATaskChanges(testCase)
            st0 = isolation_state();
            testCase.changeProcessState();
            changes = strjoin(isolation_state_diff(st0, isolation_state()), newline);
            for item = {'threads', 'warnings: isolation:test', 'rng', 'pwd', 'path', ...
                    'env: ISOLATION_TEST_VAR', 'globals: isolationTestGlobal', 'open files'}
                testCase.verifySubstring(changes, item{1});
            end
        end

        function fileDescriptorsAreThoseOfFilesOutsideMATLAB(testCase)
            % only descriptors of files outside MATLAB count, as a task can
            % leave them open: MATLAB's own (sockets, pipes, event and device
            % descriptors) come and go in the background, and their count
            % changed by one between two looks (fullResetRestoresTheProcess
            % on GitHub's R2026b Linux runner). A device descriptor (here
            % /dev/null) doesn't count, a file left open does.
            testCase.assumeTrue(isunix && ~ismac, 'File descriptors are counted on Linux only.');
            st0 = isolation_state();
            fidDev = fopen('/dev/null', 'r');
            testCase.addTeardown(@() fclose(fidDev));
            st1 = isolation_state();
            testCase.verifyEqual(st1.fds, st0.fds, 'A device descriptor was counted.');
            fid = fopen([testCase.createTempFolder(), 'left_open.txt'], 'w');
            testCase.addTeardown(@() fclose(fid));
            st2 = isolation_state();
            testCase.verifyEqual(st2.fds - st1.fds, 1, 'A file left open was not counted.');
            testCase.verifySubstring(strjoin(isolation_state_diff(st1, st2), newline), 'file descriptors (of files outside MATLAB)');
        end

        function fullResetRestoresTheProcess(testCase)
            % (a new process has no last warning)
            lastwarn('', '');
            st0 = isolation_state();
            testCase.changeProcessState();
            isolation_reset('full', st0);
            leftover = isolation_state_diff(st0, isolation_state());
            testCase.verifyEmpty(leftover(~startsWith(leftover, '(info)')), strjoin(leftover, newline));
        end

        function prototypeResetRestoresItsPart(testCase)
            st0 = isolation_state();
            testCase.changeProcessState();
            isolation_reset('prototype', st0);
            leftover = strjoin(isolation_state_diff(st0, isolation_state()), newline);
            for item = {'warnings', 'rng', 'pwd', 'open files'}
                testCase.verifyFalse(contains(leftover, item{1}), sprintf('%s was not restored: %s', item{1}, leftover));
            end
            % the rest is left
            testCase.verifySubstring(leftover, 'threads');
            testCase.verifySubstring(leftover, 'env: ISOLATION_TEST_VAR');
        end

        function runRecordsLeaksAndResetsUndoThem(testCase)
            % two runs of a task that saves random numbers of the global
            % stream: in one process without a reset, the second saves other
            % numbers than the first; with the prototype's reset, the same
            dataPath = testCase.createTempFolder();
            for reset = {'none', 'prototype'}
                rng('default');
                tasks = testCase.exampleTasks(dataPath, reset{1}, {'rng', 'rng'});
                records = isolation_run_tasks(tasks, 'reset', reset{1});
                a = load(tasks(1).output);
                b = load(tasks(2).output);
                testCase.verifyEqual(isequal(a.value, b.value), strcmp(reset{1}, 'prototype'), reset{1});
                testCase.verifySubstring(strjoin(records(1).changes, newline), 'rng');
                testCase.verifyTrue(all([records.ok]));
            end
        end

        function fullResetClearsPersistentVariables(testCase)
            dataPath = testCase.createTempFolder();
            for reset = {'prototype', 'full'}
                clear isolation_example_task
                tasks = testCase.exampleTasks(dataPath, reset{1}, {'persistent', 'persistent'});
                isolation_run_tasks(tasks, 'reset', reset{1});
                a = load(tasks(1).output);
                b = load(tasks(2).output);
                testCase.verifyEqual(isequal(a.value, b.value), strcmp(reset{1}, 'full'), reset{1});
            end
        end

        function compareFindsDifferentOutputs(testCase)
            root = testCase.createTempFolder();
            a = [root, 'a/'];
            b = [root, 'b/'];
            im = uint16(reshape(1 : 4000, 10, 20, 20));
            for d = {a, b}
                mkdir([d{1}, 'sub']);
                writetiff(im, [d{1}, 'same.tif']);
                writezarr(im, [d{1}, 'same.zarr'], 'blockSize', [8, 8, 8]);
                value = 1;
                save([d{1}, 'sub/same.mat'], 'value');
            end
            % the same data in other chunks
            writezarr(im, [a, 'chunks.zarr'], 'blockSize', [8, 8, 8]);
            writezarr(im, [b, 'chunks.zarr'], 'blockSize', [16, 16, 16]);
            im2 = im;
            im2(5) = 0;
            writetiff(im, [a, 'diff.tif']);
            writetiff(im2, [b, 'diff.tif']);
            writezarr(im, [a, 'diff.zarr'], 'blockSize', [8, 8, 8]);
            writezarr(im2, [b, 'diff.zarr'], 'blockSize', [8, 8, 8]);
            value = 1;
            save([a, 'diff.mat'], 'value');
            value = 2;
            save([b, 'diff.mat'], 'value');
            % uuid names are matched
            uuids = {'0d8c8a76-1b6e-4c3e-9a51-0f5f2b6c1e11', 'ffa1b2c3-1111-2222-3333-444455556666'};
            fclose(fopen(sprintf('%stmp_%s.txt', a, uuids{1}), 'w'));
            fclose(fopen(sprintf('%stmp_%s.txt', b, uuids{2}), 'w'));
            fclose(fopen([a, 'onlyA.txt'], 'w'));

            res = isolation_compare(a, b);
            status = containers.Map({res.path}, {res.status});
            testCase.verifyEqual(status('same.tif'), 'equal');
            testCase.verifyEqual(status('same.zarr'), 'equal');
            testCase.verifyEqual(status('sub/same.mat'), 'equal');
            testCase.verifyEqual(status('chunks.zarr'), 'differ', 'Zarr metadata (chunks) is compared.');
            testCase.verifyEqual(status('diff.tif'), 'differ');
            testCase.verifyEqual(status('diff.zarr'), 'differ');
            testCase.verifyEqual(status('diff.mat'), 'differ');
            testCase.verifyEqual(status('tmp_<uuid>.txt'), 'equal');
            testCase.verifyEqual(status('onlyA.txt'), 'onlyA');
            detail = res(strcmp({res.path}, 'diff.tif')).detail;
            testCase.verifySubstring(detail, '1 of 4000 voxels differ');
        end
    end

    methods (Test, TestTags = {'Integration'})
        function capturedTasksGiveTheSameOutputsInAnyOrder(testCase)
            % the crop and MIP scenarios: their folders are left as they are
            % before the tasks, and the tasks give the same outputs one per
            % process (here: after a full reset) and in one process in
            % reverse order
            root = testCase.createTempFolder();
            work = [root, 'work/'];
            [tasks, scenarios] = evalc_outputs(@() isolation_prepare(work, 'scenarios', {'crop', 'mip'}));
            testCase.assertEqual({scenarios.name}, {'crop', 'mip'});
            testCase.assertEqual(numel(tasks), 4);
            testCase.verifyFalse(any(cellfun(@(x) exist(x, 'file') > 0, {tasks.output})), ...
                'Outputs exist before the tasks ran.');
            copyfile(work, [root, 'template/']);

            evalc('isolation_run_tasks(tasks, ''order'', 1 : 4, ''reset'', ''full'')');
            copyfile(work, [root, 'fwd/']);
            rmdir(work, 's');
            copyfile([root, 'template/'], work);
            evalc('isolation_run_tasks(tasks, ''order'', 4 : -1 : 1, ''reset'', ''none'')');

            testCase.verifyTrue(all(cellfun(@(x) exist(x, 'file') > 0, {tasks.output})), 'Outputs are missing.');
            res = isolation_compare([root, 'fwd/'], work);
            testCase.verifyTrue(all(strcmp({res.status}, 'equal')), strjoin(strcat({res.path}, {': '}, {res.status}), newline));
        end
    end

    methods (Access = private)
        function [] = changeProcessState(testCase)
            % change each part of the process state the resets restore (the
            % base class and the teardowns restore it as well)
            n0 = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(n0));
            maxNumCompThreads(max(1, n0 - 1) + (n0 == 1));
            warning('on', 'isolation:test');
            warning('off', 'isolation:test');
            rand(10, 1);
            cd(tempdir);
            pathEntry = testCase.createTempFolder();
            p = path;
            testCase.addTeardown(@() path(p));
            addpath(pathEntry, '-end');
            testCase.addTeardown(@() setenv('ISOLATION_TEST_VAR', ''));
            setenv('ISOLATION_TEST_VAR', 'x');
            testCase.addTeardown(@() clearvars('-global', 'isolationTestGlobal'));
            global isolationTestGlobal %#ok<GVMIS>
            isolationTestGlobal = 1;
            fn = [testCase.createTempFolder(), 'open.txt'];
            fid = fopen(fn, 'w');
            testCase.addTeardown(@() fclose_if_open(fid));
        end

        function [tasks] = exampleTasks(~, dataPath, name, leaks)
            tasks = struct('id', {}, 'scenario', {}, 'func', {}, 'funcStr', {}, 'output', {}, 'gpu', {}, 'threads', {});
            for i = 1 : numel(leaks)
                t.id = sprintf('%s_%d', name, i);
                t.scenario = 'example';
                t.func = 'isolation_example_task';
                t.output = sprintf('%s%s_%d.mat', dataPath, name, i);
                t.funcStr = sprintf('isolation_example_task(''%s'',''%s'')', t.output, leaks{i});
                t.gpu = false;
                t.threads = [];
                tasks(end + 1) = t; %#ok<AGROW>
            end
        end
    end
end


function [] = fclose_if_open(fid)

if any(isolation_open_fids() == fid)
    fclose(fid);
end

end


function [varargout] = evalc_outputs(fcn)
% call fcn without printing its output

varargout = cell(1, nargout);
evalc('[varargout{:}] = fcn();');

end
