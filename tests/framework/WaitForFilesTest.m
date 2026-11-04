classdef WaitForFilesTest < pktest.TestCase
    % wait_for_files, the wait of the callers of single jobs between their
    % rounds: at most waitTime, until all files of one group exist. pause is
    % a spy (it doesn't wait); files appear at its calls (pktest.atCall).

    methods (Test, TestTags = {'Unit'})
        function withoutFilesItIsOnePause(testCase)
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause'));
            testCase.verifyFalse(wait_for_files({}, 4.5));
            testCase.verifyFalse(wait_for_files({{}, {}}, 3));
            testCase.verifyEqual(testCase.waits(pauses), [4.5, 3]);
        end

        function filesThatDoNotComeTakeTheWholeWait(testCase)
            % a look after each second, none at the end (the caller's round
            % looks then)
            folder = testCase.createTempFolder();
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause'));
            testCase.verifyFalse(wait_for_files({{[folder, 'a.tif']}}, 3.5));
            testCase.verifyEqual(testCase.waits(pauses), [1, 1, 1, 0.5], 'AbsTol', 1e-12);
        end

        function waitEndsOnceTheFilesOfAGroupExist(testCase)
            % groups {a, b} and {c}: b comes during the second second, a
            % during the third
            folder = testCase.createTempFolder();
            fns = strcat(folder, {'a.tif', 'b.tif', 'c.tif'});
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [folder, 'calls.txt'], [2, 3], {@() fclose(fopen(fns{2}, 'w')), @() fclose(fopen(fns{1}, 'w'))})));
            testCase.verifyTrue(wait_for_files({fns(1 : 2), fns(3)}, 30));
            testCase.verifyEqual(testCase.waits(pauses), [1, 1, 1]);
        end

        function waitEndsOnceAllFilesOfOneGroupExist(testCase)
            % one group (the callers' results): a during the second second,
            % b during the third
            folder = testCase.createTempFolder();
            fns = strcat(folder, {'a.tif', 'b.tif'});
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [folder, 'calls.txt'], [2, 3], {@() fclose(fopen(fns{1}, 'w')), @() fclose(fopen(fns{2}, 'w'))})));
            testCase.verifyTrue(wait_for_files({fns}, 30));
            testCase.verifyEqual(testCase.waits(pauses), [1, 1, 1]);
        end

        function resultsCanBeFolders(testCase)
            % (zarr results)
            folder = testCase.createTempFolder();
            resultFn = [folder, 'result.zarr'];
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [folder, 'calls.txt'], 2, @() mkdir(resultFn))));
            testCase.verifyTrue(wait_for_files({{resultFn}}, 30));
            testCase.verifyEqual(testCase.waits(pauses), [1, 1]);
        end

        function patternsCountOnlyFinishedNames(testCase)
            % a pattern, as the pipeline looks for a stitched result: the
            % result's temporary name (with a suffix) doesn't count, its
            % final name does
            folder = testCase.createTempFolder();
            resultFn = [folder, 'Scan_Iter_0000_msec_0000000001msecAbs.zarr'];
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [folder, 'calls.txt'], [2, 3], {@() mkdir([resultFn, '_1234']), @() mkdir(resultFn)})));
            testCase.verifyTrue(wait_for_files({{[folder, 'Scan_Iter_0000_msec*Abs.zarr']}}, 30));
            testCase.verifyEqual(testCase.waits(pauses), [1, 1, 1]);
        end

        function slowLooksAreFewer(testCase)
            % a look that takes 0.3 s (slow storage, many files): the next
            % ones are 3 s or more apart, so that looking takes at most a
            % tenth of the wait, which stays the same
            testCase.assumeTrue(usejava('jvm'), 'The slow stub sleeps with Java.');
            stubFolder = testCase.createTempFolder();
            fid = fopen([stubFolder, 'batch_file_exist.m'], 'w');
            fprintf(fid, 'function [tf] = batch_file_exist(paths, varargin)\njava.lang.Thread.sleep(300);\ntf = false(numel(paths), 1);\nend\n');
            fclose(fid);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(stubFolder));
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause'));
            testCase.verifyFalse(wait_for_files({{[stubFolder, 'a.tif']}}, 10));
            waits = testCase.waits(pauses);
            testCase.verifyEqual(waits(1), 1);
            testCase.verifyGreaterThanOrEqual(waits(2 : end - 1), 3);
            testCase.verifyEqual(sum(waits), 10, 'AbsTol', 1e-9);
        end
    end

    methods (Access = private)
        function [w] = waits(~, pauses)
            % the duration of each pause call
            w = cellfun(@(x) x{1}, pauses.getCalls());
        end
    end
end
