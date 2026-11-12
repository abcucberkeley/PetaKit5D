classdef WaitForFilesTest < pktest.TestCase
    % wait_for_files, the wait of the callers of single jobs between their
    % rounds: at most waitTime, until all files of one group exist. pause is
    % a spy (it doesn't wait); files appear at its calls (pktest.atCall).

    properties (TestParameter)
        % how long a look takes (ms): slow, and as slow as the first look
        % on a loaded Windows runner (0.95 s: the next pause is the rest of
        % the wait)
        lookTime = struct('slow', 300, 'verySlow', 950)
    end

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
            % (1, 1, 1, 0.5 s; a look over 0.1 s, as on a slow machine,
            % stretches the next interval)
            waits = testCase.waits(pauses);
            testCase.verifyLooksAtLeastASecondApart(waits);
            testCase.verifyEqual(sum(waits), 3.5, 'AbsTol', 1e-9);
        end

        function waitEndsOnceTheFilesOfAGroupExist(testCase)
            % groups {a, b} and {c}: b comes during the second second, a
            % during the third
            folder = testCase.createTempFolder();
            fns = strcat(folder, {'a.tif', 'b.tif', 'c.tif'});
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [folder, 'calls.txt'], [2, 3], {@() fclose(fopen(fns{2}, 'w')), @() fclose(fopen(fns{1}, 'w'))})));
            testCase.verifyTrue(wait_for_files({fns(1 : 2), fns(3)}, 30));
            testCase.verifyLooksAtLeastASecondApart(testCase.waits(pauses), 3);
        end

        function waitEndsOnceAllFilesOfOneGroupExist(testCase)
            % one group (the callers' results): a during the second second,
            % b during the third
            folder = testCase.createTempFolder();
            fns = strcat(folder, {'a.tif', 'b.tif'});
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [folder, 'calls.txt'], [2, 3], {@() fclose(fopen(fns{1}, 'w')), @() fclose(fopen(fns{2}, 'w'))})));
            testCase.verifyTrue(wait_for_files({fns}, 30));
            testCase.verifyLooksAtLeastASecondApart(testCase.waits(pauses), 3);
        end

        function resultsCanBeFolders(testCase)
            % (zarr results)
            folder = testCase.createTempFolder();
            resultFn = [folder, 'result.zarr'];
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [folder, 'calls.txt'], 2, @() mkdir(resultFn))));
            testCase.verifyTrue(wait_for_files({{resultFn}}, 30));
            testCase.verifyLooksAtLeastASecondApart(testCase.waits(pauses), 2);
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
            testCase.verifyLooksAtLeastASecondApart(testCase.waits(pauses), 3);
        end

        function slowLooksAreFewer(testCase, lookTime)
            % a look that takes 0.3 s or more (slow storage, many files):
            % the next ones are 10 times as far apart (3 s or more), so that
            % looking takes at most a tenth of the wait, which stays the
            % same. The stub records how long each look took, so a slower
            % machine (longer looks, fewer of them) passes too.
            testCase.assumeTrue(usejava('jvm'), 'The slow stub sleeps with Java.');
            stubFolder = testCase.createTempFolder();
            fid = fopen([stubFolder, 'batch_file_exist.m'], 'w');
            fprintf(fid, '%s\n', 'function [tf] = batch_file_exist(paths, varargin)', 't0 = tic;', ...
                sprintf('java.lang.Thread.sleep(%d);', lookTime), ...
                'fid = fopen([fileparts(mfilename(''fullpath'')), ''/looks.txt''], ''a'');', ...
                'fprintf(fid, ''%.6f\n'', toc(t0));', 'fclose(fid);', 'tf = false(numel(paths), 1);', 'end');
            fclose(fid);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(stubFolder));
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause'));
            testCase.verifyFalse(wait_for_files({{[stubFolder, 'a.tif']}}, 10));
            waits = testCase.waits(pauses);
            fid = fopen([stubFolder, 'looks.txt'], 'r');
            looks = fscanf(fid, '%f');
            fclose(fid);
            % a look after each pause but the last
            testCase.assertNumElements(looks, numel(waits) - 1);
            testCase.verifyEqual(waits(1), 1);
            for k = 1 : numel(looks)
                % at least 10 times the look before it, or the rest of the wait
                left = 10 - sum(waits(1 : k));
                testCase.verifyGreaterThanOrEqual(waits(k + 1), min(10 * looks(k), left) - 1e-9, ...
                    sprintf('The pause after look %d (%.3f s) is %.3f s.', k, looks(k), waits(k + 1)));
            end
            % (looks of under 0.1 s: 9)
            testCase.verifyLessThanOrEqual(numel(looks), 3);
            testCase.verifyEqual(sum(waits), 10, 'AbsTol', 1e-9);
        end
    end

    methods (Access = private)
        function [w] = waits(~, pauses)
            % the duration of each pause call
            w = cellfun(@(x) x{1}, pauses.getCalls());
        end

        function verifyLooksAtLeastASecondApart(testCase, waits, n)
            % a look after each second or more (more after a look that took
            % over 0.1 s, as on a slow machine): the first pause is a second,
            % none but the last is shorter; n: the number of pauses (the
            % files come at the pause calls of the test, whatever their length)
            if nargin > 2
                testCase.verifyNumElements(waits, n);
            end
            testCase.verifyEqual(waits(1), 1);
            testCase.verifyGreaterThanOrEqual(waits(1 : end - 1), 1 - 1e-12);
        end
    end
end
