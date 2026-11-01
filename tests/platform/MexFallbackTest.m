classdef MexFallbackTest < pktest.TestCase
    % mex_fallback, which records where the MATLAB code runs instead of a
    % mex (warning once per mex, error and session if the mex exists), and the check
    % of pktest.TestCase that fails a test in which that happened unless the
    % test expects it (expectMexFallback), a known platform gap or a
    % rejection by design. crop3d is the example: a bbox column makes
    % crop3d_mex fail ('Input range for bbox is not 6', as wrong arguments
    % would), int16 data is a type it rejects by design, and the MATLAB code
    % gives the same crop.

    properties (Constant)
        Image = reshape(uint16(1 : 6 * 5 * 4), 6, 5, 4)
        Bbox = [2, 1, 2, 5, 4, 3]
    end

    methods (Test, TestTags = {'Unit'})
        function failedMexWarnsOnceAndIsRecorded(testCase)
            testCase.assumeMexAvailable('crop3d_mex');
            testCase.expectMexFallback('crop3d_mex');
            b = testCase.Bbox;
            expected = testCase.Image(b(1) : b(4), b(2) : b(5), b(3) : b(6));
            out = testCase.verifyWarning(@() crop3d(testCase.Image, b'), 'PetaKit5D:mexFallback');
            testCase.verifyEqual(out, expected);
            % once per mex, error and session
            out = testCase.verifyWarningFree(@() crop3d(testCase.Image, b'));
            testCase.verifyEqual(out, expected);

            r = mex_fallback('-log');
            testCase.verifyEqual(numel(r), 1);
            testCase.verifyEqual({r.mex, r.kind, r.identifier, r.message, r.count}, ...
                {'crop3d_mex', 'failed', 'crop:inputError', 'Input range for bbox is not 6', 2});
        end

        function warningIsOneLineWithTheMexItsErrorAndTheCaller(testCase)
            testCase.assumeMexAvailable('crop3d_mex');
            testCase.expectMexFallback('crop3d_mex');
            lastwarn('', '');
            out = evalc('crop3d(testCase.Image, testCase.Bbox'');');
            [msg, id] = lastwarn;
            testCase.verifyEqual(id, 'PetaKit5D:mexFallback');
            testCase.verifyMatches(msg, ['^crop3d_mex failed in crop3d \(line \d+\): Input range for bbox is not 6; ', ...
                'the MATLAB code runs instead\.$']);
            % shown without the backtrace (the display wraps long lines)
            testCase.verifySubstring(out, 'crop3d_mex failed in crop3d');
            testCase.verifyFalse(contains(out, '> In '), out);
        end

        function eachMexWarnsOnce(testCase)
            testCase.expectMexFallback('crop3d_mex', 'crop4d_mex');
            testCase.assumeMexAvailable('crop3d_mex', 'crop4d_mex');
            ME = MException('pktest:failed', 'Failed on purpose.');
            testCase.verifyWarning(@() mex_fallback('crop3d_mex', ME), 'PetaKit5D:mexFallback');
            testCase.verifyWarning(@() mex_fallback('crop4d_mex', ME), 'PetaKit5D:mexFallback');
            testCase.verifyWarningFree(@() mex_fallback('crop3d_mex', ME));
            testCase.verifyWarningFree(@() mex_fallback('crop4d_mex', ME));
        end

        function eachErrorOfAMexWarnsOnce(testCase)
            % a failure that doesn't matter (e.g., a missing file) doesn't
            % hide a later failure of the same mex with another error
            testCase.expectMexFallback('crop3d_mex');
            testCase.assumeMexAvailable('crop3d_mex');
            ME1 = MException('pktest:missingFile', 'No such file.');
            ME2 = MException('pktest:loadError', 'The library did not load.');
            testCase.verifyWarning(@() mex_fallback('crop3d_mex', ME1), 'PetaKit5D:mexFallback');
            testCase.verifyWarning(@() mex_fallback('crop3d_mex', ME2), 'PetaKit5D:mexFallback');
            testCase.verifyWarningFree(@() mex_fallback('crop3d_mex', ME1));
            testCase.verifyWarningFree(@() mex_fallback('crop3d_mex', ME2));
        end

        function unusedMexWithoutAnErrorWarns(testCase)
            % code that checks that a mex exists and runs the MATLAB code
            % anyway (e.g., a wrong check) calls it without an error
            testCase.assumeMexAvailable('crop3d_mex');
            testCase.expectMexFallback('crop3d_mex');
            testCase.verifyWarning(@() mex_fallback('crop3d_mex'), 'PetaKit5D:mexFallback');
            [msg, ~] = lastwarn;
            testCase.verifySubstring(msg, 'crop3d_mex was not used');
            r = mex_fallback('-log');
            testCase.verifyEqual({r.kind, r.identifier}, {'failed', ''});
        end

        function missingMexIsRecordedWithoutWarning(testCase)
            % e.g. a mex not built for the platform: calling it errors
            testCase.expectMexFallback('pktest_no_such_mex');
            try
                pktest_no_such_mex(1);
            catch ME
            end
            testCase.verifyWarningFree(@() mex_fallback('pktest_no_such_mex', ME));
            r = mex_fallback('-log');
            testCase.verifyEqual({r.mex, r.kind, r.identifier}, {'pktest_no_such_mex', 'missing', 'MATLAB:UndefinedFunction'});
        end

        function designedRejectionIsRecordedWithoutWarning(testCase)
            testCase.assumeMexAvailable('crop3d_mex');
            im = int16(testCase.Image);
            b = testCase.Bbox;
            out = testCase.verifyWarningFree(@() crop3d(im, b));
            testCase.verifyEqual(out, im(b(1) : b(4), b(2) : b(5), b(3) : b(6)));
            r = mex_fallback('-log');
            testCase.verifyEqual({r.mex, r.kind, r.identifier}, {'crop3d_mex', 'designed', 'crop:dataTypeError'});
            % a designed rejection doesn't use up the warning of the mex
            testCase.expectMexFallback('crop3d_mex');
            testCase.verifyWarning(@() crop3d(testCase.Image, b'), 'PetaKit5D:mexFallback');
        end

        function recordsAreCountedPerMexKindAndIdentifier(testCase)
            testCase.expectMexFallback('crop3d_mex', 'pktest_no_such_mex');
            testCase.assumeMexAvailable('crop3d_mex');
            evalc('mex_fallback(''crop3d_mex'', MException(''pktest:a'', ''first a''))');
            mex_fallback('crop3d_mex', MException('pktest:a', 'second a'));
            mex_fallback('crop3d_mex', MException('pktest:b', 'b'));
            mex_fallback('crop3d_mex', MException('pktest:c', 'c'), {'pktest:c'});
            mex_fallback('pktest_no_such_mex');
            mex_fallback('pktest_no_such_mex');
            r = mex_fallback('-log');
            testCase.verifyEqual({r.kind}, {'failed', 'failed', 'designed', 'missing'});
            testCase.verifyEqual({r.identifier}, {'pktest:a', 'pktest:b', 'pktest:c', ''});
            testCase.verifyEqual([r.count], [2, 1, 1, 2]);
            % the first message of each
            testCase.verifyEqual(r(1).message, 'first a');
        end

        function resetClearsTheRecordsAndWarnings(testCase)
            testCase.assumeMexAvailable('crop3d_mex');
            testCase.expectMexFallback('crop3d_mex');
            ME = MException('pktest:failed', 'Failed on purpose.');
            testCase.verifyWarning(@() mex_fallback('crop3d_mex', ME), 'PetaKit5D:mexFallback');
            testCase.verifyNotEmpty(mex_fallback('-log'));
            mex_fallback('-reset');
            testCase.verifyEmpty(mex_fallback('-log'));
            testCase.verifyWarning(@() mex_fallback('crop3d_mex', ME), 'PetaKit5D:mexFallback');
        end

        function warningsTurnedOffDoNotUseUpTheWarning(testCase)
            % callers may turn warnings off (getImageSize hides tiff tag
            % warnings), the next fallback with warnings on still warns
            testCase.assumeMexAvailable('crop3d_mex');
            testCase.expectMexFallback('crop3d_mex');
            ME = MException('pktest:failed', 'Failed on purpose.');
            state = warning('off', 'all');
            lastwarn('', '');
            mex_fallback('crop3d_mex', ME);
            warning(state);
            testCase.verifyEqual(lastwarn, '');
            testCase.verifyWarning(@() mex_fallback('crop3d_mex', ME), 'PetaKit5D:mexFallback');
            testCase.verifyEqual([mex_fallback('-log').count], 2);
        end

        function functionsInPackagesAreFound(testCase)
            % (e.g., MATLAB's built-in resize mex of imresize3_fast)
            name = 'matlab.images.internal.resize.imresizemex';
            testCase.assumeNotEmpty(which(name), 'No imresizemex in this MATLAB.');
            testCase.expectMexFallback(name, 'matlab.images.internal.resize.pktest_no_such_mex');
            evalc('mex_fallback(name, MException(''pktest:failed'', ''Failed on purpose.''))');
            testCase.verifyEqual(mex_fallback('-log').kind, 'failed');
            mex_fallback('-reset');
            mex_fallback('matlab.images.internal.resize.pktest_no_such_mex');
            testCase.verifyEqual(mex_fallback('-log').kind, 'missing');
        end

        function unexpectedFallbacksAreThoseOfExistingAndUntrackedMissingMex(testCase)
            records = struct('mex', {'crop3d_mex', 'crop4d_mex', 'replace_nan_with_value', 'min_bbox_3d_mex', 'tile_blending_3d_mex'}, ...
                'kind', {'failed', 'designed', 'missing', 'missing', 'missing'}, 'identifier', '', 'message', '', 'count', 1);
            names = @(r) {r.mex};
            % Linux has all mex functions: a missing one is a problem
            testCase.verifyEqual(names(pktest.unexpectedMexFallbacks(records, {}, 'mexa64')), ...
                {'crop3d_mex', 'replace_nan_with_value', 'min_bbox_3d_mex', 'tile_blending_3d_mex'});
            % known gaps of Windows and ARM macOS don't count
            testCase.verifyEqual(names(pktest.unexpectedMexFallbacks(records, {}, 'mexw64')), ...
                {'crop3d_mex', 'replace_nan_with_value', 'tile_blending_3d_mex'});
            testCase.verifyEqual(names(pktest.unexpectedMexFallbacks(records, {}, 'mexmaca64')), ...
                {'crop3d_mex', 'tile_blending_3d_mex'});
            % on platforms whose builds aren't tracked, no missing mex counts
            testCase.verifyEqual(names(pktest.unexpectedMexFallbacks(records, {}, 'mexmaci64')), {'crop3d_mex'});
            % expected ones don't count
            testCase.verifyEmpty(pktest.unexpectedMexFallbacks(records, {'crop3d_mex', 'tile_blending_3d_mex'}, 'mexmaca64'));
            testCase.verifyEmpty(pktest.unexpectedMexFallbacks(records([]), {}, 'mexa64'));
        end

        function testsFailWithAnUnexpectedFallback(testCase)
            % a test class of pktest.TestCase with fallbacks, run here
            testCase.assumeMexAvailable('crop3d_mex');
            results = testCase.runEnforcementTests({});
            outcome = containers.Map({results.Name}, num2cell([results.Failed]));
            cls = 'MexFallbackEnforcementTests';
            testCase.verifyTrue(outcome([cls, '/failedMex']), 'A failed mex call did not fail the test.');
            testCase.verifyFalse(outcome([cls, '/expectedFailedMex']), 'An expected fallback failed the test.');
            testCase.verifyFalse(outcome([cls, '/designedRejection']), 'A rejection by design failed the test.');
            [~, tracked] = pktest.mexPlatformGaps();
            if tracked
                testCase.verifyTrue(outcome([cls, '/missingMex']), 'A mex missing from a tracked platform did not fail the test.');
            else
                testCase.verifyFalse(outcome([cls, '/missingMex']), 'A missing mex failed the test on an untracked platform.');
            end
            testCase.verifyFalse(outcome([cls, '/noFallback']), 'A test without fallbacks failed.');
            % the diagnostic names the mex and the opt-out
            r = results(strcmp({results.Name}, [cls, '/failedMex']));
            report = strjoin(arrayfun(@(x) x.Report, r.Details.DiagnosticRecord, 'unif', 0), newline);
            testCase.verifySubstring(report, 'crop3d_mex (failed, 2 times): crop:inputError Input range for bbox is not 6');
            testCase.verifySubstring(report, 'testCase.expectMexFallback');
        end

        function fallbackInAClassSetupFailsTheFirstTest(testCase)
            testCase.assumeMexAvailable('crop3d_mex');
            lines = {
                'classdef MexFallbackClassSetupTests < pktest.TestCase'
                '    % written by MexFallbackTest'
                '    methods (TestClassSetup)'
                '        function fallBack(testCase)'
                '            crop3d(uint16(ones(4, 4, 4)), [1; 1; 1; 2; 2; 2]);'
                '        end'
                '    end'
                '    methods (Test)'
                '        function first(testCase)'
                '        end'
                '        function second(testCase)'
                '        end'
                '    end'
                'end'};
            plugin = pktest.MexFallbackPlugin;
            results = testCase.runNestedTests('MexFallbackClassSetupTests', lines, {plugin});
            testCase.verifyEqual([results.Failed], [true, false]);
            report = strjoin(arrayfun(@(x) x.Report, results(1).Details.DiagnosticRecord, 'unif', 0), newline);
            testCase.verifySubstring(report, 'outside the tests (in a class setup or teardown)');
            % the summary of run_tests names them so
            testCase.verifyEqual({plugin.Records.test}, {sprintf('(class setup or teardown before %s)', results(1).Name)});
        end

        function pluginSummarizesTheFallbacksOfARun(testCase)
            testCase.assumeMexAvailable('crop3d_mex');
            plugin = pktest.MexFallbackPlugin;
            testCase.runEnforcementTests({plugin});
            cls = 'MexFallbackEnforcementTests';
            r = plugin.Records;
            testCase.verifyEqual(sort({r.test}), sort(strcat([cls, '/'], {'designedRejection', 'expectedFailedMex', ...
                'failedMex', 'missingMex'})));
            lines = plugin.summary();
            testCase.verifyEqual(numel(lines), 3);
            testCase.verifyTrue(any(startsWith(lines, sprintf('crop3d_mex failed (crop:inputError): 3 times in 2 tests: %s/', cls))), ...
                strjoin(lines, newline));
            testCase.verifyTrue(any(startsWith(lines, 'crop3d_mex designed (crop:dataTypeError): 1 times in 1 tests')), ...
                strjoin(lines, newline));
            testCase.verifyTrue(any(startsWith(lines, 'pktest_no_such_mex missing (MATLAB:UndefinedFunction): 1 times in 1 tests')), ...
                strjoin(lines, newline));
        end
    end

    methods (Access = private)
        function [results] = runEnforcementTests(testCase, plugins)
            % a test class with fallbacks, run with a runner of its own
            % (and plugins)
            lines = {
                'classdef MexFallbackEnforcementTests < pktest.TestCase'
                '    % written by MexFallbackTest'
                '    methods (Test)'
                '        function failedMex(testCase)'
                '            crop3d(uint16(ones(4, 4, 4)), [1; 1; 1; 2; 2; 2]);'
                '            crop3d(uint16(ones(4, 4, 4)), [1; 1; 1; 3; 3; 3]);'
                '        end'
                '        function expectedFailedMex(testCase)'
                '            testCase.expectMexFallback(''crop3d_mex'');'
                '            crop3d(uint16(ones(4, 4, 4)), [1; 1; 1; 2; 2; 2]);'
                '        end'
                '        function designedRejection(testCase)'
                '            crop3d(int16(ones(4, 4, 4)), [1, 1, 1, 2, 2, 2]);'
                '        end'
                '        function missingMex(testCase)'
                '            try'
                '                pktest_no_such_mex(1);'
                '            catch ME'
                '                mex_fallback(''pktest_no_such_mex'', ME);'
                '            end'
                '        end'
                '        function noFallback(testCase)'
                '            crop3d(uint16(ones(4, 4, 4)), [1, 1, 1, 2, 2, 2]);'
                '        end'
                '    end'
                'end'};
            results = testCase.runNestedTests('MexFallbackEnforcementTests', lines, plugins);
        end

        function [results] = runNestedTests(testCase, name, lines, plugins)
            % write a test class and run it with a runner of its own (and
            % plugins), without output
            folder = testCase.createTempFolder();
            fid = fopen([folder, name, '.m'], 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
            addpath(folder);
            testCase.addTeardown(@rmpath, folder);
            suite = matlab.unittest.TestSuite.fromFile([folder, name, '.m']);
            runner = matlab.unittest.TestRunner.withNoPlugins;
            runner.addPlugin(matlab.unittest.plugins.DiagnosticsRecordingPlugin);
            for i = 1 : numel(plugins)
                runner.addPlugin(plugins{i});
            end
            evalc('results = runner.run(suite);');
            % the fallbacks of these tests are not this test's
            mex_fallback('-reset');
        end
    end
end
