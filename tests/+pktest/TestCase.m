classdef (SharedTestFixtures = {pktest.PathFixture, pktest.NoSlurmFixture}) TestCase < matlab.unittest.TestCase
    % base class for PetaKit5D tests. It provides
    %   PetaKit5D on the path (shared pktest.PathFixture)
    %   no Slurm, so tests never submit jobs even on the cluster and run
    %   the same way as on CI (shared pktest.NoSlurmFixture)
    %   a seeded random number generator for every test
    %   restoring the current folder and closing figures opened by the test
    %   a check that the mex functions the test reaches are used: a test
    %   fails if the MATLAB code ran instead of a mex this platform has
    %   (mex_fallback), unless it calls expectMexFallback
    %   helpers for temporary folders, platform/GPU assumptions, and
    %   tolerant image comparison with informative diagnostics

    properties (Access = private)
        FiguresBefore
        ExpectedMexFallbacks = {}
    end

    methods (TestMethodSetup)
        function seedRandomNumberGenerator(testCase)
            prevRng = rng(0, 'twister');
            testCase.addTeardown(@rng, prevRng);
        end

        function restoreCurrentFolder(testCase)
            % the job framework cd's to the repo root as a side effect
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(pwd));
        end

        function restoreWarningState(testCase)
            % getImageSize and getZarrInfo call warning('off', 'all') and
            % never turn warnings back on
            warningState = warning;
            testCase.addTeardown(@warning, warningState);
        end

        function clearPersistentCaches(testCase)
            % the decon functions cache the OTF in persistent variables
            % keyed only on size and sum of squares of the PSF
            clear decon_lucy_function decon_lucy_omw_function
        end

        function closeFiguresOpenedByTest(testCase)
            testCase.FiguresBefore = findall(groot, 'Type', 'figure');
            testCase.addTeardown(@() close(setdiff(findall(groot, 'Type', 'figure'), testCase.FiguresBefore), 'force'));
        end

        function checkMexFallbacks(testCase)
            % code that falls back to MATLAB code when a mex fails records it
            % with mex_fallback; outputs alone can't tell, since the MATLAB
            % code gives the same results (only slower). Fallbacks since the
            % last test (in a class setup or teardown) fail this test.
            testCase.ExpectedMexFallbacks = {};
            testCase.verifyNoUnexpectedMexFallback(' outside the tests (in a class setup or teardown)');
            testCase.addTeardown(@() testCase.verifyNoUnexpectedMexFallback(''));
        end
    end

    methods (Access = protected)
        function [dataPath] = createTempFolder(testCase)
            % new folder that is deleted after the test, with forward
            % slashes and a trailing slash as in the demos on all OSes
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            dataPath = [strrep(fixture.Folder, '\', '/'), '/'];
        end

        function [] = assumeMexAvailable(testCase, varargin)
            % skip (instead of fail) if a mex is not built for this
            % platform. MexInventoryTest tracks the missing ones.
            for i = 1 : numel(varargin)
                testCase.assumeEqual(exist(varargin{i}, 'file'), 3, ...
                    sprintf('%s.%s is not available on this platform.', varargin{i}, mexext));
            end
        end

        function [] = expectMexFallback(testCase, varargin)
            % the test runs the MATLAB code instead of these mex functions
            % on purpose (e.g., it hides a mex with pktest.FunctionStubFixture
            % or makes it fail), so their fallbacks don't fail it
            testCase.ExpectedMexFallbacks = [testCase.ExpectedMexFallbacks, varargin];
        end

        function [] = assumeNoKnownIssueOn(testCase, platforms, description)
            % skip on platforms (mexext values, e.g. 'mexmaca64' or
            % 'mexw64') where the test fails because of a known bug that
            % only affects them. List it in tests/README.md and remove the
            % platform when the bug is fixed.
            testCase.assumeFalse(ismember(mexext, cellstr(platforms)), ...
                sprintf('KNOWN ISSUE on %s: %s', mexext, description));
        end

        function [] = assumeGPUAvailable(testCase)
            hasGPU = false;
            try
                hasGPU = gpuDeviceCount('available') > 0;
            catch
            end
            testCase.assumeTrue(hasGPU, 'No CUDA GPU available.');
        end

        function [] = verifyImagesClose(testCase, actual, expected, varargin)
            % compare images with relative L2 error (robust to the rare
            % ULP or +/-1 rounding differences across platforms) and,
            % optionally, max absolute error.
            %
            %   relTol : a number (default: 1e-6). Relative L2 error tolerance.
            %   absTol : a number (default: inf). Max absolute error tolerance.
            %     name : char (default: 'image'). Name used in diagnostics.

            ip = inputParser;
            ip.addParameter('relTol', 1e-6, @isnumeric);
            ip.addParameter('absTol', inf, @isnumeric);
            ip.addParameter('name', 'image', @ischar);
            ip.parse(varargin{:});
            pr = ip.Results;

            testCase.verifyEqual(size(actual), size(expected), sprintf('%s has the wrong size.', pr.name));
            if ~isequal(size(actual), size(expected)) || isempty(expected)
                return;
            end

            a = double(actual(:));
            e = double(expected(:));
            relErr = norm(a - e) / max(norm(e), eps);
            [maxErr, ind] = max(abs(a - e));
            sub = cell(1, ndims(expected));
            [sub{:}] = ind2sub(size(expected), ind);
            diagStr = sprintf(['%s: relative L2 error %.3g (tolerance %.3g); max abs error %.4g ', ...
                '(tolerance %.4g) at %s, actual %.6g vs expected %.6g.'], pr.name, relErr, ...
                pr.relTol, maxErr, pr.absTol, mat2str([sub{:}]), a(ind), e(ind));

            testCase.verifyLessThanOrEqual(relErr, pr.relTol, diagStr);
            if isfinite(pr.absTol)
                testCase.verifyLessThanOrEqual(maxErr, pr.absTol, diagStr);
            end
        end
    end

    methods (Access = private)
        function [] = verifyNoUnexpectedMexFallback(testCase, where)
            % fallbacks of a mex this platform has (or should have: not a
            % known gap), except designed ones and those the test expects;
            % then clears the records
            r = pktest.unexpectedMexFallbacks(mex_fallback('-log'), testCase.ExpectedMexFallbacks);
            mex_fallback('-reset');
            lines = arrayfun(@(x) sprintf('%s (%s, %d times): %s %s', x.mex, x.kind, x.count, x.identifier, x.message), ...
                r, 'unif', 0);
            testCase.verifyEmpty(lines, sprintf(['The MATLAB code ran instead of mex functions this platform has%s:\n  %s\n', ...
                'If the test runs the MATLAB code on purpose, call testCase.expectMexFallback with the mex names.'], ...
                where, strjoin(lines, '\n  ')));
        end
    end
end
