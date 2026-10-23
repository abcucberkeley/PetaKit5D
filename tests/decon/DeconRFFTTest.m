classdef DeconRFFTTest < pktest.TestCase
    % the RL and OMW iterations with half-spectrum FFTs (decon_lucy_rfft_mex
    % on the CPU) against the MATLAB loops of decon_lucy_function and
    % decon_lucy_omw_function, the choice of the path (useMex, fallback
    % without the mex), the half OTF, and the groups of lines MATLAB's FFTs
    % transform together. The mex runs the same FFTW transforms as MATLAB's
    % fftn and ifftn (on the half spectrum), so the results are equal bit for
    % bit, whether this MATLAB session runs its FFTs in parallel or not.

    properties (Constant)
        ImSize = [48, 46, 36]
        Background = 100
        PSFSigma = [1.5, 2]
        BlurSigma = sqrt(0.8 ^ 2 + [1.5, 2] .^ 2)
    end

    properties (TestParameter)
        method = {'simplified', 'omw'}
        inputClass = {'uint16', 'single', 'uint8'}
        % even and odd sizes, 2D; length 64 along y, x or z (FFTW's
        % transforms of four and of one line of 64 differ in the last bits),
        % sizes whose lines MATLAB splits among threads (64 x 230) or just
        % doesn't (64 x 210), and two rows along y
        dataSize = {[48, 46, 36], [47, 45, 33], [66, 64, 33], [64, 50, 40], [10, 23, 64], [64, 50], [30, 64], ...
            [64, 230], [64, 210], [10, 64, 23], [2, 64, 33]}
        options = {{'dampFactor', 1.2}, {'EdgeErosion', 3}, {'deconOffset', 5}, {'scaleFactor', 2}, ...
            {'bbox', [3, 4, 5, 40, 41, 30]}, {'save16bit', true}, {'Background', 100.3}, ...
            {'EdgeErosion', 2, 'scaleFactor', 1.5, 'deconOffset', 2, 'bbox', [2, 2, 2, 44, 40, 33], 'save16bit', true, 'dampFactor', 1.3}}
    end

    properties (SetAccess = private)
        PSF
        BackProjector
    end

    methods (TestClassSetup)
        function createPSFAndBackProjector(testCase)
            testCase.PSF = pktest.deconPSF([31, 31, 31], testCase.PSFSigma);
            [~, testCase.BackProjector] = evalc('omw_backprojector_generation(testCase.PSF, 0.005, false)');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function mexEqualsMatlabLoop(testCase, method, inputClass)
            testCase.assumeCpuMex();
            im = testCase.blurredBeads(testCase.ImSize);
            bg = testCase.Background;
            if strcmp(inputClass, 'uint8')
                im = im / 8;
                bg = bg / 8;
            end
            im = cast(im, inputClass);
            args = {'Background', bg, 'save16bit', false};
            ref = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', false);
            out = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', true);
            testCase.verifyClass(out, 'single');
            testCase.verifyEqual(out, ref);
        end

        function defaultIsTheMexOnTheCpu(testCase, method)
            % both methods use the CPU mex by default (it gives the same
            % results as the MATLAB loop); true or false forces a path
            testCase.assumeCpuMex();
            im = testCase.blurredBeads(testCase.ImSize);
            [useMex, canUseMex] = decon_use_rfft_mex([], method, im, false, false, 100);
            testCase.verifyTrue(useMex);
            testCase.verifyTrue(canUseMex);
            testCase.verifyTrue(decon_use_rfft_mex(true, method, im, false, false, 100));
            testCase.verifyFalse(decon_use_rfft_mex(false, method, im, false, false, 100));
        end

        function mexEqualsMatlabLoopForOtherSizes(testCase, method, dataSize)
            testCase.assumeCpuMex();
            im = testCase.blurredBeads([dataSize, ones(1, 3 - numel(dataSize))]);
            psf = testCase.PSF;
            if numel(dataSize) == 2
                % a 2D PSF for 2D data
                psf = psf(:, :, 16);
            end
            % the PSF is cropped to the data as RLdecon does
            psf = testCase.cropToSize(psf, size(im, 1 : 3));
            args = {'Background', testCase.Background, 'save16bit', false};
            ref = testCase.runDecon(method, im, psf, args{:}, 'useMex', false);
            out = testCase.runDecon(method, im, psf, args{:}, 'useMex', true);
            testCase.verifyEqual(out, ref);
        end

        function rlAccelerationIsClampedAsInTheMatlabLoop(testCase)
            % a single noise-free bead drives the acceleration factor of the
            % RL loop above 1 (to about 1.09 in iterations 12 and 14), where
            % it is clamped to 1; realistic data rarely gets there
            testCase.assumeCpuMex();
            im = uint16(pktest.beads(testCase.ImSize, [24, 23, 18], testCase.PSFSigma, 20000, testCase.Background));
            nIter = 15;
            lambdas = testCase.doubleLambdas(im, testCase.PSF, nIter);
            testCase.assertGreaterThan(max(lambdas), 1.05, 'The input no longer drives the acceleration factor above 1.');
            args = {'Background', testCase.Background, 'save16bit', false};
            ref = decon_lucy_function(im, testCase.PSF, nIter, args{:}, 'useGPU', false, 'useMex', false);
            out = decon_lucy_function(im, testCase.PSF, nIter, args{:}, 'useGPU', false, 'useMex', true);
            testCase.verifyEqual(out, ref);
        end

        function optionsAreAppliedAsInTheMatlabLoop(testCase, method, options)
            testCase.assumeCpuMex();
            im = testCase.blurredBeads(testCase.ImSize);
            args = [{'Background', testCase.Background, 'save16bit', false}, options];
            ref = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', false);
            out = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', true);
            testCase.verifyEqual(out, ref);
            % RL with single input (the background subtraction in single)
            if strcmp(method, 'simplified')
                ims = single(im) + 0.3;
                ref = testCase.runDecon(method, ims, testCase.PSF, args{:}, 'useMex', false);
                out = testCase.runDecon(method, ims, testCase.PSF, args{:}, 'useMex', true);
                testCase.verifyEqual(out, ref);
            end
        end

        function zeroIterationsGiveTheBackgroundSubtractedInput(testCase, method)
            % also the semantics of the background subtraction (in double,
            % rounded to single, as MATLAB does for single - double)
            testCase.assumeCpuMex();
            im = single(testCase.blurredBeads(testCase.ImSize)) + 0.37;
            bg = 100.3;
            [out, errMat] = testCase.runDecon(method, im, testCase.PSF, 'Background', bg, 'save16bit', false, 'useMex', true, 'NUMIT', 0);
            testCase.verifyEqual(out, max(im - bg, 0));
            testCase.verifyEqual(errMat, zeros(0, 4));
            u = uint16(im);
            out = testCase.runDecon(method, u, testCase.PSF, 'Background', bg, 'save16bit', false, 'useMex', true, 'NUMIT', 0);
            testCase.verifyEqual(out, max(single(u) - bg, 0));
        end

        function iterationOutputsMatchTheMatlabLoop(testCase)
            % err_mat and the iteration count of decon_lucy_function
            testCase.assumeCpuMex();
            im = testCase.blurredBeads(testCase.ImSize);
            for nIter = [0, 3]
                [~, errRef, kRef] = decon_lucy_function(im, testCase.PSF, nIter, 'Background', testCase.Background, ...
                    'useGPU', false, 'useMex', false);
                [~, errMex, kMex] = decon_lucy_function(im, testCase.PSF, nIter, 'Background', testCase.Background, ...
                    'useGPU', false, 'useMex', true);
                testCase.verifyEqual(errMex, errRef);
                testCase.verifyEqual(kMex, kRef);
            end
        end

        function halfOTFIsTheFirstHalfOfTheOTF(testCase, dataSize)
            % the elements up to the middle of the last dimension
            testCase.assumeCpuMex();
            psf = testCase.PSF;
            if numel(dataSize) == 2
                psf = psf(:, :, 16);
            end
            psf = testCase.cropToSize(psf, [dataSize, ones(1, 3 - numel(dataSize))]);
            psf = psf ./ sum(psf(:));
            H = decon_lucy_rfft(false, 'otf', psf, dataSize);
            otf = decon_psf2otf(psf, dataSize);
            if numel(dataSize) == 3
                otf = otf(:, :, 1 : floor(dataSize(3) / 2) + 1);
            else
                otf = otf(:, 1 : floor(dataSize(2) / 2) + 1);
            end
            testCase.verifyEqual(H, otf);
            testCase.verifyError(@() decon_lucy_rfft(false, 'otf', psf, size(psf) - 1), 'decon_lucy_rfft_mex:outSizeIsSmallerThanPsfSize');
        end

        function resultsDoNotDependOnTheThreadCount(testCase, method)
            % as MATLAB's FFTs: the groups of lines don't depend on the threads
            testCase.assumeCpuMex();
            im = testCase.blurredBeads([66, 64, 33]);
            psf = testCase.PSF ./ sum(testCase.PSF(:));
            psfB = testCase.BackProjector ./ sum(testCase.BackProjector(:));
            for grain = [32, 0]
                out = cell(1, 2);
                nThreads = [1, 4];
                for i = 1 : 2
                    H = decon_lucy_rfft_mex('otf', psf, size(im), nThreads(i), grain);
                    switch method
                        case 'simplified'
                            Y = max(single(im) - 100, 0);
                            out{i} = decon_lucy_rfft_mex('rlstep', Y + 1, Y, H, nThreads(i), grain);
                        case 'omw'
                            Hb = decon_lucy_rfft_mex('otf', psfB, size(im), nThreads(i), grain);
                            out{i} = decon_lucy_rfft_mex('omw', im, H, Hb, 2, testCase.Background, nThreads(i), grain);
                    end
                end
                testCase.verifyEqual(out{2}, out{1});
            end
        end

        function linesAreGroupedAsMatlabDoes(testCase)
            % 'fft1' transforms the columns as MATLAB's fft(X, [], 1) does: in
            % its parallel mode, MATLAB halves the columns until a piece has at
            % most 32, then runs each piece in groups of four from its start
            % and the rest alone; without it (or when n * N * (bits of n) <=
            % 98304), the whole range so. MATLAB's own fft of a piece of at
            % most 32 columns, or of four or one, runs as one group (or alone),
            % which gives the reference.
            testCase.assumeCpuMex();
            for N = [2000, 222, 210]
                X = rand(64, N, 'single');
                pieces = testCase.halve(0, N, 32);
                if 64 * N * 7 <= 98304
                    pieces = [0, N];
                end
                refSplit = testCase.fftByPieces(X, pieces);
                refSerial = testCase.fftByPieces(X, [0, N]);
                if N == 2000
                    testCase.assumeFalse(isequal(refSplit, refSerial), ...
                        'FFTW gives the same bits for groups of four and single lines of length 64 on this CPU.');
                end
                testCase.verifyEqual(decon_lucy_rfft_mex('fft1', X, 4, 32), refSplit, sprintf('%d columns, split', N));
                testCase.verifyEqual(decon_lucy_rfft_mex('fft1', X, 4, 0), refSerial, sprintf('%d columns, no split', N));
                % this session's MATLAB FFT equals one of them
                Y = fft(X, [], 1);
                Y = Y(1 : 33, :);
                testCase.verifyTrue(isequal(Y, refSplit) || isequal(Y, refSerial));
            end
        end

        function withoutTheMexTheMatlabLoopRuns(testCase, method)
            % a .m file of the same name first on the path hides the mex
            % (exist(..., 'file') is 2), as on a platform without the mex
            testCase.applyFixture(pktest.FunctionStubFixture('decon_lucy_rfft_mex'));
            testCase.verifyEqual(exist('decon_lucy_rfft_mex', 'file'), 2);
            im = testCase.blurredBeads(testCase.ImSize);
            args = {'Background', testCase.Background, 'save16bit', false};
            out = testCase.runDecon(method, im, testCase.PSF, args{:});
            clear decon_lucy_function decon_lucy_omw_function
            ref = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', false);
            testCase.verifyEqual(out, ref);
            [~, canUseMex] = decon_use_rfft_mex([], method, im, false, false, testCase.Background);
            testCase.verifyFalse(canUseMex);
            testCase.verifyError(@() testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', true), ?MException);
        end

        function unsupportedCallsRunTheMatlabLoop(testCase, method)
            % double data, debug output, a non-scalar background, a single
            % row along y or x
            testCase.assumeCpuMex();
            im = testCase.blurredBeads(testCase.ImSize);
            [~, canDouble] = decon_use_rfft_mex([], method, double(im), false, false, 100);
            [~, canDebug] = decon_use_rfft_mex([], method, im, false, true, 100);
            [~, canVoxelBackground] = decon_use_rfft_mex([], method, im, false, false, [100, 100]);
            [~, canRow] = decon_use_rfft_mex([], method, im(1, :, :), false, false, 100);
            [~, canColumn] = decon_use_rfft_mex([], method, im(:, 1, :), false, false, 100);
            [~, canUint16] = decon_use_rfft_mex([], method, im, false, false, 100);
            testCase.verifyEqual([canDouble, canDebug, canVoxelBackground, canRow, canColumn, canUint16], ...
                [false, false, false, false, false, true]);

            args = {'Background', testCase.Background, 'save16bit', false};
            out = testCase.runDecon(method, double(im), testCase.PSF, args{:});
            clear decon_lucy_function decon_lucy_omw_function
            testCase.verifyEqual(out, testCase.runDecon(method, double(im), testCase.PSF, args{:}, 'useMex', false));
            testCase.verifyError(@() testCase.runDecon(method, double(im), testCase.PSF, args{:}, 'useMex', true), ?MException);
            % a background per voxel
            bg = 100 * ones(size(im), 'single');
            out = testCase.runDecon(method, im, testCase.PSF, 'Background', bg, 'save16bit', false);
            clear decon_lucy_function decon_lucy_omw_function
            testCase.verifyEqual(out, testCase.runDecon(method, im, testCase.PSF, 'Background', bg, 'save16bit', false, 'useMex', false));
            % a single row along y
            row = testCase.blurredBeads([1, 64, 33]);
            psf = testCase.cropToSize(testCase.PSF, [1, 64, 33]);
            out = testCase.runDecon(method, row, psf, args{:});
            clear decon_lucy_function decon_lucy_omw_function
            testCase.verifyEqual(out, testCase.runDecon(method, row, psf, args{:}, 'useMex', false));
            testCase.verifyError(@() testCase.runDecon(method, row, psf, args{:}, 'useMex', true), ?MException);
        end

        function switchingPathsRecomputesTheOTF(testCase, method)
            % the cached OTF is a half OTF for the mex and the full OTF for
            % the MATLAB loop
            testCase.assumeCpuMex();
            im = testCase.blurredBeads(testCase.ImSize);
            args = {'Background', testCase.Background, 'save16bit', false};
            mex1 = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', true);
            loop1 = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', false);
            mex2 = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', true);
            clear decon_lucy_function decon_lucy_omw_function
            loop2 = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', false);
            testCase.verifyEqual(mex2, mex1);
            testCase.verifyEqual(loop1, loop2);
            testCase.verifyEqual(mex1, loop1);
        end

        function mexRejectsInvalidInputs(testCase)
            testCase.assumeCpuMex();
            im = testCase.blurredBeads([16, 12, 10]);
            Y = single(im);
            H = decon_lucy_rfft_mex('otf', single(ones(3, 3, 3) / 27), [16, 12, 10], 1, 32);
            testCase.verifyEqual(size(H), [16, 12, 6]);
            id = 'decon_lucy_rfft_mex:inputError';
            testCase.verifyError(@() decon_lucy_rfft_mex('omw', im, H(:, :, 1 : end - 1), H, 2, 100, 1, 32), id);
            testCase.verifyError(@() decon_lucy_rfft_mex('omw', double(im), H, H, 2, 100, 1, 32), id);
            testCase.verifyError(@() decon_lucy_rfft_mex('omw', im, real(H), H, 2, 100, 1, 32), id);
            testCase.verifyError(@() decon_lucy_rfft_mex('omw', im, H, H, 2, 100, 1), id);
            testCase.verifyError(@() decon_lucy_rfft_mex('rlstep', Y, im, H, 1, 32), id);
            testCase.verifyError(@() decon_lucy_rfft_mex('rlstep', Y, Y(:, :, 1 : end - 1), H, 1, 32), id);
            testCase.verifyError(@() decon_lucy_rfft_mex('rlstep', Y(1, :, :), Y(1, :, :), H, 1, 32), id);
            testCase.verifyError(@() decon_lucy_rfft_mex('otf', single(1), [1, 12, 10], 1, 32), id);
            testCase.verifyError(@() decon_lucy_rfft_mex('fft1', Y, 1, 32), id);
            testCase.verifyError(@() decon_lucy_rfft_mex('unknown', im, H, 2, 100, 1, 32), id);
            testCase.verifyError(@() decon_lucy_rfft_mex('omw', zeros(4, 4, 4, 2, 'uint16'), H, H, 2, 100, 1, 32), id);
        end
    end

    methods (Access = private)
        function [] = assumeCpuMex(testCase)
            % the mex gives the same bits as MATLAB's FFTs on Linux, where it
            % is the default (on macOS they differ in the last bits)
            testCase.assumeMexAvailable('decon_lucy_rfft_mex');
            testCase.assumeTrue(isunix && ~ismac, 'The mex equals the MATLAB loops bit for bit on Linux only.');
        end

        function [im] = blurredBeads(testCase, sz)
            % uint16 beads blurred by the test PSF plus a little noise, on
            % a background of 100
            n = max(3, round(prod(sz) / 4000));
            pos = (sz - 1) .* rand(n, 3) + 1;
            im = pktest.beads(sz, pos, testCase.BlurSigma, 1000, testCase.Background);
            im = uint16(im + 3 * rand(sz, 'single'));
        end

        function [lambdas] = doubleLambdas(testCase, im, psf, nIter)
            % RL's acceleration factors (before clamping) of the MATLAB loop's
            % iterations in double precision
            lambdas = zeros(1, nIter);
            sz = size(im);
            H = decon_psf2otf(double(psf) ./ sum(double(psf(:))), sz);
            I = max(double(im) - testCase.Background, 0);
            J_2 = I;
            J_3 = zeros(sz);
            J_4 = zeros(numel(I), 1);
            lambda = 0;
            for k = 1 : nIter
                if k > 2
                    lambda = sum((J_2(:) - Y(:)) .* J_4) / (sum(J_4 .* J_4) + eps);
                    lambdas(k) = lambda;
                    lambda = max(min(lambda, 1), 0);
                end
                if k > 1
                    J_4 = J_2(:) - Y(:);
                end
                Y = max(J_2 + lambda * (J_2 - J_3), 0);
                R = I ./ max(real(ifftn(H .* fftn(Y))), eps);
                J_3 = J_2;
                J_2 = max(Y .* real(ifftn(conj(H) .* fftn(R))), 0);
            end
        end

        function [out, errMat] = runDecon(testCase, method, im, psf, varargin)
            % NUMIT can be given as a name-value pair (default 10 for RL, 2 for OMW)
            nIter = 10;
            if strcmp(method, 'omw')
                nIter = 2;
            end
            ind = find(strcmp(varargin(1 : 2 : end), 'NUMIT'));
            if ~isempty(ind)
                nIter = varargin{2 * ind};
                varargin(2 * ind - 1 : 2 * ind) = [];
            end
            switch method
                case 'simplified'
                    [out, errMat] = decon_lucy_function(im, psf, nIter, 'useGPU', false, varargin{:});
                case 'omw'
                    b = testCase.cropToSize(testCase.BackProjector, size(im, 1 : 3));
                    if ismatrix(psf)
                        b = b(:, :, ceil(end / 2));
                    end
                    [out, errMat] = decon_lucy_omw_function(im, psf, b, nIter, 'useGPU', false, varargin{:});
            end
        end
    end

    methods (Static, Access = private)
        function [psf] = cropToSize(psf, sz)
            % central crop to at most sz, as RLdecon crops PSFs larger than the data
            pSize = size(psf, 1 : 3);
            s = max(0, floor((pSize - sz) / 2)) + 1;
            t = s + min(pSize, sz) - 1;
            psf = psf(s(1) : t(1), s(2) : t(2), s(3) : t(3));
        end

        function [b] = halve(a, e, grain)
            % piece boundaries of MATLAB's split of the lines a .. e - 1
            if e - a > grain
                m = a + floor((e - a) / 2);
                b = [DeconRFFTTest.halve(a, m, grain), DeconRFFTTest.halve(m, e, grain)];
                b = unique(b);
                return;
            end
            b = [a, e];
        end

        function [Y] = fftByPieces(X, bounds)
            % rows 1 .. n / 2 + 1 of the FFT of the columns, each piece of
            % columns in groups of four from its start, the rest alone, by
            % MATLAB's fft of a group (or a column) at a time
            n = size(X, 1);
            Y = complex(zeros(n / 2 + 1, size(X, 2), 'single'));
            for p = 1 : numel(bounds) - 1
                a = bounds(p);
                e = bounds(p + 1);
                c = a;
                while c < e
                    m = 4;
                    if e - c < 4
                        m = 1;
                    end
                    F = fft(X(:, c + 1 : c + m), [], 1);
                    Y(:, c + 1 : c + m) = F(1 : n / 2 + 1, :);
                    c = c + m;
                end
            end
        end
    end
end
