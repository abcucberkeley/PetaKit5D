classdef RegistrationXcorrTest < pktest.TestCase
    % the normalized cross-correlation of the registration in the window of
    % the max shifts (normxcorr3_max_shift -> normxcorr3_fast with a bbox ->
    % normxcorr3_window_mex on Linux) against the computation of the whole
    % output: the mex runs the same FFTW transforms as MATLAB's fftn and ifftn
    % and the cumulative sums of the window sums in double, only for the
    % window, so every value is the same bit for bit (also signed zeros),
    % whether this MATLAB session runs its FFTs in parallel or not. The NCC is
    % that of double precision also for data with a large offset and flat
    % regions.

    properties (TestParameter)
        % template and image sizes: even and odd, a length of 64 along y, x
        % or z of the output (FFTW's transforms of four and of one line of 64
        % differ in the last bits), lines MATLAB splits among threads or
        % doesn't, a template of one row or slice, 2D, equal sizes
        xcorrSize = struct('evenOdd', [20, 17, 9, 31, 33, 15], 'y64', [32, 20, 5, 33, 25, 12], ...
            'x64', [10, 30, 7, 20, 35, 9], 'z64', [7, 9, 31, 9, 13, 34], 'split', [40, 64, 6, 60, 101, 20], ...
            'thinTemplate', [1, 5, 3, 2, 7, 4], 'singleSlice', [30, 20, 1, 50, 60, 3], 'twoD', [30, 40, 1, 50, 60, 1], ...
            'equal', [16, 16, 16, 16, 16, 16])
    end

    methods (Test, TestTags = {'Unit'})
        function windowEqualsTheWholeOutput(testCase, xcorrSize)
            % the window of normxcorr3_fast (whole output, corners, single
            % elements, random parts) equals the part of the whole output
            [T, A] = testCase.pair(xcorrSize(1 : 3), xcorrSize(4 : 6));
            C = normxcorr3_fast(T, A);
            sz = size(T, 1 : 3) + size(A, 1 : 3) - 1;
            for bbox = testCase.windows(sz)
                Cw = normxcorr3_fast(T, A, 'full', bbox');
                testCase.verifyBitEqual(Cw, C(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)), mat2str(bbox'));
            end
        end

        function mexPartsEqualMatlabs(testCase, xcorrSize)
            % the mex's correlation and window sums equal those of the MATLAB
            % code (fftn/ifftn, cumulative sums in double) in the window
            testCase.assumeExactMex();
            [T, A] = testCase.pair(xcorrSize(1 : 3), xcorrSize(4 : 6));
            szT = size(T, 1 : 3);
            sz = szT + size(A, 1 : 3) - 1;
            corr = real(ifftn(fftn(A, sz) .* fftn(T(end : -1 : 1, end : -1 : 1, end : -1 : 1), sz)));
            grain = testCase.matlabGrain();
            for bbox = testCase.windows(sz)
                b = bbox';
                [sA, sA2, sT] = testCase.windowSums(T, A, b);
                [c, s1, s2, s3] = normxcorr3_window_mex('ncc', T, A, b, 4, grain);
                testCase.verifyBitEqual(c, corr(b(1) : b(4), b(2) : b(5), b(3) : b(6)), ['corr ', mat2str(b)]);
                testCase.verifyBitEqual(s1, sA, ['sA ', mat2str(b)]);
                testCase.verifyBitEqual(s2, sA2, ['sA2 ', mat2str(b)]);
                testCase.verifyBitEqual(s3, sT, ['sT ', mat2str(b)]);
            end
        end

        function maxShiftEqualsTheWholeOutputsMax(testCase)
            % normxcorr3_max_shift finds the shift of a moved texture and gives
            % the same shift, maximum and window as the whole output did
            im = pktest.texture([90, 80, 40]);
            A = single(im(11 : 70, 6 : 65, 6 : 35));
            T = single(im(25 : 54, 18 : 47, 13 : 27));
            maxShifts = [-20, -15, -8; 25, 20, 9];
            [off, mx, C] = normxcorr3_max_shift(T, A, maxShifts);
            [offRef, mxRef, CRef] = testCase.maxShiftOfWholeOutput(T, A, maxShifts);
            testCase.verifyEqual(off, offRef);
            testCase.verifyBitEqual(mx, mxRef, 'max');
            testCase.verifyBitEqual(C, CRef, 'window');
            % T starts at (15, 13, 8) in A
            testCase.verifyEqual(off, [14, 12, 7]);
            % a symmetric max shift (one row)
            [off, mx] = normxcorr3_max_shift(T, A, [20, 15, 8]);
            [offRef, mxRef] = testCase.maxShiftOfWholeOutput(T, A, [20, 15, 8]);
            testCase.verifyEqual(off, offRef);
            testCase.verifyBitEqual(mx, mxRef, 'max');
        end

        function fastFFTIsCloseInEveryWindow(testCase, xcorrSize)
            % fastFFT's smaller FFTs (no wrap-around in the window) change the
            % correlation only by rounding; the integral images stay the same
            testCase.assumeExactMex();
            [T, A] = testCase.pair(xcorrSize(1 : 3), xcorrSize(4 : 6));
            szT = size(T, 1 : 3);
            sz = szT + size(A, 1 : 3) - 1;
            corr = real(ifftn(fftn(A, sz) .* fftn(T(end : -1 : 1, end : -1 : 1, end : -1 : 1), sz)));
            tol = 1e-5 * max(abs(corr(:)));
            for bbox = testCase.windows(sz)
                b = bbox';
                sA = testCase.windowSums(T, A, b);
                [c, s1] = normxcorr3_window_mex('ncc', T, A, b, 4, 32, true);
                ref = corr(b(1) : b(4), b(2) : b(5), b(3) : b(6));
                testCase.verifyEqual(size(c), size(ref), ['corr ', mat2str(b)]);
                testCase.verifyLessThanOrEqual(max(abs(c(:) - ref(:))), tol, ['corr ', mat2str(b)]);
                testCase.verifyBitEqual(s1, sA, ['sA ', mat2str(b)]);
            end
        end

        function fastFFTFindsTheSameShift(testCase)
            im = pktest.texture([90, 80, 40]);
            A = single(im(11 : 70, 6 : 65, 6 : 35));
            T = single(im(25 : 54, 18 : 47, 13 : 27));
            maxShifts = [-20, -15, -8; 25, 20, 9];
            [off, mx] = normxcorr3_max_shift(T, A, maxShifts);
            [offFast, mxFast] = normxcorr3_max_shift(T, A, maxShifts, true);
            testCase.verifyEqual(offFast, off);
            testCase.verifyEqual(mxFast, mx, 'AbsTol', single(1e-5));
        end

        function lowContrastOnALargeOffset(testCase)
            % data with a large offset (16-bit data shifted by 2^15) and a
            % contrast of a few counts: in single precision the correlation
            % and the window sums of the NCC kept the offset, whose rounding
            % was larger than the NCC itself (|NCC| up to about 10, the peak
            % at a wrong shift). The NCC is now that of double precision.
            im = round(pktest.texture([90, 80, 40], 2, [2 ^ 15 - 6, 2 ^ 15 + 6]));
            A = single(im(11 : 70, 6 : 65, 6 : 35));
            T = single(im(25 : 54, 18 : 47, 13 : 27));
            maxShifts = [-20, -15, -8; 25, 20, 9];
            [off, mx, C] = normxcorr3_max_shift(T, A, maxShifts);
            testCase.verifyEqual(off, [14, 12, 7]);
            testCase.verifyEqual(mx, single(1), 'AbsTol', single(1e-5));
            testCase.verifyLessThanOrEqual(max(abs(C(:))), 1 + 1e-5);
            [s, t] = testCase.maxShiftWindow(T, A, maxShifts);
            Cref = testCase.directNCC(T, A, s, t);
            testCase.verifyLessThanOrEqual(max(abs(double(C(:)) - Cref(:))), 1e-4);
        end

        function flatWindowsAreZero(testCase)
            % windows of an exactly flat region (as the background of
            % flat-field corrected data) have no NCC: 0, not the rounding
            % errors of their sums divided by each other (here, with a wide
            % range of intensities, up to about 4.5). The others are the NCC
            % of double precision.
            im = pktest.texture([90, 80, 40]);
            im = single(exp(12 * (im - 100) / 1900));
            im(:, 31 : end, :) = single(1000.1);
            A = im(11 : 70, 6 : 65, 6 : 35);
            T = im(25 : 54, 18 : 47, 13 : 27);
            maxShifts = [-20, -15, -8; 25, 30, 9];
            [off, mx, C] = normxcorr3_max_shift(T, A, maxShifts);
            testCase.verifyEqual(off, [14, 12, 7]);
            testCase.verifyEqual(mx, single(1), 'AbsTol', single(1e-5));
            testCase.verifyLessThanOrEqual(max(abs(C(:))), 1 + 1e-5);
            [s, t] = testCase.maxShiftWindow(T, A, maxShifts);
            [Cref, flat] = testCase.directNCC(T, A, s, t);
            testCase.assertGreaterThan(nnz(flat), 100);
            testCase.verifyEqual(nnz(C(flat)), 0, 'The NCC of flat windows is not 0.');
            testCase.verifyLessThanOrEqual(max(abs(double(C(:)) - Cref(:))), 1e-4);
            % also without the mex
            testCase.applyFixture(pktest.FunctionStubFixture('normxcorr3_window_mex'));
            testCase.expectMexFallback('normxcorr3_window_mex');
            [~, ~, Cm] = normxcorr3_max_shift(T, A, maxShifts);
            testCase.verifyEqual(Cm, C, 'AbsTol', single(1e-5));
            testCase.verifyEqual(nnz(Cm(flat)), 0);
        end

        function faintRegionsHaveTheirNCC(testCase)
            % windows of a faint region (a contrast of a few counts next to
            % hundreds) are not flat to single precision: their NCC is that of
            % double precision, not 0
            im = pktest.texture([90, 80, 40]);
            f = pktest.texture([90, 80, 40], 2, [0, 6]);
            im(:, 31 : end, :) = 1000 + f(:, 31 : end, :);
            A = im(11 : 70, 6 : 65, 6 : 35);
            T = im(25 : 54, 18 : 47, 13 : 27);
            maxShifts = [-20, -15, -8; 25, 30, 9];
            [off, ~, C] = normxcorr3_max_shift(T, A, maxShifts);
            testCase.verifyEqual(off, [14, 12, 7]);
            [s, t] = testCase.maxShiftWindow(T, A, maxShifts);
            Cref = testCase.directNCC(T, A, s, t);
            testCase.verifyLessThanOrEqual(max(abs(double(C(:)) - Cref(:))), 1e-4);
        end

        function validAndSameShapes(testCase)
            % 'valid' started along x at the template's size along y, which
            % gave a wrong part (and size) for templates of different sizes
            % along y and x
            [T, A] = testCase.pair([6, 9, 4], [20, 22, 11]);
            szT = size(T, 1 : 3);
            C = normxcorr3_fast(T, A);
            Cv = normxcorr3_fast(T, A, 'valid');
            testCase.verifyEqual(size(Cv), size(A, 1 : 3) - szT + 1);
            testCase.verifyEqual(Cv, C(szT(1) : end - szT(1) + 1, szT(2) : end - szT(2) + 1, szT(3) : end - szT(3) + 1));
            p = fix((szT - 1) / 2);
            testCase.verifyEqual(normxcorr3_fast(T, A, 'same'), C(p(1) + 1 : p(1) + 20, p(2) + 1 : p(2) + 22, p(3) + 1 : p(3) + 11));
        end

        function maxShiftOf2DImages(testCase)
            % 2D images errored: the window's upper bound was capped with
            % size(C), which has 2 elements for a 2D output
            im = pktest.texture([120, 100, 1]);
            A = single(im(11 : 90, 6 : 85));
            T = single(im(31 : 70, 21 : 60));
            [off, mx] = normxcorr3_max_shift(T, A, [-30, -30, -1; 30, 30, 1]);
            testCase.verifyEqual(off, [20, 15, 0]);
            testCase.verifyEqual(mx, single(1), 'AbsTol', single(1e-3));
        end

        function resultsDoNotDependOnTheThreadCount(testCase)
            % as MATLAB's FFTs: the groups of lines don't depend on the threads
            testCase.assumeExactMex();
            [T, A] = testCase.pair([40, 64, 6], [60, 101, 20]);
            bbox = [10, 20, 3, 80, 150, 22];
            for grain = [32, 0]
                out = cell(1, 2);
                nThreads = [1, 4];
                for i = 1 : 2
                    r = cell(1, 3);
                    [r{:}] = normxcorr3_window_mex('ncc', T, A, bbox, nThreads(i), grain);
                    out{i} = r;
                end
                testCase.verifyEqual(out{2}, out{1}, sprintf('grain %d', grain));
            end
        end

        function nccDoesNotDependOnTheThreadCount(testCase)
            % the means, the sum and the variances of the template and the
            % image (89000 or more elements each) were MATLAB's, which change
            % in the last bits with the number of threads, and with them the
            % NCC of double data in most values (of single data rarely: at
            % windows at the flatness threshold). Now the NCC is the same at
            % 2, 4 and 24 threads as at 1: the MATLAB code (double) and the
            % mex (single, on Linux).
            nt = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(nt));
            T = rand([50, 45, 41]) * 100;
            A = rand([70, 66, 50]) * 100;
            bbox = [30, 25, 20, 90, 85, 70];
            threads = [1, 2, 4, 24];
            for cls = {'double', 'single'}
                C = cell(1, numel(threads));
                for k = 1 : numel(threads)
                    maxNumCompThreads(threads(k));
                    C{k} = normxcorr3_fast(cast(T, cls{1}), cast(A, cls{1}), 'full', bbox);
                end
                for k = 2 : numel(threads)
                    testCase.verifyBitEqual(C{k}, C{1}, sprintf('%s data, %d threads', cls{1}, threads(k)));
                end
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
            testCase.assumeExactMex();
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
                testCase.verifyEqual(normxcorr3_window_mex('fft1', X, 4, 32), refSplit, sprintf('%d columns, split', N));
                testCase.verifyEqual(normxcorr3_window_mex('fft1', X, 4, 0), refSerial, sprintf('%d columns, no split', N));
                % this session's MATLAB FFT equals one of them
                Y = fft(X, [], 1);
                Y = Y(1 : 33, :);
                testCase.verifyTrue(isequal(Y, refSplit) || isequal(Y, refSerial));
            end
        end

        function withoutTheMexTheMatlabCodeRuns(testCase)
            % where the mex is missing (or not used), the whole output is
            % computed and cropped
            [T, A] = testCase.pair([20, 17, 9], [31, 33, 15]);
            C = normxcorr3_fast(T, A);
            bbox = [5, 8, 3, 40, 37, 20];
            testCase.applyFixture(pktest.FunctionStubFixture('normxcorr3_window_mex'));
            testCase.expectMexFallback('normxcorr3_window_mex');
            Cw = normxcorr3_fast(T, A, 'full', bbox);
            testCase.verifyBitEqual(Cw, C(5 : 40, 8 : 37, 3 : 20), 'window');
            % fastFFT then gives the same as without it
            Cw = normxcorr3_fast(T, A, 'full', bbox, true);
            testCase.verifyBitEqual(Cw, C(5 : 40, 8 : 37, 3 : 20), 'window, fastFFT');
        end

        function unsupportedCallsRunTheMatlabCode(testCase)
            % double data runs the MATLAB code (by design, without a fallback
            % record of the mex)
            [T, A] = testCase.pair([12, 10, 5], [20, 22, 9]);
            T = double(T);
            A = double(A);
            C = normxcorr3_fast(T, A);
            mex_fallback('-reset');
            Cw = normxcorr3_fast(T, A, 'full', [2, 3, 1, 30, 31, 13]);
            testCase.verifyClass(Cw, 'double');
            testCase.verifyEqual(Cw, C(2 : 30, 3 : 31, 1 : 13));
            records = mex_fallback('-log');
            testCase.verifyFalse(any(strcmp({records.mex}, 'normxcorr3_window_mex')));
            % SHAPE must be 'full' with a bbox
            testCase.verifyError(@() normxcorr3_fast(T, A, 'same', [1, 1, 1, 2, 2, 2]), ?MException);
        end

        function mexRejectsInvalidInputs(testCase)
            testCase.assumeMexAvailable('normxcorr3_window_mex');
            [T, A] = testCase.pair([5, 6, 4], [9, 11, 7]);
            id = 'normxcorr3_window_mex:inputError';
            bbox = [1, 1, 1, 13, 16, 10];
            import matlab.unittest.constraints.Throws
            testCase.verifyThat(@() normxcorr3_window_mex('ncc', double(T), A, bbox, 1, 32), Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('ncc', A, T, bbox, 1, 32), Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('ncc', T, A, [0, 1, 1, 13, 16, 10], 1, 32), Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('ncc', T, A, [1, 1, 1, 14, 16, 10], 1, 32), Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('ncc', T, A, [3, 1, 1, 2, 16, 10], 1, 32), Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('ncc', T, A, [1.5, 1, 1, 13, 16, 10], 1, 32), Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('ncc', T, A, bbox(1 : 5), 1, 32), Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('ncc', T, cat(4, A, A), bbox, 1, 32), Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('ncc', T, A, bbox, 1), Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('ncc', T(1, 1, :), A(1 : 1, 1, :), [1, 1, 1, 1, 1, 10], 1, 32), ...
                Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('fft1', T, 1, 32), Throws(id, 'WhenNargoutIs', 1));
            testCase.verifyThat(@() normxcorr3_window_mex('nccc', T, A, bbox, 1, 32), Throws(id, 'WhenNargoutIs', 1));
        end
    end

    methods (Access = private)
        function [] = assumeExactMex(testCase)
            % the mex gives the same bits as MATLAB's FFTs on Linux, where it
            % is used (as decon_lucy_rfft_mex; on macOS they differ in the last
            % bits)
            testCase.assumeMexAvailable('normxcorr3_window_mex');
            testCase.assumeTrue(isunix && ~ismac, 'The mex equals MATLAB''s FFTs bit for bit on Linux only.');
        end

        function [] = verifyBitEqual(testCase, actual, expected, name)
            testCase.verifyEqual(size(actual), size(expected), [name, ': size']);
            testCase.verifyClass(actual, class(expected));
            if isequal(size(actual), size(expected)) && isa(actual, class(expected))
                if isa(expected, 'double')
                    nDiff = nnz(typecast(actual(:), 'uint64') ~= typecast(expected(:), 'uint64'));
                else
                    nDiff = nnz(typecast(actual(:), 'uint32') ~= typecast(expected(:), 'uint32'));
                end
                testCase.verifyEqual(nDiff, 0, sprintf('%s: %d of %d values differ in their bits.', name, nDiff, numel(expected)));
            end
        end
    end

    methods (Static, Access = private)
        function [T, A] = pair(szT, szA)
            % image and template with a little structure and zero borders
            % (zero lines in the padding and in the data), single
            A = single(rand(szA) * 100);
            A(1 : min(2, end), :, :) = 0;
            T = single(rand(szT) * 50);
            T(:, end, :) = 0;
            if szT(3) > 2
                T(:, :, 1) = 0;
            end
        end

        function [bboxes] = windows(sz)
            % whole output, first and last element, edges and random parts
            % (columns)
            bboxes = [[ones(1, 3), sz]; [1, 1, 1, 1, 1, 1]; [sz, sz]; [1, sz(2), 1, sz(1), sz(2), sz(3)]];
            s = max(1, floor(sz / 3));
            bboxes = [bboxes; [s, min(sz, s + floor(sz / 4))]];
            for k = 1 : 3
                a = arrayfun(@(n) randi(n), sz);
                b = arrayfun(@(lo, n) randi([lo, n]), a, sz);
                bboxes = [bboxes; [a, b]];
            end
            bboxes = bboxes';
        end

        function [grain] = matlabGrain()
            % the split of this session's FFTs, as normxcorr3_fast finds it
            grain = 32;
            X = rand(64, 2000, 'single');
            Y = fft(X, [], 1);
            if isequal(normxcorr3_window_mex('fft1', X, 4, 0), Y(1 : 33, :)) ...
                    && ~isequal(normxcorr3_window_mex('fft1', X, 4, 32), Y(1 : 33, :))
                grain = 0;
            end
        end

        function [sA, sA2, sT] = windowSums(T, A, bbox)
            % the window sums of normxcorr3_fast's MATLAB code in the window
            % bbox: of A, A .* A (zero outside A) and of T where it overlaps A,
            % in double, by cumulative sums along y, x and z
            szT = size(T, 1 : 3);
            szA = size(A, 1 : 3);
            [loT, hiT, loA, hiA] = deal(cell(1, 3));
            for d = 1 : 3
                k = bbox(d) : bbox(d + 3);
                loT{d} = max(1, szT(d) - k + 1);
                hiT{d} = min(szT(d), szT(d) - k + szA(d));
                loA{d} = max(1, k - szT(d) + 1);
                hiA{d} = min(k, szA(d));
            end
            X = double(A);
            sA = RegistrationXcorrTest.boxSums(X, loA, hiA);
            sA2 = RegistrationXcorrTest.boxSums(X .* X, loA, hiA);
            sT = RegistrationXcorrTest.boxSums(double(T), loT, hiT);
        end

        function [S] = boxSums(S, lo, hi)
            % sums over the boxes lo{d}(k) .. hi{d}(k) along each dimension
            for d = 1 : 3
                P = cumsum(S, d);
                idx = repmat({':'}, 1, 3);
                idx{d} = hi{d};
                W = P(idx{:});
                m = lo{d} > 1;
                if any(m)
                    idx{d} = lo{d}(m) - 1;
                    L = P(idx{:});
                    idx{d} = find(m);
                    W(idx{:}) = W(idx{:}) - L;
                end
                S = W;
            end
        end

        function [s, t] = maxShiftWindow(T, A, maxShifts)
            % the window of the full output that normxcorr3_max_shift takes
            sz_t = size(T, 1 : 3);
            s = max(1, ceil(maxShifts(1, :) + sz_t));
            t = min(floor(maxShifts(2, :) + sz_t), sz_t + size(A, 1 : 3) - 1);
        end

        function [C, flat] = directNCC(T, A, s, t)
            % the NCC of the windows s .. t of the full output, each from its
            % elements in double precision (the image zero-padded); flat:
            % windows of one value (0 there)
            T = double(T);
            szT = size(T, 1 : 3);
            szA = size(A, 1 : 3);
            Ap = zeros(szA + 2 * szT);
            Ap(szT(1) + 1 : szT(1) + szA(1), szT(2) + 1 : szT(2) + szA(2), szT(3) + 1 : szT(3) + szA(3)) = double(A);
            T0 = T(:) - mean(T(:));
            sT = std(T(:));
            n = numel(T);
            C = zeros(t - s + 1);
            flat = false(t - s + 1);
            for z = s(3) : t(3)
                for x = s(2) : t(2)
                    for y = s(1) : t(1)
                        w = Ap(y + 1 : y + szT(1), x + 1 : x + szT(2), z + 1 : z + szT(3));
                        sw = std(w(:));
                        k = {y - s(1) + 1, x - s(2) + 1, z - s(3) + 1};
                        flat(k{:}) = sw == 0;
                        if sw > 0
                            C(k{:}) = sum(T0 .* w(:)) / ((n - 1) * sT * sw);
                        end
                    end
                end
            end
        end

        function [off, mx, C] = maxShiftOfWholeOutput(T, A, maxShifts)
            % normxcorr3_max_shift of 2a851bf: the whole output, then cropped
            if size(maxShifts, 1) == 1
                maxShifts = [-maxShifts; maxShifts];
            end
            C = normxcorr3_fast(T, A);
            sz_t = size(T, 1 : 3);
            s = max(1, ceil(maxShifts(1, :) + sz_t));
            t = min(floor(maxShifts(2, :) + sz_t), size(C));
            C = C(s(1) : t(1), s(2) : t(2), s(3) : t(3));
            [mx, ind] = max(C(:));
            [y, x, z] = ind2sub(size(C), ind);
            off = [y, x, z] + (s - 1) - sz_t;
        end

        function [b] = halve(a, e, grain)
            % piece boundaries of MATLAB's split of the lines a .. e - 1
            if e - a > grain
                m = a + floor((e - a) / 2);
                b = [RegistrationXcorrTest.halve(a, m, grain), RegistrationXcorrTest.halve(m, e, grain)];
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
