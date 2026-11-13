classdef DeconAccelTest < pktest.TestCase
    % the acceleration step of the CPU 'simplified' RL (decon_lucy_accel:
    % decon_lucy_accel_mex from 89000 elements on, else the MATLAB code) and
    % the sums in double in an order that doesn't depend on the number of
    % threads (decon_sum_double: MATLAB's order at one thread), against the
    % MATLAB expressions, bit for bit (also -0 and NaN), and decon results at
    % different thread counts. MATLAB splits sums of 89000 or more elements
    % among its threads, so lambda, and with it every voxel of the result,
    % changed in the last bits with a job's CPU count.

    properties (Constant)
        Background = 100
        PSFSigma = [1.5, 2]
        % above MATLAB's 89000 elements, with a last block of 4096 that isn't
        % full and a last group of four lanes that isn't either
        Sz = [47, 45, 49]
    end

    properties (TestParameter)
        dataClass = {'single', 'double'}
        % around MATLAB's switches (44 elements: lanes, 89000: blocks)
        sumSize = {1, 43, 44, 4097, 88999, 89000, 89001, 4096 * 30, 4096 * 30 + 3, 250047}
        lambdaCase = {'aboveOne', 'negative', 'nan', 'negativeZero', 'zeroDenominator', 'inf', 'inside'}
    end

    methods (Test, TestTags = {'Unit'})
        function sumIsMatlabsAtOneThread(testCase, dataClass, sumSize)
            % data with a wide range of exponents, so that another order
            % gives other bits. Below 89000 elements decon_sum_double is
            % MATLAB's sum; above, MATLAB's order at one thread (as found for
            % R2024b and R2025b on Linux)
            x = DeconAccelTest.wideRange(sumSize, dataClass);
            nt = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(nt));
            maxNumCompThreads(1);
            ref = sum(x, 'double');
            for t = unique([1, 2, min(4, feature('numcores'))])
                maxNumCompThreads(t);
                s = decon_sum_double(x);
                if sumSize < 89000
                    testCase.verifyEqual(DeconAccelTest.bits(s), DeconAccelTest.bits(sum(x, 'double')), sprintf('%d threads', t));
                end
                if isunix && ~ismac
                    testCase.verifyEqual(DeconAccelTest.bits(s), DeconAccelTest.bits(ref), ...
                        sprintf('%d threads: not MATLAB''s sum at one thread', t));
                end
            end
            % 3D arrays are summed as x(:)
            testCase.verifyEqual(DeconAccelTest.bits(decon_sum_double(reshape(x, 1, 1, []))), DeconAccelTest.bits(decon_sum_double(x)));
        end

        function sumDoesNotDependOnTheThreadCount(testCase, dataClass)
            % MATLAB's own sum of these elements changes with the thread count
            testCase.assumeGreaterThanOrEqual(feature('numcores'), 2, 'Needs at least 2 cores.');
            nt = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(nt));
            x = DeconAccelTest.wideRange(250047, dataClass);
            threads = unique([1, 2, min(4, feature('numcores'))]);
            s = zeros(size(threads));
            for k = 1 : numel(threads)
                maxNumCompThreads(threads(k));
                s(k) = decon_sum_double(x);
            end
            testCase.verifyEqual(DeconAccelTest.bits(s), repmat(DeconAccelTest.bits(s(1)), numel(s), 1));
        end

        function lambdaDoesNotDependOnTheThreadCount(testCase, dataClass)
            % products over a wide range of exponents, whose sums MATLAB's
            % sum gives with other last bits at each thread count (so did
            % the acceleration factor of 70be771's decon_lucy_function)
            testCase.assumeGreaterThanOrEqual(feature('numcores'), 2, 'Needs at least 2 cores.');
            nt = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(nt));
            n = 250047;
            J_2 = cast(rand(n, 1) * 1000, dataClass);
            J_4 = DeconAccelTest.wideRange(n, dataClass);
            Y = J_2 - 0.5 * J_4 - 0.1 * DeconAccelTest.wideRange(n, dataClass);
            threads = unique([1, 2, min(4, feature('numcores'))]);
            lambdas = zeros(size(threads));
            for k = 1 : numel(threads)
                maxNumCompThreads(threads(k));
                lambdas(k) = decon_lucy_accel('lambda', J_2, Y, J_4);
            end
            testCase.verifyEqual(DeconAccelTest.bits(lambdas), repmat(DeconAccelTest.bits(lambdas(1)), numel(lambdas), 1));
            testCase.verifyGreaterThan(lambdas(1), 0);
            testCase.verifyLessThan(lambdas(1), 1);
        end

        function mexSumEqualsTheMatlabCode(testCase, dataClass)
            % special values: -0 everywhere (the sum is -0), a +0 among -0s
            % (+0), NaN, Inf, -Inf, subnormal numbers
            testCase.assumeMexAvailable('decon_lucy_accel_mex');
            xs = DeconAccelTest.specialVectors(prod(testCase.Sz), dataClass);
            ref = zeros(1, numel(xs));
            for i = 1 : numel(xs)
                ref(i) = decon_sum_double(xs{i});
            end
            testCase.applyFixture(pktest.FunctionStubFixture('decon_lucy_accel_mex'));
            testCase.expectMexFallback('decon_lucy_accel_mex');
            for i = 1 : numel(xs)
                s = decon_sum_double(xs{i});
                testCase.verifyTrue(isequal(DeconAccelTest.bits(s), DeconAccelTest.bits(ref(i))) || (isnan(s) && isnan(ref(i))), ...
                    sprintf('Vector %d: mex %.17g, MATLAB code %.17g.', i, ref(i), s));
            end
            testCase.verifyEqual(DeconAccelTest.bits(ref(1 : 3)), DeconAccelTest.bits([-0, 0, 0]));
        end

        function lambdaStepIsTheMatlabExpressions(testCase, dataClass)
            % the products in the class of the data (single), summed in
            % double, then the clamp; the new J_4 is J_2(:) - Y(:), a column
            % (lambda of finite inputs, also over a wide range of exponents,
            % whose sums depend on the order, and with large products in the
            % last three elements, three lanes of the last group of four; NaN
            % and Inf in J_4's difference)
            testCase.assumeMexAvailable('decon_lucy_accel_mex');
            n = prod(testCase.Sz);
            for c = 1 : 4
                nonFinite = c == 3;
                [J_2, Y, J_4] = DeconAccelTest.accelInputs(n, dataClass, nonFinite);
                if c == 2
                    J_4 = DeconAccelTest.wideRange(n, dataClass);
                    Y = J_2 - 0.5 * J_4 - 0.1 * DeconAccelTest.wideRange(n, dataClass);
                elseif c == 4
                    J_4(n - 2 : n) = cast([3e9; -5e9; 7e9], dataClass);
                    Y(n - 2 : n) = J_2(n - 2 : n) - 0.5 * J_4(n - 2 : n);
                end
                J_2 = reshape(J_2, testCase.Sz);
                [lambda, J4n] = decon_lucy_accel('lambda', J_2, Y, J_4);
                ref = max(min(decon_sum_double((J_2(:) - Y(:)) .* J_4) / (decon_sum_double(J_4 .* J_4) + eps), 1), 0);
                testCase.verifyEqual(DeconAccelTest.bits(lambda), DeconAccelTest.bits(ref));
                testCase.verifyEqual(size(J4n), [numel(J_2), 1]);
                testCase.verifyClass(J4n, dataClass);
                testCase.verifyEqual(DeconAccelTest.bits(J4n), DeconAccelTest.bits(J_2(:) - Y(:)));
                if ~nonFinite
                    testCase.verifyGreaterThan(lambda, 0);
                    testCase.verifyLessThan(lambda, 1);
                    % and the sums of the mex itself
                    [~, ~, s] = decon_lucy_accel_mex('lambda', J_2, Y, J_4, 3);
                    testCase.verifyEqual(DeconAccelTest.bits(s), DeconAccelTest.bits([decon_sum_double((J_2(:) - Y(:)) .* J_4), ...
                        decon_sum_double(J_4 .* J_4)]));
                end
            end
        end

        function predictIsTheMatlabExpression(testCase, dataClass)
            % max(J_2 + lambda * (J_2 - J_3), 0): lambda * x in double, then
            % rounded to single; max is 0 for NaN and keeps -0
            testCase.assumeMexAvailable('decon_lucy_accel_mex');
            [J_2, ~, ~, J_3] = DeconAccelTest.accelInputs(prod(testCase.Sz), dataClass, true);
            J_2 = reshape(J_2, testCase.Sz);
            J_3 = reshape(J_3, testCase.Sz);
            for lambda = [0, -0, 1, 0.5, 0.123456789012345678, 1 - eps, 2 ^ -60]
                Y = decon_lucy_accel('predict', J_2, J_3, lambda);
                ref = max(J_2 + lambda * (J_2 - J_3), 0);
                testCase.verifyClass(Y, dataClass);
                testCase.verifySize(Y, size(J_2));
                testCase.verifyEqual(DeconAccelTest.bits(Y), DeconAccelTest.bits(ref), sprintf('lambda %.17g', lambda));
            end
            % the first iteration: J_3 and lambda zero
            Y = decon_lucy_accel('predict', J_2, zeros(size(J_2), dataClass), 0);
            testCase.verifyEqual(DeconAccelTest.bits(Y), DeconAccelTest.bits(max(J_2 + 0 * (J_2 - zeros(size(J_2), dataClass)), 0)));
        end

        function lambdaIsClampedAsMatlabDoes(testCase, lambdaCase)
            % max(min(lambda, 1), 0): MATLAB's min gives 1 for NaN, and min
            % and max keep -0
            testCase.assumeMexAvailable('decon_lucy_accel_mex');
            n = prod(testCase.Sz);
            J_2 = rand(n, 1, 'single') * 100;
            J_4 = rand(n, 1, 'single') - 0.5;
            switch lambdaCase
                case 'aboveOne'
                    J_4 = ones(n, 1, 'single');
                    Y = J_2 - 2;
                    expected = 1;
                case 'negative'
                    J_4 = ones(n, 1, 'single');
                    Y = J_2 + 1;
                    expected = 0;
                case 'nan'
                    % 0 * Inf
                    Y = J_2;
                    J_4(7) = Inf;
                    expected = 1;
                case 'negativeZero'
                    % all products are 0 .* -1 = -0
                    Y = J_2;
                    J_4 = -ones(n, 1, 'single');
                    expected = -0;
                case 'zeroDenominator'
                    Y = J_2 - 1;
                    J_4 = zeros(n, 1, 'single');
                    expected = 0;
                case 'inf'
                    J_4 = ones(n, 1, 'single');
                    J_2(5) = Inf;
                    Y = J_2 - 1;
                    Y(5) = 0;
                    expected = 1;
                case 'inside'
                    Y = J_2 - 0.5 * J_4;
                    expected = [];
            end
            [lambda, J4n] = decon_lucy_accel('lambda', J_2, Y, J_4);
            ref = max(min(decon_sum_double((J_2(:) - Y(:)) .* J_4) / (decon_sum_double(J_4 .* J_4) + eps), 1), 0);
            testCase.verifyEqual(DeconAccelTest.bits(lambda), DeconAccelTest.bits(ref));
            if ~isempty(expected)
                testCase.verifyEqual(DeconAccelTest.bits(lambda), DeconAccelTest.bits(expected));
            else
                testCase.verifyGreaterThan(lambda, 0.4);
                testCase.verifyLessThan(lambda, 0.6);
            end
            testCase.verifyEqual(DeconAccelTest.bits(J4n), DeconAccelTest.bits(J_2(:) - Y(:)));
            % the MATLAB code without the mex
            testCase.applyFixture(pktest.FunctionStubFixture('decon_lucy_accel_mex'));
            testCase.expectMexFallback('decon_lucy_accel_mex');
            [lambda, J4n] = decon_lucy_accel('lambda', J_2, Y, J_4);
            testCase.verifyEqual(DeconAccelTest.bits(lambda), DeconAccelTest.bits(ref), 'MATLAB code');
            testCase.verifyEqual(DeconAccelTest.bits(J4n), DeconAccelTest.bits(J_2(:) - Y(:)), 'MATLAB code');
        end

        function loopIsTheMatlabLoop(testCase)
            % decon_lucy_function's iterations (k == 1 and 2: lambda 0, J_4
            % from the second) against the MATLAB loop with the same sums, bit
            % for bit; with the MATLAB FFT steps, whose max(., 0) keeps the -0s
            % of products with zeros of the data
            testCase.assumeMexAvailable('decon_lucy_accel_mex');
            [im, psf] = testCase.texturedData([47, 45, 43]);
            for nIter = [1, 2, 3, 6]
                for useMex = [false, true]
                    if useMex && ~DeconAccelTest.rfftMexDefault()
                        continue;
                    end
                    out = decon_lucy_function(im, psf, nIter, 'Background', testCase.Background, 'useGPU', false, ...
                        'save16bit', false, 'useMex', useMex);
                    ref = testCase.matlabLoop(im, psf, nIter, useMex);
                    testCase.verifyEqual(DeconAccelTest.bits(out), DeconAccelTest.bits(ref), ...
                        sprintf('%d iterations, FFT mex %d', nIter, useMex));
                end
            end
        end

        function withoutTheMexTheMatlabCodeRuns(testCase)
            % a .m file of the same name first on the path hides the mex, as on
            % a platform without it; the same bits (also -0)
            testCase.assumeMexAvailable('decon_lucy_accel_mex');
            [im, psf] = testCase.texturedData([47, 45, 43]);
            args = {'Background', testCase.Background, 'useGPU', false, 'save16bit', false, 'useMex', false};
            ref = decon_lucy_function(im, psf, 6, args{:});
            clear decon_lucy_function
            testCase.applyFixture(pktest.FunctionStubFixture('decon_lucy_accel_mex'));
            testCase.expectMexFallback('decon_lucy_accel_mex');
            out = decon_lucy_function(im, psf, 6, args{:});
            testCase.verifyEqual(DeconAccelTest.bits(out), DeconAccelTest.bits(ref));
            records = mex_fallback('-log');
            testCase.verifyTrue(any(strcmp({records.mex}, 'decon_lucy_accel_mex')));
        end

        function smallArraysRunTheMatlabCode(testCase, dataClass)
            % below 89000 elements, where MATLAB's sums don't depend on the
            % number of threads, the step is MATLAB's own code (the same bits
            % on every platform)
            % (exponents over a wide range: other orders give other bits)
            [J_2, ~, ~, J_3] = DeconAccelTest.accelInputs(88999, dataClass, false);
            J_4 = DeconAccelTest.wideRange(88999, dataClass);
            Y = J_2 - 0.5 * J_4 - 0.1 * DeconAccelTest.wideRange(88999, dataClass);
            [lambda, J4n] = decon_lucy_accel('lambda', J_2, Y, J_4);
            ref = max(min(sum((J_2(:) - Y(:)) .* J_4, 'double') / (sum(J_4 .* J_4, 'double') + eps), 1), 0);
            testCase.verifyEqual(DeconAccelTest.bits(lambda), DeconAccelTest.bits(ref));
            testCase.verifyEqual(DeconAccelTest.bits(J4n), DeconAccelTest.bits(J_2(:) - Y(:)));
            Y = decon_lucy_accel('predict', J_2, J_3, 0.3);
            testCase.verifyEqual(DeconAccelTest.bits(Y), DeconAccelTest.bits(max(J_2 + 0.3 * (J_2 - J_3), 0)));
        end

        function psfNormalizationsUseTheFixedOrderSum(testCase)
            % (the PSF sums of these sizes give the same normalized PSFs at
            % any thread count, so a spy checks the call sites): RLdecon and
            % decon_lucy_function normalize the PSF, decon_lucy_omw_function
            % the PSF and the back projector, omw_backprojector_generation the
            % PSF and sums the OTF magnitudes with decon_sum_double
            im = testCase.texturedData([48, 48, 40]);
            psf = pktest.deconPSF([31, 31, 31], testCase.PSFSigma);
            [~, b] = evalc('omw_backprojector_generation(psf, 0.005, false)');
            root = testCase.createTempFolder();
            psfFn = [root, 'psf.tif'];
            writetiff(psf, psfFn);
            fid = fopen([root, 'decon_sum_double.m'], 'w');
            fprintf(fid, '%s\n', 'function [s] = decon_sum_double(x)', ...
                '% test spy written by DeconAccelTest: records numel(x), returns MATLAB''s sum', ...
                'persistent n', 'if ischar(x)', '    s = n;', '    n = [];', '    return;', 'end', ...
                'n(end + 1) = numel(x);', 's = sum(x(:), ''double'');', 'end');
            fclose(fid);
            addpath(root, '-begin');
            testCase.addTeardown(@() rmpath(root));
            testCase.addTeardown(@() clear('decon_sum_double'));
            clear decon_lucy_function decon_lucy_omw_function
            np = numel(psf);
            evalc(['RLdecon('''', '''', psfFn, 0.108, 0.2, 0.2, ''rawdata'', im, ''RLMethod'', ''simplified'', ', ...
                '''psfGen'', false, ''useGPU'', false, ''Background'', 100, ''DeconIter'', 1, ''save16bit'', false)']);
            testCase.verifyEqual(nnz(decon_sum_double('-calls') == np), 2, 'RLdecon and decon_lucy_function');
            decon_lucy_omw_function(im, psf, b, 1, 'Background', testCase.Background, 'useGPU', false);
            testCase.verifyEqual(nnz(decon_sum_double('-calls') == np), 2, 'decon_lucy_omw_function');
            evalc('omw_backprojector_generation(psf, 0.005, false)');
            testCase.verifyEqual(nnz(decon_sum_double('-calls') == np), 2, 'omw_backprojector_generation');
        end

        function mexRejectsInvalidInputs(testCase)
            testCase.assumeMexAvailable('decon_lucy_accel_mex');
            a = rand(10, 1, 'single');
            id = 'decon_lucy_accel_mex:inputError';
            testCase.verifyError(@() decon_lucy_accel_mex('lambda', a, double(a), a, 1), id);
            testCase.verifyError(@() decon_lucy_accel_mex('lambda', a, a(1 : 9), a, 1), id);
            testCase.verifyError(@() decon_lucy_accel_mex('lambda', complex(a), a, a, 1), id);
            testCase.verifyError(@() decon_lucy_accel_mex('lambda', uint16(a), uint16(a), uint16(a), 1), id);
            testCase.verifyError(@() decon_lucy_accel_mex('lambda', a, a, a), id);
            testCase.verifyError(@() decon_lucy_accel_mex('predict', a, a, single(0.5), 1), id);
            testCase.verifyError(@() decon_lucy_accel_mex('predict', a, a, [0.5, 0.5], 1), id);
            testCase.verifyError(@() decon_lucy_accel_mex('predict', a, double(a), 0.5, 1), id);
            testCase.verifyError(@() decon_lucy_accel_mex('sum', sparse(double(a)), 1), id);
            d = sparse(double(a));
            testCase.verifyError(@() decon_lucy_accel_mex('lambda', d, d, d, 1), id);
            testCase.verifyError(@() decon_lucy_accel_mex('predict', d, d, 0.5, 1), id);
            testCase.verifyError(@() decon_lucy_accel_mex('sum', a, 'x'), id);
            testCase.verifyError(@() decon_lucy_accel_mex('other', a, 1), id);
            testCase.verifyError(@() decon_lucy_accel_mex(1), id);
            % the empty sum is 0
            testCase.verifyEqual(DeconAccelTest.bits(decon_lucy_accel_mex('sum', zeros(0, 1, 'single'), 2)), DeconAccelTest.bits(0));
        end
    end

    methods (Test, TestTags = {'Integration'})
        function deconDoesNotDependOnTheThreadCount(testCase)
            % with MATLAB's sums the acceleration factor changed in the last
            % bits with the thread count, and with it most voxels of large
            % volumes (a demo tile: 40% of them); at this size a change of a
            % few ulps rarely reaches a voxel, so lambdaDoesNotDependOn...
            % checks the sums; also with the MATLAB FFT steps
            testCase.assumeGreaterThanOrEqual(feature('numcores'), 2, 'Needs at least 2 cores.');
            nt = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(nt));
            [im, psf] = testCase.texturedData([128, 128, 64]);
            threads = [1, min(4, feature('numcores'))];
            outs = cell(2, 2);
            for k = 1 : 2
                maxNumCompThreads(threads(k));
                clear decon_lucy_function
                outs{k, 1} = decon_lucy_function(im, psf, 10, 'Background', testCase.Background, 'useGPU', false, 'save16bit', false);
                outs{k, 2} = decon_lucy_function(im(:, :, 1 : 32), psf, 5, 'Background', testCase.Background, 'useGPU', false, ...
                    'save16bit', false, 'useMex', false);
            end
            testCase.verifyEqual(DeconAccelTest.bits(outs{2, 1}), DeconAccelTest.bits(outs{1, 1}), sprintf('%d and %d threads differ.', threads));
            testCase.verifyEqual(DeconAccelTest.bits(outs{2, 2}), DeconAccelTest.bits(outs{1, 2}), ...
                sprintf('%d and %d threads differ with the MATLAB FFT steps.', threads));
        end

        function debugErrorDoesNotDependOnTheThreadCount(testCase)
            % the error metric of the debug output (the mean squared
            % difference of the result and the data, 89000 or more voxels) was
            % MATLAB's mean, which changes in the last bits with the number of
            % threads (here in 2 to 4 of the 4 iterations); now the same at 2,
            % 4 and 24 threads as at 1: decon_lucy_function, and
            % decon_lucy_omw_function with the MATLAB loop and with the FFT
            % mex in pieces
            nt = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(nt));
            [im, psf] = testCase.texturedData([48, 48, 40]);
            [~, b] = evalc('omw_backprojector_generation(psf, 0.005, false)');
            root = testCase.createTempFolder();
            threads = [1, 2, 4, 24];
            useMex = [false, DeconAccelTest.rfftMexDefault()];
            errs = cell(numel(threads), 3);
            for k = 1 : numel(threads)
                maxNumCompThreads(threads(k));
                clear decon_lucy_function decon_lucy_omw_function
                args = {'Background', testCase.Background, 'useGPU', false, 'save16bit', false, 'debug', true};
                [~, errs{k, 1}] = decon_lucy_function(im, psf, 4, args{:}, 'debug_folder', sprintf('%srl_%d', root, k));
                for m = 1 : 2
                    if m == 2 && ~useMex(m)
                        continue;
                    end
                    [~, errs{k, m + 1}] = decon_lucy_omw_function(im, psf, b, 4, args{:}, 'useMex', useMex(m), ...
                        'debug_folder', sprintf('%somw_%d_%d', root, m, k));
                end
            end
            names = {'decon_lucy_function', 'OMW, MATLAB loop', 'OMW, FFT mex'};
            for i = 1 : 3
                if isempty(errs{1, i})
                    continue;
                end
                testCase.verifyEqual(size(errs{1, i}), [4, 4], names{i});
                for k = 2 : numel(threads)
                    testCase.verifyEqual(DeconAccelTest.bits(errs{k, i}), DeconAccelTest.bits(errs{1, i}), ...
                        sprintf('%s: %d and 1 threads differ.', names{i}, threads(k)));
                end
            end
        end

        function psfNormalizationDoesNotDependOnTheThreadCount(testCase)
            % a PSF of 63^3 elements: RLdecon's normalization, the OTF of
            % decon_lucy_function and decon_lucy_omw_function, and the back
            % projector (normalization and the sum of the OTF magnitudes)
            testCase.assumeGreaterThanOrEqual(feature('numcores'), 2, 'Needs at least 2 cores.');
            nt = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(nt));
            psf = pktest.deconPSF([63, 63, 63], testCase.PSFSigma);
            psf = psf .* (1 + 0.1 * rand(size(psf), 'single'));
            [im, ~] = testCase.texturedData([64, 64, 64]);
            psfFn = [testCase.createTempFolder(), 'psf.tif'];
            writetiff(psf, psfFn);
            threads = [1, min(4, feature('numcores'))];
            outs = cell(2, 3);
            for k = 1 : 2
                maxNumCompThreads(threads(k));
                clear decon_lucy_function decon_lucy_omw_function
                [~, outs{k, 1}] = evalc(['RLdecon('''', '''', psfFn, 0.108, 0.2, 0.2, ''rawdata'', im, ''RLMethod'', ''simplified'', ', ...
                    '''psfGen'', false, ''useGPU'', false, ''Background'', 100, ''DeconIter'', 2, ''save16bit'', false)']);
                [~, outs{k, 2}] = evalc('omw_backprojector_generation(psf, 0.005, false)');
                outs{k, 3} = decon_lucy_omw_function(im, psf, outs{k, 2}, 1, 'Background', testCase.Background, 'useGPU', false, ...
                    'save16bit', false);
            end
            names = {'RLdecon', 'back projector', 'OMW'};
            for i = 1 : 3
                testCase.verifyEqual(DeconAccelTest.bits(outs{2, i}), DeconAccelTest.bits(outs{1, i}), ...
                    sprintf('%s: %d and %d threads differ.', names{i}, threads));
            end
        end
    end

    methods (Access = private)
        function [im, psf] = texturedData(testCase, sz)
            % a smooth texture blurred by the test PSF, on the background, with
            % Poisson-like noise
            psf = pktest.deconPSF([31, 31, 31], testCase.PSFSigma);
            obj = pktest.texture(sz, 2, [0, 2000]);
            im = convn(obj, psf / sum(psf(:)), 'same') + testCase.Background;
            im = uint16(im + sqrt(im) .* randn(sz, 'single'));
        end

        function [J_2] = matlabLoop(testCase, im, psf, nIter, useMex)
            % decon_lucy_function's CPU loop as MATLAB code, with the sums of
            % decon_sum_double
            sz = size(im);
            psf = single(psf);
            psf = psf ./ decon_sum_double(psf);
            I = max(single(im) - testCase.Background, 0);
            if useMex
                H = decon_lucy_rfft(false, 'otf', psf, double(sz));
            else
                H = decon_psf2otf(psf, double(sz));
            end
            J_2 = I;
            J_3 = zeros(sz, 'single');
            J_4 = zeros(prod(sz), 1, 'single');
            lambda = 0;
            for k = 1 : nIter
                if k > 2
                    lambda = decon_sum_double((J_2(:) - Y(:)) .* J_4) / (decon_sum_double(J_4 .* J_4) + eps);
                    lambda = max(min(lambda, 1), 0);
                    J_4 = J_2(:) - Y(:);
                elseif k == 2
                    J_4 = J_2(:) - Y(:);
                end
                Y = max(J_2 + lambda * (J_2 - J_3), 0);
                J_3 = J_2;
                if useMex
                    J_2 = decon_lucy_rfft(false, 'rlstep', Y, I, H);
                else
                    R = cast(I, 'like', J_2) ./ max(real(ifftn(H .* fftn(Y))), eps);
                    J_2 = max(Y .* real(ifftn(conj(H) .* fftn(R))), 0);
                end
            end
        end
    end

    methods (Static, Access = private)
        function [b] = bits(x)
            % the bits of x, to compare -0 and NaN too
            x = x(:);
            if isa(x, 'single')
                b = typecast(x, 'uint32');
            else
                b = typecast(double(x), 'uint64');
            end
        end

        function [x] = wideRange(n, cls)
            % random values with exponents from -40 to 40
            x = cast(randn(n, 1) .* 2 .^ randi([-40, 40], n, 1), cls);
        end

        function [xs] = specialVectors(n, cls)
            z = -zeros(n, 1, cls);
            p = z;
            p(4097) = 0;
            x = DeconAccelTest.wideRange(n, cls);
            xn = x;
            xn(5000) = NaN;
            xi = x;
            xi(77) = Inf;
            xi(9000) = Inf;
            xm = x;
            xm(3) = -Inf;
            % subnormal numbers (and some zeros)
            xd = cast(randn(n, 1), cls) * realmin(cls) / 64;
            xs = {z, p, zeros(n, 1, cls), xn, xi, xm, xd, x};
        end

        function [J_2, Y, J_4, J_3] = accelInputs(n, cls, nonFinite)
            % values in the range of the iterations (J_2(:) - Y(:) correlated
            % with J_4, so that lambda is inside (0, 1)), with -0, +0 and
            % subnormal numbers among them, and NaN, Inf and -Inf if nonFinite
            J_2 = cast(rand(n, 1) * 1000, cls);
            J_3 = cast(rand(n, 1) * 1000, cls);
            J_4 = cast(randn(n, 1) * 10, cls);
            Y = J_2 - 0.6 * J_4 + cast(randn(n, 1), cls);
            sp = cast([0, -0, NaN, Inf, -Inf, 1e-40, -1e-40, 5, -0, 0], cls);
            if ~nonFinite
                sp(3 : 5) = cast([0, -0, 3], cls);
            end
            J_2(1 : 10) = sp;
            J_3(1 : 10) = sp([2, 1, 4, 5, 3, 1, 2, 9, 2, 1]);
            Y(11 : 20) = sp;
            J_4(21 : 30) = sp;
            J_2(n - 9 : n) = sp;
            J_3(n - 9 : n) = -0;
        end

        function [tf] = rfftMexDefault()
            % decon_lucy_rfft_mex is used where it gives MATLAB's FFTs (Linux)
            tf = exist('decon_lucy_rfft_mex', 'file') == 3 && isunix && ~ismac;
        end
    end
end
