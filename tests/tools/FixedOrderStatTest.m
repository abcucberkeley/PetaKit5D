classdef FixedOrderStatTest < pktest.TestCase
    % fixed_order_stat: MATLAB's sum, mean, var and std of all elements in
    % its order at one thread, bit for bit (also -0), at any number of
    % threads (here 1, 2, 4 and 24). MATLAB splits sums of 89000 or more
    % elements among its threads, so its own results there change in the
    % last bits with the thread count (in 778 to 821 of 1476 cases of a
    % probe on Linux at 2 to 24 threads); below 89000 elements, and for other
    % arrays and options, fixed_order_stat is MATLAB's own function.

    properties (Constant)
        Threads = [1, 2, 4, 24]
    end

    properties (TestParameter)
        dataClass = {'single', 'double'}
        statCall = struct('sum', {{'sum'}}, 'sumDouble', {{'sum', 'double'}}, 'mean', {{'mean'}}, ...
            'meanDouble', {{'mean', 'double'}}, 'var', {{'var'}}, 'var0', {{'var', 0}}, 'var1', {{'var', 1}}, ...
            'std', {{'std'}}, 'std1', {{'std', 1}})
        % values over a wide range of exponents, and large values with a small
        % spread (single sums of them cancel), whose sums other orders round
        % otherwise
        dataKind = {'wideRange', 'offset'}
    end

    methods (Test, TestTags = {'Unit'})
        function isMatlabsAtOneThread(testCase, dataClass, statCall, dataKind)
            % around MATLAB's switches (44 elements: lanes, 89000: blocks of
            % 4096), with a last block that isn't full and a last group of four
            % lanes that isn't either, and a volume: from 89000 elements the
            % same at every thread count, on Linux MATLAB's order at one thread
            % (as found for R2024b); below, MATLAB's own sums, which don't
            % depend on the thread count on Linux and macOS (as measured)
            nt = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(nt));
            sizes = {1, 43, 44, 88999, 89000, 89001, 4096 * 30, 4096 * 30 + 3, [61, 47, 89]};
            for i = 1 : numel(sizes)
                x = FixedOrderStatTest.data(dataKind, sizes{i}, dataClass);
                name = sprintf('%s of %s', statCall{1}, mat2str(sizes{i}));
                maxNumCompThreads(1);
                ref = feval(statCall{1}, x(:), statCall{2 : end});
                v = cell(1, numel(testCase.Threads));
                for k = 1 : numel(testCase.Threads)
                    maxNumCompThreads(testCase.Threads(k));
                    v{k} = fixed_order_stat(statCall{1}, x, statCall{2 : end});
                    if numel(x) >= 89000 || isunix
                        testCase.verifyEqual(FixedOrderStatTest.bits(v{k}), FixedOrderStatTest.bits(v{1}), ...
                            sprintf('%s: %d threads differ from 1.', name, testCase.Threads(k)));
                    end
                end
                if numel(x) < 89000 || (isunix && ~ismac)
                    testCase.verifyEqual(FixedOrderStatTest.bits(v{1}), FixedOrderStatTest.bits(ref), ...
                        sprintf('%s: not MATLAB''s at one thread.', name));
                end
                testCase.verifyClass(v{1}, class(ref));
            end
        end

        function signedZerosAndNonFiniteValues(testCase, dataClass)
            % a sum of -0s is -0 (also the mean), a +0 among them +0; NaN and
            % Inf give NaN or Inf as in MATLAB (above 89000 elements)
            n = 250047;
            z = -zeros(n, 1, dataClass);
            testCase.verifyEqual(FixedOrderStatTest.bits(fixed_order_stat('sum', z)), FixedOrderStatTest.bits(-zeros(1, dataClass)));
            testCase.verifyEqual(FixedOrderStatTest.bits(fixed_order_stat('sum', z, 'double')), FixedOrderStatTest.bits(-0));
            testCase.verifyEqual(FixedOrderStatTest.bits(fixed_order_stat('mean', z)), FixedOrderStatTest.bits(-zeros(1, dataClass)));
            testCase.verifyEqual(fixed_order_stat('var', z), zeros(1, dataClass));
            p = z;
            p(4097) = 0;
            testCase.verifyEqual(FixedOrderStatTest.bits(fixed_order_stat('sum', p)), FixedOrderStatTest.bits(zeros(1, dataClass)));
            x = FixedOrderStatTest.data('wideRange', n, dataClass);
            x(5000) = NaN;
            for f = {'sum', 'mean', 'var', 'std'}
                testCase.verifyTrue(isnan(fixed_order_stat(f{1}, x)), [f{1}, ' with NaN']);
            end
            x(5000) = Inf;
            testCase.verifyEqual(fixed_order_stat('sum', x), cast(Inf, dataClass));
            testCase.verifyTrue(isnan(fixed_order_stat('var', x)), 'var with Inf');
            if isunix && ~ismac
                nt = maxNumCompThreads;
                testCase.addTeardown(@() maxNumCompThreads(nt));
                maxNumCompThreads(1);
                testCase.verifyEqual(FixedOrderStatTest.bits(fixed_order_stat('mean', z)), FixedOrderStatTest.bits(mean(z)));
                testCase.verifyEqual(FixedOrderStatTest.bits(fixed_order_stat('std', p, 1)), FixedOrderStatTest.bits(std(p, 1)));
            end
        end

        function callersAskForMatlabsExpressions(testCase)
            % the registration's NCC, the crop of its region and decon's debug
            % error metric take MATLAB's expressions of before (mean(T(:),
            % 'double'), ..., var(A(:), 1)) from fixed_order_stat: a spy records
            % its calls and returns MATLAB's own value; at one thread their
            % results are the same as with fixed_order_stat (on Linux, where
            % its order is MATLAB's)
            nt = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(nt));
            maxNumCompThreads(1);
            T = single(rand([50, 45, 41]) * 100);
            A = single(rand([70, 66, 50]) * 100);
            bbox = [30, 25, 20, 90, 85, 70];
            region = uint16(randi([100, 120], [40, 40, 120]));
            region(20 : 22, 10 : 12, 80 : 82) = 5000;
            psf = pktest.deconPSF([15, 15, 15], [1.5, 2]);
            im = uint16(100 + 1000 * rand([48, 48, 40]));
            root = testCase.createTempFolder();
            compute = @() {normxcorr3_fast(T, A, 'full', bbox), ...
                FixedOrderStatTest.nthOut(2, @crop_subregion_by_intensity, region, 50), ...
                FixedOrderStatTest.nthOut(2, @decon_lucy_function, im, psf, 3, 'Background', 100, 'useGPU', false, 'debug', true, ...
                'debug_folder', [root, sprintf('rl_%d', numel(dir(root)))])};
            ref = compute();
            spyDir = [root, 'spy/'];
            mkdir(spyDir);
            fid = fopen([spyDir, 'fixed_order_stat.m'], 'w');
            fprintf(fid, '%s\n', 'function [v] = fixed_order_stat(fun, x, opt)', ...
                '% test spy written by FixedOrderStatTest: records the calls, returns MATLAB''s own value', ...
                'persistent calls', 'if strcmp(fun, ''-calls'')', '    v = calls;', '    calls = {};', '    return;', 'end', ...
                'if nargin < 3', '    opt = [];', 'end', 'calls(end + 1, :) = {fun, class(x), numel(x), opt};', ...
                'if isempty(opt)', '    v = feval(fun, x(:));', 'else', '    v = feval(fun, x(:), opt);', 'end', 'end');
            fclose(fid);
            addpath(spyDir, '-begin');
            testCase.addTeardown(@() clear('fixed_order_stat'));
            testCase.addTeardown(@() rmpath(spyDir));
            clear fixed_order_stat normxcorr3_fast crop_subregion_by_intensity decon_lucy_function
            fixed_order_stat('-calls');
            out = compute();
            calls = fixed_order_stat('-calls');
            nT = numel(T);
            nA = numel(A);
            expected = [{'mean', 'single', nT, 'double'}; {'mean', 'single', nA, 'double'}; {'sum', 'single', nT, 'double'}; ...
                {'std', 'double', nT, []}; {'var', 'single', nA, 1}; {'mean', 'uint16', numel(region), []}; ...
                {'std', 'single', numel(region), []}; repmat({'mean', 'double', numel(im), []}, 3, 1)];
            testCase.verifyEqual(calls, expected);
            if isunix && ~ismac
                names = {'normxcorr3_fast', 'crop_subregion_by_intensity', 'decon_lucy_function''s err_mat'};
                for i = 1 : 3
                    testCase.verifyEqual(FixedOrderStatTest.bits(out{i}), FixedOrderStatTest.bits(ref{i}), names{i});
                end
            end
        end

        function otherArraysAndOptionsAreMatlabs(testCase)
            % integer (exact in double), logical and complex data, options it
            % doesn't take: MATLAB's own function of x(:)
            n = 250047;
            u = uint16(randi([0, 65535], [61, 41, 100]));
            testCase.verifyEqual(fixed_order_stat('mean', u), mean(u(:)));
            testCase.verifyClass(fixed_order_stat('mean', u), 'double');
            testCase.verifyEqual(fixed_order_stat('sum', u), sum(u(:)));
            b = rand(n, 1) > 0.5;
            testCase.verifyEqual(fixed_order_stat('sum', b), sum(b));
            c = complex(rand(n, 1, 'single'), rand(n, 1, 'single'));
            testCase.verifyEqual(fixed_order_stat('sum', c), sum(c));
            x = rand(n, 1, 'single');
            testCase.verifyEqual(fixed_order_stat('sum', x, 'native'), sum(x, 'native'));
            w = rand(n, 1);
            testCase.verifyEqual(fixed_order_stat('var', x, w), var(x, w));
            testCase.verifyError(@() fixed_order_stat('median', x), ?MException);
        end
    end

    methods (Static, Access = private)
        function [b] = bits(x)
            % the bits of x, to compare -0 too (NaN as NaN)
            x = x(:);
            if all(isnan(x))
                b = 'NaN';
            elseif isa(x, 'single')
                b = typecast(x, 'uint32');
            else
                b = typecast(double(x), 'uint64');
            end
        end

        function [v] = nthOut(n, f, varargin)
            % the n-th output of f(varargin{:})
            out = cell(1, n);
            [out{:}] = f(varargin{:});
            v = out{n};
        end

        function [x] = data(kind, sz, cls)
            if isscalar(sz)
                sz = [sz, 1];
            end
            switch kind
                case 'wideRange'
                    % exponents from -40 to 40
                    x = cast(randn(sz) .* 2 .^ randi([-40, 40], sz), cls);
                case 'offset'
                    x = cast(1e7 + round(16 * randn(sz)), cls);
            end
        end
    end
end
