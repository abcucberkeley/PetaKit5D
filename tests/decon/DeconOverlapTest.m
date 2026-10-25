classdef DeconOverlapTest < pktest.TestCase
    % with overlapFcn, decon_lucy_omw_function and decon_lucy_function run a
    % function on the CPU while the GPU computes (XR_RLdeconFrame3D reads
    % the next file of a GPU job then) and return its output, with the same
    % results

    properties (Constant)
        Background = 100
        ImSize = [40, 44, 36]
    end

    properties (SetAccess = private)
        PSF
        BackProjector
    end

    methods (TestClassSetup)
        function createPSFAndBackProjector(testCase)
            testCase.PSF = pktest.deconPSF([31, 31, 31], [1.5, 2]);
            [~, testCase.BackProjector] = evalc('omw_backprojector_generation(testCase.PSF, 0.005, false)');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function overlapFunctionRunsOnceAndReturnsItsOutput(testCase)
            im = testCase.beads();
            calls = containers.Map({'n'}, {0});
            overlap = @() DeconOverlapTest.countCall(calls);
            expected = decon_lucy_omw_function(im, testCase.PSF, testCase.BackProjector, 2, 'Background', testCase.Background, ...
                'useGPU', false, 'save16bit', true);
            [out, ~, ov] = decon_lucy_omw_function(im, testCase.PSF, testCase.BackProjector, 2, 'Background', testCase.Background, ...
                'useGPU', false, 'save16bit', true, 'overlapFcn', overlap);
            testCase.verifyEqual(out, expected);
            testCase.verifyEqual(ov, 42);
            testCase.verifyEqual(calls('n'), 1);
            expected = decon_lucy_function(im, testCase.PSF, 3, 'Background', testCase.Background, 'useGPU', false, 'save16bit', true);
            [out, ~, ~, ov] = decon_lucy_function(im, testCase.PSF, 3, 'Background', testCase.Background, ...
                'useGPU', false, 'save16bit', true, 'overlapFcn', overlap);
            testCase.verifyEqual(out, expected);
            testCase.verifyEqual(ov, 42);
            testCase.verifyEqual(calls('n'), 2);
            % without a function, the output is empty
            [~, ~, ov] = decon_lucy_omw_function(im, testCase.PSF, testCase.BackProjector, 2, 'useGPU', false);
            testCase.verifyEmpty(ov);
            [~, ~, ~, ov] = decon_lucy_function(im, testCase.PSF, 3, 'useGPU', false);
            testCase.verifyEmpty(ov);
        end

        function rldeconReturnsTheOverlapOutput(testCase)
            % as the large-file batches call RLdecon (data in memory, no
            % result file)
            im = testCase.beads();
            folder = testCase.createTempFolder();
            psfFn = [folder, 'psf.tif'];
            writetiff(testCase.PSF, psfFn);
            mkdir([folder, 'out']);
            outFn = [folder, 'out/result.tif'];
            args = {'rawdata', im, 'DeconIter', 2, 'Background', testCase.Background, 'useGPU', false, 'psfGen', false, ...
                'skewed', false, 'save3Dstack', [false, false, false], 'mipAxis', [0, 0, 0]};
            for method = {'omw', 'simplified'}
                [~, out, ~, ~, ov] = evalc('RLdecon('''', outFn, psfFn, 0.108, 0.2, 0.2, args{:}, ''RLMethod'', method{1}, ''overlapFcn'', @() ''next'')');
                testCase.verifyEqual(ov, 'next', method{1});
                [~, expected] = evalc('RLdecon('''', outFn, psfFn, 0.108, 0.2, 0.2, args{:}, ''RLMethod'', method{1})');
                testCase.verifyEqual(out, expected, method{1});
            end
        end
    end

    methods (Test, TestTags = {'Unit', 'GPU'})
        function gpuOverlapFunctionRunsWhileTheGPUComputes(testCase)
            testCase.assumeGPUAvailable();
            im = testCase.beads();
            expected = decon_lucy_omw_function(im, testCase.PSF, testCase.BackProjector, 2, 'Background', testCase.Background, ...
                'useGPU', true);
            [out, ~, ov] = decon_lucy_omw_function(im, testCase.PSF, testCase.BackProjector, 2, 'Background', testCase.Background, ...
                'useGPU', true, 'overlapFcn', @() 'read');
            testCase.verifyEqual(out, expected, 'omw');
            testCase.verifyEqual(ov, 'read', 'omw');
            expected = decon_lucy_function(im, testCase.PSF, 5, 'Background', testCase.Background, 'useGPU', true);
            [out, ~, ~, ov] = decon_lucy_function(im, testCase.PSF, 5, 'Background', testCase.Background, ...
                'useGPU', true, 'overlapFcn', @() 'read');
            testCase.verifyEqual(out, expected, 'simplified');
            testCase.verifyEqual(ov, 'read', 'simplified');
        end
    end

    methods (Static)
        function [out] = countCall(calls)
            calls('n') = calls('n') + 1;
            out = 42;
        end
    end

    methods (Access = private)
        function [im] = beads(testCase)
            sz = testCase.ImSize;
            pos = [4 + (sz(1) - 8) * rand(5, 1), 4 + (sz(2) - 8) * rand(5, 1), 4 + (sz(3) - 8) * rand(5, 1)];
            im = uint16(pktest.beads(sz, pos, sqrt(0.8 ^ 2 + [1.5, 2] .^ 2), 1000, testCase.Background));
        end
    end
end
