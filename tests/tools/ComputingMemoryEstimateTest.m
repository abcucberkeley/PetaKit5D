classdef ComputingMemoryEstimateTest < pktest.TestCase
    % XR_estimateComputingMemory: memory (GB) for each processing step of
    % an image, rawImageSize * memFactors(step) (deskew: also times
    % (nz / 300) ^ 2 above 300 slices), and the GPU memory of deconvolution

    methods (Test, TestTags = {'Unit'})
        function eachStepOfSeveralGetsItsEstimate(testCase)
            % a step was estimated only if every step matched it, so with
            % the default three steps every estimate was 0 and the GPU
            % estimate was not assigned
            sz = [100, 200, 600];
            raw = prod(sz) * 4 / 1024 ^ 3;
            [est, gpuEst, rawSize, imSize] = XR_estimateComputingMemory('', 'imSize', sz);
            testCase.verifyEqual(est, [raw * 15 * 4; raw * 5; raw * 15], 'RelTol', 1e-12);
            testCase.verifyEqual(gpuEst, raw * 1.5, 'RelTol', 1e-12);
            testCase.verifyEqual(rawSize, raw, 'RelTol', 1e-12);
            testCase.verifyEqual(imSize, sz);

            est = XR_estimateComputingMemory('', {'rotate', 'deskew'}, 'imSize', sz);
            testCase.verifyEqual(est, [raw * 5; raw * 15 * 4], 'RelTol', 1e-12);
        end

        function dataSizeWithSeveralSteps(testCase)
            % with dataSize, several steps failed at || (a condition of
            % two values)
            est = XR_estimateComputingMemory('', {'rotate', 'deconvolution'}, 'dataSize', 2 ^ 30, 'cudaDecon', false);
            testCase.verifyEqual(est, [2 * 5; 2 * 15], 'RelTol', 1e-12);
            % with deskew among them, the image size counts
            sz = [100, 200, 600];
            raw = prod(sz) * 4 / 1024 ^ 3;
            est = XR_estimateComputingMemory('', {'deskew', 'rotate'}, 'imSize', sz, 'dataSize', 2 ^ 30);
            testCase.verifyEqual(est, [raw * 15 * 4; raw * 5], 'RelTol', 1e-12);
        end

        function singleStepsAsTheCallersPassThem(testCase)
            % the callers pass one step: their estimates stay as they were
            sz = [64, 48, 301];
            raw = prod(sz) * 4 / 1024 ^ 3;
            testCase.verifyEqual(XR_estimateComputingMemory('', 'steps', {'deskew'}, 'imSize', sz), ...
                raw * 15 * (301 / 300) ^ 2, 'RelTol', 1e-12);
            testCase.verifyEqual(XR_estimateComputingMemory('', {'rotate'}, 'imSize', sz), raw * 5, 'RelTol', 1e-12);
            [est, gpuEst] = XR_estimateComputingMemory('', {'deconvolution'}, 'imSize', sz, 'cudaDecon', true);
            testCase.verifyEqual([est, gpuEst], [raw * 15, raw * 1.5], 'RelTol', 1e-12);
            [~, gpuEst] = XR_estimateComputingMemory('', {'deconvolution'}, 'imSize', sz, 'cudaDecon', false);
            testCase.verifyTrue(isnan(gpuEst));

            % a tiff file, as XR_cppDeconFrame3D passes it
            fn = [testCase.createTempFolder(), 'im.tif'];
            writetiff(zeros(16, 12, 5, 'uint16'), fn);
            [est, ~, rawSize, imSize] = XR_estimateComputingMemory(fn, {'deconvolution'}, 'cudaDecon', false);
            testCase.verifyEqual(imSize, [16, 12, 5]);
            testCase.verifyEqual(rawSize, 16 * 12 * 5 * 4 / 1024 ^ 3, 'RelTol', 1e-12);
            testCase.verifyEqual(est, rawSize * 15, 'RelTol', 1e-12);
        end
    end
end
