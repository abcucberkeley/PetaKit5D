classdef PunctaRemovalTest < pktest.TestCase
    % XR_ExM_PunctaRemoval_zarr (expansion microscopy puncta removal): a
    % large stained structure and a few small, very bright puncta on a dark
    % background. The puncta must be removed, the structure kept, and the
    % output must keep the input size.

    properties (Constant)
        Puncta = [20, 20, 20; 40, 44, 30; 30, 50, 45]
        PunctumAmplitude = 20000
        % y, x, z ranges of the stained structure
        Structure = {45 : 60, 5 : 30, 5 : 35}
    end

    methods (Test, TestTags = {'Integration'})
        function removesPunctaWithDefaultThresholds(testCase)
            % the task string used to be invalid with the default (empty)
            % thresholds, and the block function couldn't handle them
            [zarrFn, vol] = testCase.punctaVolume();
            XR_ExM_PunctaRemoval_zarr(zarrFn, 'batchSize', [64, 64, 64], 'blockSize', [32, 32, 32], ...
                'parseCluster', false, 'masterCompute', true);
            testCase.verifyPunctaRemoved(zarrFn, vol);
        end

        function removesPunctaWithScalarThresholds(testCase)
            % the task string used to pass parameters the block function
            % doesn't take, and the block called a misspelled function
            [zarrFn, vol] = testCase.punctaVolume();
            XR_ExM_PunctaRemoval_zarr(zarrFn, 'batchSize', [64, 64, 64], 'blockSize', [32, 32, 32], ...
                'MinThreshold', 0, 'MaxThreshold', 3000, 'BaseThreshold', 1000, ...
                'parseCluster', false, 'masterCompute', true);
            testCase.verifyPunctaRemoved(zarrFn, vol);
        end
    end

    methods (Access = private)
        function [zarrFn, vol] = punctaVolume(testCase)
            dataPath = testCase.createTempFolder();
            vol = 100 + pktest.beads([64, 64, 64], testCase.Puncta, 1.5, testCase.PunctumAmplitude, 0);
            st = testCase.Structure;
            vol(st{:}) = vol(st{:}) + 4900 + 200 * randn(cellfun(@numel, st));
            zarrFn = [dataPath, 'expanded.zarr'];
            writezarr(single(vol), zarrFn, blockSize=[32, 32, 32]);
        end

        function [] = verifyPunctaRemoved(testCase, zarrFn, vol)
            outFn = [fileparts(zarrFn), '/Puncta_Removed/expanded.zarr'];
            testCase.assertTrue(exist(outFn, 'dir') == 7, 'No puncta removal output.');
            out = readzarr(outFn);
            testCase.verifySize(out, size(vol));
            for i = 1 : size(testCase.Puncta, 1)
                p = testCase.Puncta(i, :);
                testCase.verifyLessThan(out(p(1), p(2), p(3)), 0.1 * testCase.PunctumAmplitude, ...
                    sprintf('Punctum at %s was not removed.', mat2str(p)));
            end
            % the interior of the structure (its edges are smoothed)
            inner = cellfun(@(r) r(4 : end - 3), testCase.Structure, 'unif', 0);
            testCase.verifyEqual(mean(out(inner{:}), 'all'), mean(vol(inner{:}), 'all'), 'RelTol', 0.01, ...
                'The stained structure was not kept.');
        end
    end
end
