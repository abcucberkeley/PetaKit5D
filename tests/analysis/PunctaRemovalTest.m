classdef PunctaRemovalTest < pktest.TestCase
    % XR_ExM_PunctaRemoval_zarr (expansion microscopy puncta removal): a
    % smooth stained background with a few very bright puncta. The puncta
    % must be removed and the output must keep the input size.

    properties (Constant)
        Puncta = [20, 20, 20; 40, 44, 30; 30, 50, 45]
        PunctumAmplitude = 20000
    end

    methods (Test, TestTags = {'Integration', 'KnownIssue'})
        function removesPunctaWithDefaultThresholds(testCase)
            % KNOWN ISSUE: XR_ExM_PunctaRemoval_zarr.m builds the 'updated'
            % (default) task string with scalar %0.20d conversions for
            % MinThreshold, MaxThreshold and BaseThreshold. The defaults
            % ([0, 0], [] and [], i.e., Otsu thresholds) shift all later
            % sprintf arguments, so the task string is invalid MATLAB
            % ("Invalid expression") and the function always fails with
            % default settings. Fix: format them with mat2str.
            [zarrFn, vol] = testCase.punctaVolume();
            XR_ExM_PunctaRemoval_zarr(zarrFn, 'batchSize', [64, 64, 64], 'blockSize', [32, 32, 32], ...
                'parseCluster', false, 'masterCompute', true);
            testCase.verifyPunctaRemoved(zarrFn, vol);
        end

        function removesPunctaWithScalarThresholds(testCase)
            % KNOWN ISSUE: even with scalar thresholds, the 'updated' task
            % string passes 'uuid' to XR_ExM_PunctaRemoval_updated_block,
            % which has no uuid parameter ("'uuid' is not a recognized
            % parameter"). Fix: add the parameter to the block function or
            % drop it from the task string.
            [zarrFn, vol] = testCase.punctaVolume();
            XR_ExM_PunctaRemoval_zarr(zarrFn, 'batchSize', [64, 64, 64], 'blockSize', [32, 32, 32], ...
                'MinThreshold', 0, 'MaxThreshold', 2000, 'BaseThreshold', 300, 'volThrsh', 100, ...
                'parseCluster', false, 'masterCompute', true);
            testCase.verifyPunctaRemoved(zarrFn, vol);
        end
    end

    methods (Access = private)
        function [zarrFn, vol] = punctaVolume(testCase)
            dataPath = testCase.createTempFolder();
            vol = pktest.texture([64, 64, 64], 4, [500, 1000]) + ...
                pktest.beads([64, 64, 64], testCase.Puncta, 1.5, testCase.PunctumAmplitude, 0);
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
        end
    end
end
