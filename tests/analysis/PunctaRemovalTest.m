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

        function parametersFileHoldsTheParameters(testCase)
            % parameters.json was written with its text as fprintf's format:
            % a % in a value (here in the folder name) cut the file off there
            [zarrFn, vol] = testCase.punctaVolume('data_100%/');
            XR_ExM_PunctaRemoval_zarr(zarrFn, 'batchSize', [64, 64, 64], 'blockSize', [32, 32, 32], ...
                'parseCluster', false, 'masterCompute', true);
            testCase.verifyPunctaRemoved(zarrFn, vol);
            pr = jsondecode(fileread([fileparts(zarrFn), '/Puncta_Removed/parameters.json']));
            testCase.verifyEqual(pr.zarrFullpath, zarrFn);
        end

        function interpolatedBackgroundOneValueDeepInZ(testCase)
            % a background estimated in batches as deep as the volume is one
            % value deep in z, which MATLAB's arrays drop: imresize3 of a
            % batch's background failed (MATLAB:badsubscript). The plane's
            % interpolation is the background at every z, the same as for the
            % plane twice in z
            [zarrFn, vol] = testCase.punctaVolume('', [64, 56, 48]);
            dataPath = fileparts(zarrFn);
            local = {'blockSize', [16, 16, 16], 'parseCluster', false, 'masterCompute', true};
            XR_ExM_PunctaRemoval_zarr(zarrFn, 'mode', 'bg_estimation', 'batchSize', [32, 32, 48], local{:});
            bgPath = [dataPath, '/Puncta_background_estimation/'];
            bgFn = [bgPath, 'expanded_batch_32_32_48_background_values.zarr'];
            bg = readzarr(bgFn);
            testCase.assertSize(bg, [2, 2]);
            twiceFn = [bgPath, 'twice_batch_32_32_24_background_values.zarr'];
            writezarr(cat(3, bg, bg), twiceFn, blockSize=[2, 2, 2]);
            outPath = [dataPath, '/Puncta_Removed/'];
            XR_ExM_PunctaRemoval_zarr(zarrFn, 'mode', 'bg_interpolated', 'bgFullpath', twiceFn, 'batchSize', [64, 64, 48], local{:});
            expected = readzarr([outPath, 'expanded.zarr']);
            rmdir(outPath, 's');

            XR_ExM_PunctaRemoval_zarr(zarrFn, 'mode', 'bg_interpolated', 'bgFullpath', bgFn, 'batchSize', [64, 64, 48], local{:});
            testCase.verifyEqual(readzarr([outPath, 'expanded.zarr']), expected);
            testCase.verifyPunctaRemoved(zarrFn, vol);
        end
    end

    methods (Test, TestTags = {'Unit'})
        function punctumAtACornerIsFilledFromTheNearestPixels(testCase)
            % the linear fill (fillmissing2) left the removed pixels of a
            % punctum at a corner of the image NaN, outside the convex hull
            % of the known pixels, and a later pass on that volume failed in
            % imregionalmax (MATLAB:imregionalmax:expectedNonNaN); they take
            % the value of the nearest pixel, so the fills stay within the
            % values of the pixels that are kept
            [y, x, z] = ndgrid(0 : 19, 0 : 20, 0 : 15);
            vol = single(300 + 2 * y + 1.5 * x + 3 * z) + pktest.beads([20, 21, 16], [1, 1, 9], 1.5, 20000, 0);
            args = {'localWinSize', [15, 15, 23], 'intThrsh', 3000};
            [~, out] = evalc('XR_ExM_PunctaRemovalPointDetection(vol, args{:})');
            testCase.verifyFalse(any(isnan(out(:))), 'The fill left NaN.');
            testCase.verifyLessThan(out(1, 1, 9), 0.1 * 20000, 'The punctum was not removed.');
            filled = out ~= vol;
            testCase.assertTrue(any(filled(:)), 'Nothing was filled.');
            kept = vol(~filled);
            testCase.verifyGreaterThanOrEqual(min(out(filled)), min(kept));
            testCase.verifyLessThanOrEqual(max(out(filled)), max(kept));
            [~, out2] = evalc('XR_ExM_PunctaRemovalPointDetection(out, args{:})');
            testCase.verifySize(out2, size(vol));
        end
    end

    methods (Access = private)
        function [zarrFn, vol] = punctaVolume(testCase, folder, sz)
            % the volume in a folder of a temporary folder (folder: its name),
            % of size sz (default [64, 64, 64])
            dataPath = testCase.createTempFolder();
            if nargin > 1 && ~isempty(folder)
                dataPath = [dataPath, folder];
                mkdir(dataPath);
            end
            if nargin < 3
                sz = [64, 64, 64];
            end
            vol = 100 + pktest.beads(sz, testCase.Puncta, 1.5, testCase.PunctumAmplitude, 0);
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
