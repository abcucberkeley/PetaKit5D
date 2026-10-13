classdef ResampleTest < pktest.TestCase
    % XR_resample_dataset / XR_resampleFrame / XR_resampleSingleZarr.
    % resampleFactor is in yxz order; a factor of 2 downsamples 2x.

    properties (Constant)
        ImSize = [40, 36, 12]
        Channels = {'CamA_ch0', 'CamB_ch0'}
    end

    properties (TestParameter)
        resampleFactor = struct('down2', [2, 2, 2], 'zOnly', [1, 1, 2], 'up2xy', [0.5, 0.5, 1], ...
            'nonInteger', [1.5, 1.5, 1.5])
        interpMethod = {'linear', 'cubic', 'nearest', 'max', 'mean'}
    end

    methods (Test, TestTags = {'Unit'})
        function linearMatchesImresize3(testCase, resampleFactor)
            [im, fn, dataPath] = testCase.writeImage('tif');
            fnout = [dataPath, 'out.tif'];
            XR_resampleFrame(fn, fnout, resampleFactor);
            outSize = round(testCase.ImSize ./ resampleFactor);
            expected = uint16(imresize3(single(im), outSize, 'Method', 'linear'));
            testCase.verifyEqual(readtiff(fnout), expected);
        end

        function constantImageStaysConstant(testCase, interpMethod)
            dataPath = testCase.createTempFolder();
            fn = [dataPath, 'im.tif'];
            writetiff(1000 * ones(testCase.ImSize, 'uint16'), fn);
            fnout = [dataPath, 'out.tif'];
            XR_resampleFrame(fn, fnout, [2, 3, 2], 'interpMethod', interpMethod);
            out = readtiff(fnout);
            if any(strcmp(interpMethod, {'max', 'mean'}))
                outSize = ceil(testCase.ImSize ./ [2, 3, 2]);
            else
                outSize = round(testCase.ImSize ./ [2, 3, 2]);
            end
            testCase.verifyEqual(size(out), outSize);
            testCase.verifyEqual(unique(out), uint16(1000));
        end

        function maxAndMeanMatchBlockReference(testCase)
            % non-divisible size: edge blocks use only the voxels inside
            [im, fn, dataPath] = testCase.writeImage('tif');
            rs = [3, 5, 5];
            XR_resampleFrame(fn, [dataPath, 'max.tif'], rs, 'interpMethod', 'max', 'save16bit', false);
            XR_resampleFrame(fn, [dataPath, 'mean.tif'], rs, 'interpMethod', 'mean', 'save16bit', false);
            [blockMax, blockMean] = testCase.blockStats(single(im), rs);
            testCase.verifyEqual(readtiff([dataPath, 'max.tif']), blockMax);
            testCase.verifyEqual(readtiff([dataPath, 'mean.tif']), blockMean, 'RelTol', 1e-6);
        end

        function image2d(testCase)
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, 40, 36));
            writetiff(im, [dataPath, 'im.tif']);
            XR_resampleFrame([dataPath, 'im.tif'], [dataPath, 'out.tif'], [2, 2]);
            expected = uint16(imresize(single(im), [20, 18], 'Method', 'bilinear'));
            testCase.verifyEqual(readtiff([dataPath, 'out.tif']), expected);
        end

        function inputBboxResamplesTheCroppedRegion(testCase)
            [im, fn, dataPath] = testCase.writeImage('tif');
            b = [5, 7, 3, 30, 33, 10];
            XR_resampleFrame(fn, [dataPath, 'out.tif'], [2, 2, 2], 'inputBbox', b);
            region = single(im(b(1) : b(4), b(2) : b(5), b(3) : b(6)));
            expected = uint16(imresize3(region, round(size(region) ./ 2), 'Method', 'linear'));
            testCase.verifyEqual(readtiff([dataPath, 'out.tif']), expected);
        end

        function singleOutputWithoutSave16bit(testCase)
            [im, fn, dataPath] = testCase.writeImage('tif');
            XR_resampleFrame(fn, [dataPath, 'out.tif'], [2, 2, 2], 'save16bit', false);
            out = readtiff([dataPath, 'out.tif']);
            testCase.verifyClass(out, 'single');
            testCase.verifyEqual(out, imresize3(single(im), round(testCase.ImSize / 2), 'Method', 'linear'));
        end

        function zarrInputAndOutput(testCase)
            [im, fn, dataPath] = testCase.writeImage('zarr');
            XR_resampleFrame(fn, [dataPath, 'out.zarr'], [2, 2, 2], 'zarrFile', true, 'saveZarr', true);
            expected = uint16(imresize3(single(im), round(testCase.ImSize / 2), 'Method', 'linear'));
            testCase.verifyEqual(readzarr([dataPath, 'out.zarr']), expected);
        end
    end

    methods (Test, TestTags = {'Integration'})
        function datasetWrapperTiff(testCase)
            [dataPath, ims, names] = testCase.writeDataset('tif');
            XR_resample_dataset(dataPath, [2, 2, 2], 'channelPatterns', testCase.Channels, ...
                'parseCluster', false, 'masterCompute', true);
            resultPath = [dataPath, 'resampled/'];
            for i = 1 : 2
                expected = uint16(imresize3(single(ims{i}), round(testCase.ImSize / 2), 'Method', 'linear'));
                testCase.verifyEqual(readtiff([resultPath, names{i}, '.tif']), expected, names{i});
            end
            testCase.verifyFalse(exist([resultPath, names{3}, '.tif'], 'file') > 0, ...
                'A file that matches no channel pattern was resampled.');
            params = jsondecode(fileread([resultPath, 'parameters.json']));
            testCase.verifyEqual(params.resampleFactor(:)', [2, 2, 2]);
        end

        function datasetWrapperZarr(testCase)
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            XR_resample_dataset(dataPath, [2, 2, 2], 'channelPatterns', testCase.Channels, 'zarrFile', true, ...
                'saveZarr', true, 'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                expected = uint16(imresize3(single(ims{i}), round(testCase.ImSize / 2), 'Method', 'linear'));
                testCase.verifyEqual(readzarr([dataPath, 'resampled/', names{i}, '.zarr']), expected, names{i});
            end
        end

        function largeZarrMatchesInMemory(testCase)
            % batch-wise resampling with border buffers vs the whole image.
            % The border is a multiple of the factor, see the known issue
            % largeZarrWithDefaultBorderSize.
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            XR_resample_dataset(dataPath, [2, 2, 2], 'channelPatterns', testCase.Channels(1), 'zarrFile', true, ...
                'largeFile', true, 'saveZarr', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 2], ...
                'borderSize', [4, 4, 4], 'parseCluster', false, 'masterCompute', true);
            out = readzarr([dataPath, 'resampled/', names{1}, '.zarr']);
            expected = uint16(imresize3(single(ims{1}), round(testCase.ImSize / 2), 'Method', 'linear'));
            % +/-1 from rounding uint16 vs single interpolation results
            testCase.verifyImagesClose(out, expected, 'relTol', 1e-4, 'absTol', 1);
        end
    end

    methods (Test, TestTags = {'Integration', 'KnownIssue'})
        function largeZarrWithDefaultBorderSize(testCase)
            % KNOWN ISSUE: resampleZarrBlock.m resizes each batch together with
            % its border and crops the result at baStart = round((ibStart_orig
            % - ibStart) ./ resampleFactor) + 1. If the border is not a
            % multiple of the factor (the default borderSize [5, 5, 5] with a
            % factor of 2), the batch is sampled on a grid shifted by half an
            % output voxel, so large-file resampling differs from in-memory
            % resampling almost everywhere (by up to ~10% of the intensity
            % range for this noise image) and has seams at batch borders.
            % Fix: round the border up to a multiple of the factor, e.g.,
            % borderSize = ceil(borderSize ./ rs) .* rs, or account for the
            % offset when computing the output grid.
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            XR_resample_dataset(dataPath, [2, 2, 2], 'channelPatterns', testCase.Channels(1), 'zarrFile', true, ...
                'largeFile', true, 'saveZarr', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 2], ...
                'parseCluster', false, 'masterCompute', true);
            out = readzarr([dataPath, 'resampled/', names{1}, '.zarr']);
            expected = uint16(imresize3(single(ims{1}), round(testCase.ImSize / 2), 'Method', 'linear'));
            testCase.verifyImagesClose(out, expected, 'relTol', 1e-4, 'absTol', 1);
        end
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function image2dWithYxzFactor(testCase)
            % KNOWN ISSUE: for 2-D images XR_resampleFrame.m expands the factor
            % with [ones(1, 3 - numel(rs)) * rs(1), rs(2 : end)], so a yxz
            % factor [a, b, c] becomes [b, c] instead of [a, b]. Dataset
            % wrappers pass yxz factors, so 2-D images are resampled with the
            % x and z factors. Fix: use rs(1 : 2) for 2-D images.
            dataPath = testCase.createTempFolder();
            writetiff(uint16(randi(60000, 40, 36)), [dataPath, 'im.tif']);
            XR_resampleFrame([dataPath, 'im.tif'], [dataPath, 'out.tif'], [1, 2, 4]);
            testCase.verifyEqual(size(readtiff([dataPath, 'out.tif'])), [40, 18], ...
                'A 2-D image was not resampled with the y and x factors.');
        end
    end

    methods (Access = private)
        function [im, fn, dataPath] = writeImage(testCase, format)
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            if strcmp(format, 'zarr')
                fn = [dataPath, 'im.zarr'];
                writezarr(im, fn, 'blockSize', [16, 16, 4]);
            else
                fn = [dataPath, 'im.tif'];
                writetiff(im, fn);
            end
        end

        function [dataPath, ims, names] = writeDataset(testCase, format)
            dataPath = testCase.createTempFolder();
            names = {'scan_CamA_ch0_t0000', 'scan_CamB_ch0_t0000', 'scan_CamA_ch1_t0000'};
            ims = cell(1, 3);
            for i = 1 : 3
                ims{i} = uint16(randi(60000, testCase.ImSize));
                if strcmp(format, 'zarr')
                    writezarr(ims{i}, [dataPath, names{i}, '.zarr'], 'blockSize', [16, 16, 4]);
                else
                    writetiff(ims{i}, [dataPath, names{i}, '.tif']);
                end
            end
        end
    end

    methods (Static, Access = private)
        function [bmax, bmean] = blockStats(im, rs)
            sz = size(im, 1 : 3);
            outSize = ceil(sz ./ rs);
            bmax = zeros(outSize, 'single');
            bmean = zeros(outSize, 'single');
            for k = 1 : outSize(3)
                for j = 1 : outSize(2)
                    for i = 1 : outSize(1)
                        block = im((i - 1) * rs(1) + 1 : min(i * rs(1), sz(1)), ...
                            (j - 1) * rs(2) + 1 : min(j * rs(2), sz(2)), (k - 1) * rs(3) + 1 : min(k * rs(3), sz(3)));
                        bmax(i, j, k) = max(block, [], 'all');
                        bmean(i, j, k) = mean(block, 'all');
                    end
                end
            end
        end
    end
end
