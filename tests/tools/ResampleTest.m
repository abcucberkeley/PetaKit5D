classdef ResampleTest < pktest.TestCase
    % XR_resample_dataset / XR_resampleFrame / XR_resampleSingleZarr.
    % resampleFactor is in yxz order; a factor of 2 downsamples 2x.

    properties (Constant)
        ImSize = [40, 36, 12]
        Channels = {'CamA_ch0', 'CamB_ch0'}
    end

    properties (TestParameter)
        % 'down3' doesn't divide the image size (40 / 3 = 13.3), which pins
        % the rounding of the output size
        resampleFactor = struct('down2', [2, 2, 2], 'zOnly', [1, 1, 2], 'up2xy', [0.5, 0.5, 1], ...
            'nonInteger', [1.5, 1.5, 1.5], 'down3', [3, 3, 3])
        interpMethod = {'linear', 'cubic', 'nearest', 'max', 'mean'}
        % imresize3_fast: dimensions that keep their size, up and down,
        % non-integer, and all orders of the passes
        fastFactor = struct('zOnly', [1, 1, 2], 'xyOnly', [2, 2, 1], 'yOnly', [3, 1, 1], 'zUp', [1, 1, 0.5], ...
            'all', [2, 2, 2], 'nonInteger', [1.5, 1.3, 2.5], 'zBetween', [2, 4, 3], 'xzDownYUp', [0.7, 1.6, 3], ...
            'none', [1, 1, 1])
        fastMethod = {'linear', 'cubic'}
        fastType = {'single', 'uint16'}
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

        function imresize3FastMatchesImresize3(testCase, fastFactor, fastMethod, fastType)
            % imresize3_fast replaces imresize3 in the resampling; imresize3
            % is replaced by a spy that errors when called with an output, so
            % the result comes from the fast path, not from the fallback
            im = cast(rand(41, 37, 23) * 60000, fastType);
            outSize = round(size(im) ./ fastFactor);
            expected = imresize3(im, outSize, fastMethod);
            testCase.applyFixture(pktest.FunctionSpyFixture('imresize3'));
            testCase.verifyEqual(imresize3_fast(im, outSize, fastMethod), expected);
        end

        function image2dWithYxzFactor(testCase)
            % 2D images used to be resampled with the x and z factors of a yxz factor
            dataPath = testCase.createTempFolder();
            writetiff(uint16(randi(60000, 40, 36)), [dataPath, 'im.tif']);
            XR_resampleFrame([dataPath, 'im.tif'], [dataPath, 'out.tif'], [1, 2, 4]);
            testCase.verifyEqual(size(readtiff([dataPath, 'out.tif'])), [40, 18], ...
                'A 2-D image was not resampled with the y and x factors.');
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

        function largeZarrKeepsResultOfAnotherCopy(testCase)
            % two copies of a task (e.g. a retry next to a slow first try)
            % share the temporary result; the copy that finished second
            % deleted the result the first one had moved into place
            [~, fn, dataPath] = testCase.writeImage('zarr');
            outFn = [dataPath, 'resampled.zarr'];
            % the other copy moves the result into place while this one runs
            testCase.applyFixture(pktest.FrameworkSpyFixture(@() movefile([outFn, '_test'], outFn)));
            XR_resampleSingleZarr(fn, outFn, [2, 2, 2], 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 2], ...
                'parseCluster', false, 'uuid', 'test');
            testCase.verifyTrue(isfolder(outFn), 'The result of the other copy was deleted.');
        end

        function largeZarrKeepsResultOfAnotherRun(testCase)
            % a copy from another run (its own uuid) that finds the result
            % published replaced it with its own, incomplete result
            [im, fn, dataPath] = testCase.writeImage('zarr');
            outFn = [dataPath, 'resampled.zarr'];
            published = im(1 : 2 : end, 1 : 2 : end, 1 : 2 : end);
            testCase.applyFixture(pktest.FrameworkSpyFixture(@() writezarr(published, outFn)));
            XR_resampleSingleZarr(fn, outFn, [2, 2, 2], 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 2], ...
                'parseCluster', false, 'uuid', 'test');
            testCase.verifyEqual(readzarr(outFn), published, 'The published result was replaced.');
            testCase.verifyFalse(isfolder([outFn, '_test']), 'The temporary result was left behind.');
        end

        function largeZarrWithDefaultBorderSize(testCase)
            % batches used to be shifted by half an output voxel when the border
            % (default [5, 5, 5]) is not a multiple of the factor
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            XR_resample_dataset(dataPath, [2, 2, 2], 'channelPatterns', testCase.Channels(1), 'zarrFile', true, ...
                'largeFile', true, 'saveZarr', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 2], ...
                'parseCluster', false, 'masterCompute', true);
            out = readzarr([dataPath, 'resampled/', names{1}, '.zarr']);
            expected = uint16(imresize3(single(ims{1}), round(testCase.ImSize / 2), 'Method', 'linear'));
            testCase.verifyImagesClose(out, expected, 'relTol', 1e-4, 'absTol', 1);
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
