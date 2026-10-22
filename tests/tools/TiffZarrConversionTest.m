classdef TiffZarrConversionTest < pktest.TestCase
    % tiff <-> zarr conversion (XR_tiffToZarr_wrapper, XR_zarrToTiff_wrapper,
    % tiffToZarr, zarrToTiff) and zarr resaving (XR_resave_zarr_wrapper).
    % XR_zarrToTiff_wrapper doesn't forward parseCluster to the job
    % framework, so on a machine with Slurm it submits jobs; the shared
    % pktest.NoSlurmFixture keeps it local.

    properties (Constant)
        ImSize = [40, 36, 12]
        Channels = {'CamA_ch0', 'CamB_ch0'}
    end

    methods (Test, TestTags = {'Integration'})
        function tiffToZarrWrapper(testCase)
            [dataPath, ims, names] = testCase.writeDataset('tif');
            XR_tiffToZarr_wrapper(dataPath, 'channelPatterns', testCase.Channels, 'blockSize', [16, 16, 4], ...
                'parseCluster', false, 'masterCompute', true);
            resultPath = [dataPath, 'zarr/'];
            for i = 1 : 2
                fn = [resultPath, names{i}, '.zarr'];
                testCase.verifyEqual(readzarr(fn), ims{i}, names{i});
                info = getZarrInfo(fn);
                testCase.verifyEqual(info.blockSize, [16, 16, 4]);
                testCase.verifyEqual(info.compressor.cname, 'zstd');
            end
            testCase.verifyFalse(exist([resultPath, names{3}, '.zarr'], 'dir') > 0, ...
                'A file that matches no channel pattern was converted.');
        end

        function tiffToZarrWrapperWithInputBbox(testCase)
            [dataPath, ims, names] = testCase.writeDataset('tif');
            b = [5, 7, 3, 30, 33, 10];
            XR_tiffToZarr_wrapper(dataPath, 'channelPatterns', testCase.Channels(1), 'inputBbox', b, ...
                'resultDirName', 'tiff_to_zarr', 'parseCluster', false, 'masterCompute', true);
            out = readzarr([dataPath, 'tiff_to_zarr/', names{1}, '.zarr']);
            testCase.verifyEqual(out, ims{1}(b(1) : b(4), b(2) : b(5), b(3) : b(6)));
        end

        function tiffToZarrWrapperWithoutReadingWholeTiff(testCase)
            % reads the tiff in blocks with blockedImage instead of at once
            [dataPath, ims, names] = testCase.writeDataset('tif');
            b = [1, 1, 1, testCase.ImSize];
            XR_tiffToZarr_wrapper(dataPath, 'channelPatterns', testCase.Channels(1), 'readWholeTiff', false, ...
                'inputBbox', b, 'blockSize', [16, 16, 4], 'parseCluster', false, 'masterCompute', true);
            testCase.verifyEqual(readzarr([dataPath, 'zarr/', names{1}, '.zarr']), ims{1});
        end

        function zarrToTiffWrapper(testCase)
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            XR_zarrToTiff_wrapper(dataPath, 'channelPatterns', testCase.Channels, ...
                'parseCluster', false, 'masterCompute', true);
            resultPath = [dataPath, 'tiffs/'];
            for i = 1 : 2
                testCase.verifyEqual(readtiff([resultPath, names{i}, '.tif']), ims{i}, names{i});
            end
            testCase.verifyFalse(exist([resultPath, names{3}, '.tif'], 'file') > 0, ...
                'A file that matches no channel pattern was converted.');
        end

        function zarrToTiffWrapperPassesJobSettings(testCase)
            % each task converts one file, so the memory per task is for one
            % file (it was the sum over all files), and parseCluster,
            % masterCompute and cpusPerTask were not passed on
            dataPath = testCase.writeDataset('zarr');
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            XR_zarrToTiff_wrapper(dataPath, 'channelPatterns', testCase.Channels, ...
                'parseCluster', false, 'masterCompute', false, 'cpusPerTask', 3);
            testCase.verifyEqual(spy.getOption(1, 'parseCluster'), false);
            testCase.verifyEqual(spy.getOption(1, 'masterCompute'), false);
            testCase.verifyEqual(spy.getOption(1, 'cpusPerTask'), 3);
            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), prod(testCase.ImSize) * 4 / 1024 ^ 3 * 2.5, ...
                'RelTol', 1e-12);
        end

        function resaveZarrChangesChunks(testCase)
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            XR_resave_zarr_wrapper(dataPath, 'blockSize', [8, 8, 4], 'channelPatterns', {'.zarr'}, ...
                'largeFile', false, 'parseCluster', false, 'masterCompute', true);
            for i = 1 : 3
                fn = [dataPath, 'resaved_zarr/', names{i}, '.zarr'];
                testCase.verifyEqual(readzarr(fn), ims{i}, names{i});
                testCase.verifyEqual(getZarrInfo(fn).blockSize, [8, 8, 4]);
            end
        end

        function resaveKeepsResultOfAnotherCopy(testCase)
            % two copies of a task (e.g. a retry next to a slow first try)
            % share the temporary result; the copy that finished second
            % deleted the result the first one had moved into place
            [dataPath, ~, names] = testCase.writeDataset('zarr');
            outFn = [dataPath, 'resaved.zarr'];
            % the other copy moves the result into place while this one runs
            testCase.applyFixture(pktest.FrameworkSpyFixture(@() movefile([outFn, '_test'], outFn)));
            XR_resaveSingleZarr([dataPath, names{1}, '.zarr'], outFn, 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'uuid', 'test');
            testCase.verifyTrue(isfolder(outFn), 'The result of the other copy was deleted.');
        end

        function resaveLargeZarrInBatches(testCase)
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            b = [5, 7, 3, 30, 33, 10];
            XR_resave_zarr_wrapper(dataPath, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], 'inputBbox', b, ...
                'channelPatterns', testCase.Channels(1), 'largeFile', true, 'parseCluster', false, 'masterCompute', true);
            out = readzarr([dataPath, 'resaved_zarr/', names{1}, '.zarr']);
            testCase.verifyEqual(out, ims{1}(b(1) : b(4), b(2) : b(5), b(3) : b(6)));
        end
    end

    methods (Test, TestTags = {'Unit'})
        function singleFileRoundTrip(testCase)
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            writetiff(im, [dataPath, 'im.tif']);
            % the third input (in-memory frame) must be given, see the known
            % issue optionalFrameCanBeOmitted
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'im.zarr'], [], 'blockSize', [16, 16, 4]);
            testCase.verifyEqual(readzarr([dataPath, 'im.zarr']), im);
            testCase.verifyEqual(getZarrInfo([dataPath, 'im.zarr']).blockSize, [16, 16, 4]);
            zarrToTiff([dataPath, 'im.zarr'], [dataPath, 'back.tif']);
            testCase.verifyEqual(readtiff([dataPath, 'back.tif']), im);
        end

        function tiffToZarrFromMemory(testCase)
            % the optional third input converts data that is already loaded
            dataPath = testCase.createTempFolder();
            im = single(randn(testCase.ImSize));
            tiffToZarr('', [dataPath, 'im.zarr'], im, 'blockSize', [16, 16, 4]);
            testCase.verifyEqual(readzarr([dataPath, 'im.zarr']), im);
        end

        function tiffToZarrProcessingOptions(testCase)
            % flip in z, then crop (tileOutBbox), resample, and the user
            % function, in that order; the result has the input's type
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            writetiff(im, [dataPath, 'im.tif']);
            b = [3, 4, 2, 38, 30, 11];
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'flip.zarr'], [], 'blockSize', [16, 16, 4], ...
                'flipZstack', true, 'tileOutBbox', b, 'usrFcn', '@(x) single(x) / 3');
            expected = flip(im, 3);
            expected = uint16(single(expected(b(1) : b(4), b(2) : b(5), b(3) : b(6))) / 3);
            testCase.verifyEqual(readzarr([dataPath, 'flip.zarr']), expected);

            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'rs.zarr'], [], 'blockSize', [16, 16, 4], ...
                'resampleFactor', [1, 1, 2]);
            testCase.verifyEqual(readzarr([dataPath, 'rs.zarr']), imresize3(im, round(testCase.ImSize ./ [1, 1, 2])));
        end

        function tiffToZarrReadsRegions(testCase)
            % readWholeTiff false reads the tiff by regions; without
            % inputBbox it failed on an undefined variable (outBbox)
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            writetiff(im, [dataPath, 'im.tif']);
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'whole.zarr'], [], 'blockSize', [16, 16, 4], ...
                'readWholeTiff', false);
            testCase.verifyEqual(readzarr([dataPath, 'whole.zarr']), im);
            b = [3, 4, 2, 38, 30, 11];
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'crop.zarr'], [], 'blockSize', [16, 16, 4], ...
                'readWholeTiff', false, 'tileOutBbox', b);
            testCase.verifyEqual(readzarr([dataPath, 'crop.zarr']), im(b(1) : b(4), b(2) : b(5), b(3) : b(6)));
        end

        function tiffToZarrPartialFilesAndFrames(testCase)
            % partial files are joined along z; an in-memory frame gets the
            % same processing as a file; a 2D image becomes one z slice
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            writetiff(im(:, :, 1 : 5), [dataPath, 'part1.tif']);
            writetiff(im(:, :, 6 : end), [dataPath, 'part2.tif']);
            tiffToZarr({[dataPath, 'part1.tif'], [dataPath, 'part2.tif']}, [dataPath, 'parts.zarr'], [], ...
                'blockSize', [16, 16, 4]);
            testCase.verifyEqual(readzarr([dataPath, 'parts.zarr']), im);

            tiffToZarr('', [dataPath, 'frame.zarr'], im, 'blockSize', [16, 16, 4], 'flipZstack', true);
            testCase.verifyEqual(readzarr([dataPath, 'frame.zarr']), flip(im, 3));

            writetiff(im(:, :, 1), [dataPath, 'im2d.tif']);
            tiffToZarr([dataPath, 'im2d.tif'], [dataPath, 'im2d.zarr'], [], 'blockSize', [16, 16, 4]);
            testCase.verifyEqual(readzarr([dataPath, 'im2d.zarr']), im(:, :, 1));
            testCase.verifyEqual(getZarrInfo([dataPath, 'im2d.zarr']).size, [testCase.ImSize(1 : 2), 1]);
        end
    end

    methods (Access = private)
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
end
