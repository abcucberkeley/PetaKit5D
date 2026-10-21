classdef CropDatasetTest < pktest.TestCase
    % XR_crop_dataset / XR_crop_frame / XR_crop_zarr on small synthetic
    % datasets (two channels plus a channel that must be ignored)

    properties (Constant)
        ImSize = [40, 36, 12]
        Bbox = [5, 7, 3, 30, 33, 10]
        Channels = {'CamA_ch0', 'CamB_ch0'}
    end

    properties
        DataPath
        Images
        Names
    end

    methods (Test, TestTags = {'Integration'})
        function tiffToTiff(testCase)
            testCase.writeDataset('tif');
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'parseCluster', false, 'masterCompute', true);

            resultPath = [testCase.DataPath, 'Cropped/'];
            for i = 1 : 2
                out = readtiff([resultPath, testCase.Names{i}, '.tif']);
                testCase.verifyEqual(out, testCase.crop(testCase.Images{i}, testCase.Bbox), testCase.Names{i});
            end
            testCase.verifyFalse(exist([resultPath, testCase.Names{3}, '.tif'], 'file') > 0, ...
                'A file that matches no channel pattern was cropped.');
            testCase.verifyTrue(exist([resultPath, 'parameters.mat'], 'file') > 0);
        end

        function zarrToZarr(testCase)
            testCase.writeDataset('zarr');
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'zarrFile', true, 'saveZarr', true, 'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                out = readzarr([testCase.DataPath, 'Cropped/', testCase.Names{i}, '.zarr']);
                testCase.verifyEqual(out, testCase.crop(testCase.Images{i}, testCase.Bbox), testCase.Names{i});
            end
        end

        function zarrToTiff(testCase)
            testCase.writeDataset('zarr');
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'zarrFile', true, 'saveZarr', false, 'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                out = readtiff([testCase.DataPath, 'Cropped/', testCase.Names{i}, '.tif']);
                testCase.verifyEqual(out, testCase.crop(testCase.Images{i}, testCase.Bbox), testCase.Names{i});
            end
        end

        function largeZarrMatchesInMemoryCrop(testCase)
            % small batches so the crop is assembled from several blocks
            testCase.writeDataset('zarr');
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'zarrFile', true, 'largeFile', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                out = readzarr([testCase.DataPath, 'Cropped/', testCase.Names{i}, '.zarr']);
                testCase.verifyEqual(out, testCase.crop(testCase.Images{i}, testCase.Bbox), testCase.Names{i});
            end
        end

        function padPastImageEdge(testCase)
            % bbox starts before and ends after the image in every axis
            testCase.writeDataset('tif');
            bbox = [-2, -1, 0, 45, 40, 15];
            XR_crop_dataset(testCase.DataPath, bbox, 'channelPatterns', testCase.Channels, 'pad', true, ...
                'parseCluster', false, 'masterCompute', true);
            out = readtiff([testCase.DataPath, 'Cropped/', testCase.Names{1}, '.tif']);
            testCase.verifyEqual(out, testCase.padCrop(testCase.Images{1}, bbox));
        end

        function padPastImageEdgeLargeZarr(testCase)
            % regression test for padding in large zarr cropping (81f3341)
            testCase.writeDataset('zarr');
            bbox = [-2, -1, 0, 45, 40, 15];
            XR_crop_dataset(testCase.DataPath, bbox, 'channelPatterns', testCase.Channels, 'pad', true, ...
                'zarrFile', true, 'largeFile', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'masterCompute', true);
            out = readzarr([testCase.DataPath, 'Cropped/', testCase.Names{1}, '.zarr']);
            testCase.verifyEqual(out, testCase.padCrop(testCase.Images{1}, bbox));
        end

        function clipPastImageEdgeWithoutPad(testCase)
            testCase.writeDataset('tif');
            bbox = [30, 20, 8, 60, 50, 20];
            XR_crop_dataset(testCase.DataPath, bbox, 'channelPatterns', testCase.Channels, 'pad', false, ...
                'parseCluster', false, 'masterCompute', true);
            out = readtiff([testCase.DataPath, 'Cropped/', testCase.Names{1}, '.tif']);
            testCase.verifyEqual(out, testCase.Images{1}(30 : end, 20 : end, 8 : end));
        end

        function singleSliceZarrCrop(testCase)
            % regression test for single slice crops of zarr files (d513b73)
            testCase.writeDataset('zarr');
            bbox = [5, 7, 6, 30, 33, 6];
            XR_crop_dataset(testCase.DataPath, bbox, 'channelPatterns', testCase.Channels(1), ...
                'zarrFile', true, 'saveZarr', true, 'parseCluster', false, 'masterCompute', true, ...
                'resultDirName', 'Cropped_small');
            XR_crop_dataset(testCase.DataPath, bbox, 'channelPatterns', testCase.Channels(1), ...
                'zarrFile', true, 'largeFile', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'masterCompute', true, 'resultDirName', 'Cropped_large');
            expected = testCase.crop(testCase.Images{1}, bbox);
            testCase.verifyEqual(readzarr([testCase.DataPath, 'Cropped_small/', testCase.Names{1}, '.zarr']), expected);
            testCase.verifyEqual(readzarr([testCase.DataPath, 'Cropped_large/', testCase.Names{1}, '.zarr']), expected);
        end

        function centerCrop(testCase)
            testCase.writeDataset('tif');
            XR_crop_dataset(testCase.DataPath, [1, 1, 1, 10, 9, 4], 'cropType', 'center', ...
                'channelPatterns', testCase.Channels(1), 'parseCluster', false, 'masterCompute', true);
            % start = round((imSize + 1) / 2) - floor(bboxSize / 2)
            start = round((testCase.ImSize + 1) / 2) - floor([10, 9, 4] / 2);
            out = readtiff([testCase.DataPath, 'Cropped/', testCase.Names{1}, '.tif']);
            testCase.verifyEqual(out, testCase.crop(testCase.Images{1}, [start, start + [10, 9, 4] - 1]));
        end

        function movingCrop(testCase)
            % the bbox moves linearly from inputBbox to lastStartCoords over
            % the time points of a channel
            dataPath = testCase.createTempFolder();
            nt = 3;
            ims = cell(nt, 1);
            for t = 1 : nt
                ims{t} = uint16(randi(60000, testCase.ImSize));
                writetiff(ims{t}, sprintf('%sscan_CamA_ch0_t%04d.tif', dataPath, t - 1));
            end
            bbox = [3, 3, 2, 12, 12, 6];
            lastStart = [21, 13, 6];
            XR_crop_dataset(dataPath, bbox, 'cropType', 'moving', 'lastStartCoords', lastStart, ...
                'channelPatterns', testCase.Channels(1), 'parseCluster', false, 'masterCompute', true);
            bboxSize = bbox(4 : 6) - bbox(1 : 3) + 1;
            for t = 1 : nt
                s = bbox(1 : 3) + round((lastStart - bbox(1 : 3)) / (nt - 1) * (t - 1));
                out = readtiff(sprintf('%sCropped/scan_CamA_ch0_t%04d.tif', dataPath, t - 1));
                testCase.verifyEqual(out, testCase.crop(ims{t}, [s, s + bboxSize - 1]), sprintf('time point %d', t));
            end
        end

        function existingOutputIsNotOverwritten(testCase)
            testCase.writeDataset('tif');
            resultPath = [testCase.DataPath, 'Cropped/'];
            mkdir(resultPath);
            sentinel = uint16(magic(4));
            writetiff(sentinel, [resultPath, testCase.Names{1}, '.tif']);
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'parseCluster', false, 'masterCompute', true);
            testCase.verifyEqual(readtiff([resultPath, testCase.Names{1}, '.tif']), sentinel);
            testCase.verifyEqual(readtiff([resultPath, testCase.Names{2}, '.tif']), testCase.crop(testCase.Images{2}, testCase.Bbox));
        end

        function tiffInputSavedAsZarr(testCase)
            % tif input saved as zarr used to get a zarr folder named *.tif
            testCase.writeDataset('tif');
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels(1), ...
                'saveZarr', true, 'parseCluster', false, 'masterCompute', true);
            fn = [testCase.DataPath, 'Cropped/', testCase.Names{1}, '.zarr'];
            testCase.verifyTrue(exist(fn, 'dir') > 0, 'Cropped zarr output was not saved as <name>.zarr.');
            if exist(fn, 'dir')
                testCase.verifyEqual(readzarr(fn), testCase.crop(testCase.Images{1}, testCase.Bbox));
            end
        end
    end

    methods (Access = private)
        function [] = writeDataset(testCase, format)
            testCase.DataPath = testCase.createTempFolder();
            testCase.Names = {'scan_CamA_ch0_t0000', 'scan_CamB_ch0_t0000', 'scan_CamA_ch1_t0000'};
            testCase.Images = cell(1, 3);
            for i = 1 : 3
                testCase.Images{i} = uint16(randi(60000, testCase.ImSize));
                fn = [testCase.DataPath, testCase.Names{i}];
                if strcmp(format, 'zarr')
                    writezarr(testCase.Images{i}, [fn, '.zarr'], 'blockSize', [16, 16, 4]);
                else
                    writetiff(testCase.Images{i}, [fn, '.tif']);
                end
            end
        end
    end

    methods (Static, Access = private)
        function [out] = crop(im, b)
            out = im(b(1) : b(4), b(2) : b(5), b(3) : b(6));
        end

        function [out] = padCrop(im, b)
            % crop with zero padding outside the image
            out = zeros(b(4 : 6) - b(1 : 3) + 1, class(im));
            s = max(b(1 : 3), 1);
            t = min(b(4 : 6), size(im));
            o = s - b(1 : 3) + 1;
            out(o(1) : o(1) + t(1) - s(1), o(2) : o(2) + t(2) - s(2), o(3) : o(3) + t(3) - s(3)) = ...
                im(s(1) : t(1), s(2) : t(2), s(3) : t(3));
        end
    end
end
