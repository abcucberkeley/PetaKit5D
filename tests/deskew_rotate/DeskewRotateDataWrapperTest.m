classdef DeskewRotateDataWrapperTest < pktest.TestCase
    % XR_deskew_rotate_data_wrapper on small synthetic datasets, run locally
    % (parseCluster=false). Outputs must equal the in-memory functions.

    properties (Constant)
        SkewAngle = 32.45
        XYPixelSize = 0.108
        Dz = 0.2
        RawSize = [64, 96, 80]
        FsnA0 = 'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
        FsnB0 = 'Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
        FsnA1 = 'Scan_Iter_0001_0000_CamA_ch0_CAM1_stack0000_488nm_0011000msec_0106071251msecAbs_000x_000y_000z_0001t'
    end

    methods (Test, TestTags = {'Integration'})
        function processesEveryMatchingFile(testCase)
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath);

            fsns = {testCase.FsnA0, testCase.FsnB0, testCase.FsnA1};
            for i = 1 : numel(fsns)
                dsr = readtiff([dataPath, 'DSR/', fsns{i}, '.tif']);
                testCase.verifyEqual(dsr, testCase.expectedDSR(frames{i}), fsns{i});
                testCase.verifyEqual(readtiff([dataPath, 'DSR/MIPs/', fsns{i}, '_MIP_z.tif']), max(dsr, [], 3));
            end
            testCase.verifyTrue(isfile([dataPath, 'DSR/parameters.json']));
            testCase.verifyTrue(isfile([dataPath, 'DSR/parameters.mat']));
            testCase.verifyEmpty(dir([dataPath, 'DSR/*.tmp']), 'Task flag files were left behind.');
            pr = jsondecode(fileread([dataPath, 'DSR/parameters.json']));
            testCase.verifyEqual(pr.dz, testCase.Dz);
            testCase.verifyEqual(pr.skewAngle, testCase.SkewAngle);
        end

        function channelPatternsSelectFiles(testCase)
            dataPath = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'channelPatterns', {'CamB'});
            dir_info = dir([dataPath, 'DSR/*.tif']);
            testCase.verifyEqual({dir_info.name}, {[testCase.FsnB0, '.tif']});
        end

        function multipleDataPathsAreProcessed(testCase)
            rootPath = testCase.createTempFolder();
            dataPaths = {[rootPath, 'tile1/'], [rootPath, 'tile2/']};
            frame = testCase.makeFrame(1);
            for i = 1 : 2
                writetiff(frame, [dataPaths{i}, testCase.FsnA0, '.tif']);
            end
            testCase.runWrapper(dataPaths);
            for i = 1 : 2
                testCase.verifyEqual(readtiff([dataPaths{i}, 'DSR/', testCase.FsnA0, '.tif']), testCase.expectedDSR(frame));
            end
        end

        function deskewOnlyWritesDS(testCase)
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'rotate', false, 'DSRCombined', false);
            ds = readtiff([dataPath, 'DS/', testCase.FsnA0, '.tif']);
            testCase.verifyEqual(ds, uint16(deskewFrame3D(single(frames{1}), testCase.SkewAngle, testCase.Dz, ...
                testCase.XYPixelSize, true, 'Crop', false)));
            testCase.verifyFalse(isfolder([dataPath, 'DSR']));
        end

        function deskewOnlyWithCombinedDefaultWritesDS(testCase)
            % DSRCombined=true (the default) used to send deskew only to the
            % combined path, which writes DSR, and then wait minutes for DS
            % before failing
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'rotate', false, 'DSRCombined', true);
            ds = readtiff([dataPath, 'DS/', testCase.FsnA0, '.tif']);
            testCase.verifyEqual(ds, uint16(deskewFrame3D(single(frames{1}), testCase.SkewAngle, testCase.Dz, ...
                testCase.XYPixelSize, true, 'Crop', false)));
            testCase.verifyFalse(isfolder([dataPath, 'DSR']));
        end

        function separateDeskewRotateWritesBoth(testCase)
            dataPath = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'DSRCombined', false);
            testCase.verifyTrue(isfile([dataPath, 'DS/', testCase.FsnA0, '.tif']));
            dsrSize = getImageSize([dataPath, 'DSR/', testCase.FsnA0, '.tif']);
            testCase.verifyEqual(dsrSize, size(testCase.expectedDSR(testCase.makeFrame(1))));
        end

        function existingOutputsAreSkippedUnlessOverwrite(testCase)
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath);
            dsrFn = [dataPath, 'DSR/', testCase.FsnA0, '.tif'];

            writetiff(frames{1} * 2, [dataPath, testCase.FsnA0, '.tif']);
            testCase.runWrapper(dataPath);
            testCase.verifyEqual(readtiff(dsrFn), testCase.expectedDSR(frames{1}));

            testCase.runWrapper(dataPath, 'overwrite', true);
            testCase.verifyEqual(readtiff(dsrFn), testCase.expectedDSR(frames{1} * 2));
        end

        function saveZarrMatchesTiff(testCase)
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'saveZarr', true);
            testCase.verifyEqual(readzarr([dataPath, 'DSR/', testCase.FsnA0, '.zarr']), testCase.expectedDSR(frames{1}));
        end

        function zarrInputMatchesTiffInput(testCase)
            dataPath = testCase.createTempFolder();
            frame = testCase.makeFrame(1);
            writezarr(frame, [dataPath, testCase.FsnA0, '.zarr']);
            testCase.runWrapper(dataPath, 'zarrFile', true);
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.FsnA0, '.tif']), testCase.expectedDSR(frame));
        end

        function largeFileMatchesRegular(testCase)
            % the large file path splits the output along y into batches
            % (4 here) and must reproduce the in-memory result exactly
            dataPath = testCase.createTempFolder();
            frame = testCase.makeFrame(1);
            writezarr(frame, [dataPath, testCase.FsnA0, '.zarr']);
            testCase.runWrapper(dataPath, 'zarrFile', true, 'largeFile', true, 'batchSize', [16, 1024, 1024], ...
                'blockSize', [16, 64, 64]);
            dsr = readzarr([dataPath, 'DSR/', testCase.FsnA0, '.zarr']);
            testCase.verifyEqual(dsr, testCase.expectedDSR(frame));
            testCase.verifyEqual(readtiff([dataPath, 'DSR/MIPs/', testCase.FsnA0, '_MIP_z.tif']), max(dsr, [], 3));
        end

        function flatFieldUsesEachChannelsImages(testCase)
            [dataPath, frames] = testCase.writeDataset();
            bg = zeros(testCase.RawSize(1 : 2), 'single') + 100;
            ffA = bg + single(linspace(500, 1000, testCase.RawSize(2)));
            ffB = bg + single(linspace(1000, 500, testCase.RawSize(2)));
            writetiff(ffA, [dataPath, 'ff/ffA.tif']);
            writetiff(ffB, [dataPath, 'ff/ffB.tif']);
            writetiff(bg, [dataPath, 'ff/bg.tif']);
            testCase.runWrapper(dataPath, 'FFCorrection', true, 'lowerLimit', 0.4, ...
                'FFImagePaths', {[dataPath, 'ff/ffA.tif'], [dataPath, 'ff/ffB.tif']}, ...
                'backgroundPaths', {[dataPath, 'ff/bg.tif'], [dataPath, 'ff/bg.tif']});

            expectedA = testCase.expectedDSR(XR_LSFlatFieldCorrection(frames{1}, ffA, bg, 'LowerLimit', 0.4));
            expectedB = testCase.expectedDSR(XR_LSFlatFieldCorrection(frames{2}, ffB, bg, 'LowerLimit', 0.4));
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.FsnA0, '.tif']), expectedA);
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.FsnB0, '.tif']), expectedB);
        end

        function partialFilesAreJoinedAlongZ(testCase)
            % acquisitions split into <name>.tif + <name>_part0001.tif are
            % processed as one volume
            dataPath = testCase.createTempFolder();
            frame = testCase.makeFrame(1);
            writetiff(frame(:, :, 1 : 40), [dataPath, testCase.FsnA0, '.tif']);
            writetiff(frame(:, :, 41 : end), [dataPath, testCase.FsnA0, '_part0001.tif']);
            testCase.runWrapper(dataPath);
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.FsnA0, '.tif']), testCase.expectedDSR(frame));
            testCase.verifyFalse(isfile([dataPath, 'DSR/', testCase.FsnA0, '_part0001.tif']));
        end

        function saveMIPFalseWritesNoMIP(testCase)
            % saveMIP and blockSize used to be dropped from the task string
            % of regular (not large) files
            dataPath = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'saveMIP', false);
            testCase.verifyFalse(isfolder([dataPath, 'DSR/MIPs']), 'MIPs were written although saveMIP is false.');
        end

        function blockSizeIsUsedForZarrOutput(testCase)
            dataPath = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'saveZarr', true, 'blockSize', [16, 32, 16]);
            info = getZarrInfo([dataPath, 'DSR/', testCase.FsnA0, '.zarr']);
            testCase.verifyEqual(info.blockSize(:)', [16, 32, 16]);
        end

        function objectiveScanLargeFileIsRejected(testCase)
            % the large-file block processing only supports sample scan; it
            % used to fail inside the tasks with an undefined variable
            dataPath = testCase.createTempFolder();
            writezarr(testCase.makeFrame(1), [dataPath, testCase.FsnA0, '.zarr']);
            try
                testCase.runWrapper(dataPath, 'zarrFile', true, 'largeFile', true, 'objectiveScan', true, ...
                    'batchSize', [16, 1024, 1024]);
                testCase.verifyFail('Large-file objective scan did not raise an error.');
            catch ME
                testCase.verifySubstring(ME.message, 'objective scan is not supported');
            end
        end
    end

    methods (Access = private)
        function [] = runWrapper(testCase, dataPaths, varargin)
            if ischar(dataPaths)
                dataPaths = {dataPaths};
            end
            XR_deskew_rotate_data_wrapper(dataPaths, 'deskew', true, 'rotate', true, 'DSRCombined', true, ...
                'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, 'skewAngle', testCase.SkewAngle, ...
                'reverse', true, 'channelPatterns', {'CamA', 'CamB'}, 'save16bit', true, ...
                'parseCluster', false, 'masterCompute', true, varargin{:});
        end

        function [frame] = makeFrame(testCase, seed)
            % uint16 beads on a background of 100, different per seed
            rs = RandStream('twister', 'Seed', seed);
            raw = [8 + 48 * rand(rs, 15, 1), 26 + 44 * rand(rs, 15, 1), 10 + 60 * rand(rs, 15, 1)];
            frame = uint16(pktest.beads(testCase.RawSize, raw, 2, 1000, 100));
        end

        function [dataPath, frames] = writeDataset(testCase)
            dataPath = testCase.createTempFolder();
            fsns = {testCase.FsnA0, testCase.FsnB0, testCase.FsnA1};
            frames = cell(1, numel(fsns));
            for i = 1 : numel(fsns)
                frames{i} = testCase.makeFrame(i);
                writetiff(frames{i}, [dataPath, fsns{i}, '.tif']);
            end
        end

        function [dsr] = expectedDSR(testCase, frame)
            dsr = deskewRotateFrame3D(frame, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                'reverse', true, 'save16bit', true);
        end
    end
end
