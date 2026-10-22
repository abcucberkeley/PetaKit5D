classdef DeskewRotateMemoryTest < pktest.TestCase
    % memory requested for deskew/rotate tasks: XR_estimateDeskewRotateMemory
    % (the arrays a code path holds at once plus MATLAB itself) and the
    % memAllocate the wrappers pass to the job framework. The wrappers used
    % to request 40 bytes per input voxel (about 4x the peak of the default
    % path, so large tiles took several CPUs for nothing) and counted only
    % the first of a group of partial files.

    properties (Constant)
        SkewAngle = 32.45
        XYPixelSize = 0.108
        RawSize = [64, 96, 80]
        Fsn = 'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
        FsnB = 'Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
    end

    properties (TestParameter)
        % options of deskewRotateFrame3D / XR_estimateDeskewRotateMemory
        % that change the output size
        geometry = struct( ...
            'dz03', {{0.3, {}}}, ...
            'dz02', {{0.2, {}}}, ...
            'dz08', {{0.8, {}}}, ...
            'objectiveScan', {{0.3, {'objectiveScan', true}}}, ...
            'resample222', {{0.3, {'resampleFactor', [2, 2, 2]}}}, ...
            'resample1115', {{0.3, {'resampleFactor', [1, 1, 1.5]}}}, ...
            'cubic', {{0.3, {'interpMethod', 'cubic'}}})

        % peak memory (GB, maximum resident set) of production-like tasks on
        % the abc cluster (2026-10-22): MATLAB -nojvm tasks, the largest over
        % 1 to 48 threads (compiled tasks need about 0.6 GB less). Input size
        % (one batch for large files), dz, options, measured peak. The demo
        % tile is 1800 x 512 x 500; others are cropped or repeated along z.
        measured = struct( ...
            'tiny', {{[64, 96, 80], 0.3, {}, 1.27}}, ...
            'default100', {{[1800, 512, 100], 0.3, {}, 1.89}}, ...
            'default500', {{[1800, 512, 500], 0.3, {}, 4.76}}, ...
            'default2000', {{[1800, 512, 2000], 0.3, {}, 13.57}}, ...
            'default5000', {{[1800, 512, 5000], 0.3, {}, 31.03}}, ...
            'wide250', {{[1800, 1024, 250], 0.3, {}, 4.26}}, ...
            'dz02', {{[1800, 512, 2000], 0.2, {}, 8.06}}, ...
            'dz05', {{[1800, 512, 500], 0.5, {}, 5.09}}, ...
            'dz08', {{[1800, 512, 500], 0.8, {}, 8.07}}, ...
            'single', {{[1800, 512, 2000], 0.3, {'save16bit', false}, 21.63}}, ...
            'zarrOutput100', {{[1800, 512, 100], 0.3, {'saveZarr', true}, 3.52}}, ...
            'zarrOutput2000', {{[1800, 512, 2000], 0.3, {'saveZarr', true}, 13.18}}, ...
            'zarrOutputSingle', {{[1800, 512, 500], 0.3, {'saveZarr', true, 'save16bit', false}, 8.81}}, ...
            'partialFiles500', {{[1800, 512, 500], 0.3, {'partialFile', true}, 4.66}}, ...
            'partialFiles2000', {{[1800, 512, 2000], 0.3, {'partialFile', true}, 13.18}}, ...
            'inputBbox', {{[900, 512, 500], 0.3, {'readSize', [1800, 512, 500]}, 2.88}}, ...
            'flatField', {{[1800, 512, 500], 0.3, {'preprocess', true}, 4.28}}, ...
            'outputAxisOrder', {{[1800, 512, 500], 0.3, {'outputAxisOrder', 'xyz'}, 4.49}}, ...
            'resample222', {{[1800, 512, 2000], 0.3, {'resampleFactor', [2, 2, 2]}, 8.15}}, ...
            'resample1115', {{[1800, 512, 500], 0.3, {'resampleFactor', [1, 1, 1.5]}, 4.96}}, ...
            'objectiveScan', {{[1800, 512, 500], 0.3, {'objectiveScan', true}, 7.86}}, ...
            'cubic100', {{[1800, 512, 100], 0.3, {'interpMethod', 'cubic'}, 26.81}}, ...
            'cubic500', {{[1800, 512, 500], 0.3, {'interpMethod', 'cubic'}, 72.59}}, ...
            'separate500', {{[1800, 512, 500], 0.3, {'DSRCombined', false}, 10.33}}, ...
            'separate1000', {{[1800, 512, 1000], 0.3, {'DSRCombined', false}, 30.04}}, ...
            'separate2000', {{[1800, 512, 2000], 0.3, {'DSRCombined', false}, 106.00}}, ...
            'deskewOnly500', {{[1800, 512, 500], 0.3, {'rotate', false, 'DSRCombined', false}, 9.79}}, ...
            'deskewOnly2000', {{[1800, 512, 2000], 0.3, {'rotate', false, 'DSRCombined', false}, 40.03}}, ...
            'largeFileBatch2000', {{[1025, 512, 2000], 0.3, {'outSize', [1024, 5984, 270], 'largeFile', true}, 8.47}}, ...
            'largeFileBatch5000', {{[1025, 512, 5000], 0.3, {'outSize', [1024, 13906, 270], 'largeFile', true}, 18.31}})
    end

    methods (Test, TestTags = {'Unit'})
        function outputSizeMatchesDeskewRotateFrame3D(testCase, geometry)
            dz = geometry{1};
            opts = geometry{2};
            vol = zeros(testCase.RawSize, 'uint16');
            dsr = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'reverse', true, ...
                'crop', true, 'save16bit', true, opts{:});
            [~, ~, outSize] = XR_estimateDeskewRotateMemory(testCase.RawSize, testCase.XYPixelSize, dz, ...
                'skewAngle', testCase.SkewAngle, opts{:});
            testCase.verifyEqual(outSize, size(dsr, 1 : 3));
        end

        function zStageScanUsesTheComplementaryAngle(testCase)
            % XR_deskewRotateFrame passes 90 - skewAngle for z stage scan
            vol = zeros(testCase.RawSize, 'uint16');
            dsr = deskewRotateFrame3D(vol, 90 - testCase.SkewAngle, 0.3, testCase.XYPixelSize, 'reverse', true, ...
                'crop', true, 'save16bit', true);
            [~, ~, outSize] = XR_estimateDeskewRotateMemory(testCase.RawSize, testCase.XYPixelSize, 0.3, ...
                'skewAngle', testCase.SkewAngle, 'zStageScan', true);
            testCase.verifyEqual(outSize, size(dsr, 1 : 3));
        end

        function mexPathHoldsInputInterpolatedSlicesAndOutput(testCase)
            % uint16, dz 0.3: the x step (2.34 pixels) is above 2, so one
            % slice is interpolated between neighboring slices
            sz = [1800, 512, 500];
            [~, arrayMem, outSize] = XR_estimateDeskewRotateMemory(sz, testCase.XYPixelSize, 0.3);
            testCase.verifyEqual(arrayMem, (prod(sz) + sz(1) * sz(2) * (sz(3) - 1) + prod(outSize)) * 2 / 2^30, ...
                'RelTol', 1e-12);
            % dz 0.2 needs no interpolation
            [~, arrayMem, outSize] = XR_estimateDeskewRotateMemory(sz, testCase.XYPixelSize, 0.2);
            testCase.verifyEqual(arrayMem, (prod(sz) + prod(outSize)) * 2 / 2^30, 'RelTol', 1e-12);
            % dz 0.8: x step 6.25 pixels, three interpolated slices per step
            [~, arrayMem, outSize] = XR_estimateDeskewRotateMemory(sz, testCase.XYPixelSize, 0.8);
            testCase.verifyEqual(arrayMem, (prod(sz) + 3 * sz(1) * sz(2) * (sz(3) - 1) + prod(outSize)) * 2 / 2^30, ...
                'RelTol', 1e-12);
        end

        function estimateCoversMeasuredPeak(testCase, measured)
            [sz, dz, opts, peak] = measured{:};
            memAllocate = XR_estimateDeskewRotateMemory(sz, testCase.XYPixelSize, dz, ...
                'skewAngle', testCase.SkewAngle, opts{:});
            testCase.verifyGreaterThanOrEqual(memAllocate, peak);
            % and not much more than needed
            testCase.verifyLessThanOrEqual(memAllocate, peak * 1.3 + 2.5);
        end
    end

    methods (Test, TestTags = {'Integration'})
        function preResamplingOutputSizeMatchesXRDeskewRotateFrame(testCase)
            % over 1000 slices, XR_deskewRotateFrame bins and resizes the
            % input before a combined deskew/rotation with resampling
            dataPath = testCase.createTempFolder();
            sz = [16, 24, 1001];
            fn = [dataPath, testCase.Fsn, '.tif'];
            writetiff(zeros(sz, 'uint16') + 100, fn);
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, 0.3, 'skewAngle', testCase.SkewAngle, ...
                'reverse', true, 'rotate', true, 'DSRCombined', true, 'resampleFactor', [2, 2, 2], ...
                'saveMIP', false);
            [~, ~, outSize] = XR_estimateDeskewRotateMemory(sz, testCase.XYPixelSize, 0.3, ...
                'skewAngle', testCase.SkewAngle, 'resampleFactor', [2, 2, 2]);
            testCase.verifyEqual(outSize, getImageSize([dataPath, 'DSR/', testCase.Fsn, '.tif']));
        end

        function wrapperRequestsTheEstimate(testCase)
            dataPath = testCase.createTempFolder();
            writetiff(zeros(testCase.RawSize, 'uint16'), [dataPath, testCase.Fsn, '.tif']);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            testCase.runWrapper(dataPath);
            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), testCase.estimate(testCase.RawSize), 'RelTol', 1e-12);

            % options that change the code path change the estimate
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            testCase.runWrapper(dataPath, 'save16bit', false, 'saveZarr', true);
            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), testCase.estimate(testCase.RawSize, ...
                'save16bit', false, 'saveZarr', true), 'RelTol', 1e-12);
        end

        function wrapperCountsAllPartialFiles(testCase)
            % <name>.tif + <name>_part0001.tif are read and concatenated in
            % one task; only the first part used to be counted
            dataPath = testCase.createTempFolder();
            frame = zeros(testCase.RawSize, 'uint16');
            writetiff(frame(:, :, 1 : 40), [dataPath, testCase.Fsn, '.tif']);
            writetiff(frame(:, :, 41 : end), [dataPath, testCase.Fsn, '_part0001.tif']);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            testCase.runWrapper(dataPath);
            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), testCase.estimate(testCase.RawSize, ...
                'partialFile', true), 'RelTol', 1e-12);
        end

        function wrapperUsesTheLargestFile(testCase)
            % one memAllocate for all tasks: it must fit the largest input
            dataPath = testCase.createTempFolder();
            writezarr(zeros(testCase.RawSize, 'uint16'), [dataPath, testCase.Fsn, '.zarr']);
            largeSize = testCase.RawSize .* [1, 1, 2];
            writezarr(zeros(largeSize, 'uint16'), [dataPath, testCase.FsnB, '.zarr']);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            testCase.runWrapper(dataPath, 'zarrFile', true);
            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), testCase.estimate(largeSize), 'RelTol', 1e-12);
        end

        function wrapperCountsTheCroppedInput(testCase)
            dataPath = testCase.createTempFolder();
            writetiff(zeros(testCase.RawSize, 'uint16'), [dataPath, testCase.Fsn, '.tif']);
            bbox = [1, 1, 1, 32, 96, 80];
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            testCase.runWrapper(dataPath, 'inputBbox', bbox);
            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), testCase.estimate([32, 96, 80], ...
                'readSize', testCase.RawSize), 'RelTol', 1e-12);
        end

        function largeFileBlocksRequestOneBatch(testCase)
            % XR_deskewRotateBlock tasks hold one batch (input rows with a
            % one-row border, interpolated slices, output region)
            dataPath = testCase.createTempFolder();
            zarrFn = [dataPath, testCase.Fsn, '.zarr'];
            writezarr(zeros(testCase.RawSize, 'uint16'), zarrFn, 'blockSize', [16, 64, 64]);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            XR_deskewRotateZarr(zarrFn, testCase.XYPixelSize, 0.3, 'skewAngle', testCase.SkewAngle, ...
                'reverse', true, 'batchSize', [16, 1024, 1024], 'blockSize', [16, 64, 64], 'saveMIP', false, ...
                'parseCluster', false, 'masterCompute', true);
            % 4 batches of 16 output rows; the middle ones read 18 input rows.
            % The output chunks are [16, 64, 64] cut to the 47 output slices.
            [~, ~, outSize] = XR_estimateDeskewRotateMemory(testCase.RawSize, testCase.XYPixelSize, 0.3, ...
                'skewAngle', testCase.SkewAngle);
            expected = XR_estimateDeskewRotateMemory([18, testCase.RawSize(2 : 3)], testCase.XYPixelSize, 0.3, ...
                'skewAngle', testCase.SkewAngle, 'outSize', [16, outSize(2 : 3)], 'largeFile', true, ...
                'blockSize', [16, 64, outSize(3)]);
            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), expected, 'RelTol', 1e-12);
        end
    end

    methods (Access = private)
        function [] = runWrapper(testCase, dataPath, varargin)
            XR_deskew_rotate_data_wrapper(dataPath, 'deskew', true, 'rotate', true, 'DSRCombined', true, ...
                'xyPixelSize', testCase.XYPixelSize, 'dz', 0.3, 'skewAngle', testCase.SkewAngle, ...
                'reverse', true, 'channelPatterns', {'CamA', 'CamB'}, 'parseCluster', false, ...
                'masterCompute', true, varargin{:});
        end

        function [memAllocate] = estimate(testCase, imSize, varargin)
            memAllocate = XR_estimateDeskewRotateMemory(imSize, testCase.XYPixelSize, 0.3, ...
                'skewAngle', testCase.SkewAngle, varargin{:});
        end
    end
end
