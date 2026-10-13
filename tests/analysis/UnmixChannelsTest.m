classdef UnmixChannelsTest < pktest.TestCase
    % channel unmixing (demo_useful_tools.m, "Channel Unmixing" sections):
    % out = max(0, sum_c factor(c) * (ch_c - constBackground(c))) + constOffset

    properties (Constant)
        Fsn = 'Scan_Iter_0000_0000_%s_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
    end

    methods (Test, TestTags = {'Integration'})
        function linearUnmixingAsInDemo(testCase)
            [dataPath, chA, chB] = testCase.twoChannelData(false);

            XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [1, -0.1], 'mode', 'linear', ...
                'channelPatterns', {'CamA', 'CamB'}, 'channelInd', 1, 'constBackground', [100, 100], ...
                'constOffset', 1, 'save16bit', true, 'parseCluster', false, 'masterCompute', true);

            out = readtiff([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.tif']);
            expected = uint16(max(0, (single(chA) - 100) - 0.1 * (single(chB) - 100)) + 1);
            testCase.verifyClass(out, 'uint16');
            testCase.verifyEqual(out, expected);
            % only the unmixed channel is written
            testCase.verifyNumElements(dir([dataPath, 'Unmixed/*.tif']), 1);
        end

        function unmixSecondChannelAsSingle(testCase)
            [dataPath, chA, chB] = testCase.twoChannelData(false);

            XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [-0.25, 1], 'channelInd', 2, ...
                'channelPatterns', {'CamA', 'CamB'}, 'save16bit', false, 'parseCluster', false, 'masterCompute', true);

            out = readtiff([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamB'), '.tif']);
            testCase.verifyClass(out, 'single');
            testCase.verifyImagesClose(out, max(0, single(chB) - 0.25 * single(chA)), 'relTol', 1e-6);
        end

        function flatFieldCorrectionBeforeUnmixing(testCase)
            % 1f40f4c added flat-field correction to unmixing; 824f476 fixed
            % flat-field/background images that are not single (as read
            % from uint16 tiffs), which made the mex reject them and fall
            % back to the slow MATLAB path
            [dataPath, chA, chB] = testCase.twoChannelData(false);
            ffPath = [dataPath, 'FF/'];
            mkdir(ffPath);
            [~, x] = ndgrid(1 : 48, 1 : 40);
            bk = uint16(100 * ones(48, 40));
            ffs = {uint16(100 + 200 + 800 * x / 40), uint16(100 + 1000 - 500 * x / 40)};
            ffFns = {[ffPath, 'ff_CamA.tif'], [ffPath, 'ff_CamB.tif']};
            bkFns = {[ffPath, 'bk_CamA.tif'], [ffPath, 'bk_CamB.tif']};
            for c = 1 : 2
                writetiff(ffs{c}, ffFns{c});
                writetiff(bk, bkFns{c});
            end

            consoleOut = evalc(['XR_unmix_channels_data_wrapper(dataPath, ''unmixFactors'', [1, -0.1], ', ...
                '''channelPatterns'', {''CamA'', ''CamB''}, ''FFCorrection'', true, ''lowerLimit'', 0.4, ', ...
                '''FFImagePaths'', ffFns, ''backgroundPaths'', bkFns, ''parseCluster'', false, ''masterCompute'', true)']);

            out = readtiff([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.tif']);
            corrected = cell(1, 2);
            raws = {chA, chB};
            for c = 1 : 2
                ls = single(ffs{c}) - single(bk);
                ls = max(ls / max(ls(:)), 0.4);
                corrected{c} = max((single(raws{c}) - single(bk)) ./ ls, 0);
            end
            expected = max(0, corrected{1} - 0.1 * corrected{2});
            % +/-1 from uint16 rounding
            testCase.verifyImagesClose(out, expected, 'absTol', 1, 'relTol', 1e-3);
            testCase.verifyFalse(contains(consoleOut, 'must be single'), ...
                'Flat-field correction fell back to the MATLAB implementation.');
        end

        function memAllocateIsInGigabytes(testCase)
            % 979fa55 regression: memAllocate (GB) was computed in bytes, so
            % every Slurm job requested the maximum number of CPUs
            dataPath = testCase.twoChannelData(false);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);

            XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [1, -0.1], 'channelPatterns', {'CamA', 'CamB'}, ...
                'parseCluster', false, 'masterCompute', true);

            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), prod([48, 40, 24]) * 4 * 2.5 / 1024 ^ 3, 'RelTol', 1e-12);
        end

        function zarrInputAndOutput(testCase)
            [dataPath, chA, chB] = testCase.twoChannelData(true);

            XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [1, -0.1], 'channelPatterns', {'CamA', 'CamB'}, ...
                'constBackground', 100, 'zarrFile', true, 'saveZarr', true, 'parseCluster', false, 'masterCompute', true);

            out = readzarr([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            testCase.verifyEqual(out, uint16(max(0, (single(chA) - 100) - 0.1 * (single(chB) - 100))));
        end

        function zarrFunctionMatchesFrameFunction(testCase)
            % XR_unmix_channels_zarr ("Channel Unmixing for Zarr Files")
            % processes in batches; small batches must give the same result
            [dataPath, chA, chB] = testCase.twoChannelData(true);
            zarrFns = cellfun(@(x) [dataPath, sprintf(testCase.Fsn, x), '.zarr'], {'CamA', 'CamB'}, 'unif', 0);

            XR_unmix_channels_zarr(zarrFns, [1, -0.1], 'resultDirName', 'Unmixed', 'mode', 'linear', ...
                'channelInd', 1, 'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16], ...
                'parseCluster', false, 'masterCompute', true);

            out = readzarr([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            expected = max(0, single(chA) - 0.1 * single(chB));
            testCase.verifyImagesClose(single(out), expected, 'absTol', 1, 'relTol', 1e-3);
        end

        function gaussianUnmixingBlursTheOtherChannel(testCase)
            % 'gaussian' mode filters each channel with its sigma (one row
            % per channel, [sy, sx, sz]) before the weighted sum
            [dataPath, chA, chB] = testCase.twoChannelData(true);
            zarrFns = cellfun(@(x) [dataPath, sprintf(testCase.Fsn, x), '.zarr'], {'CamA', 'CamB'}, 'unif', 0);
            sigmas = [0, 0, 0; 1, 1, 1];

            XR_unmix_channels_zarr(zarrFns, [1, -0.1], 'mode', 'gaussian', 'unmixSigmas', sigmas, ...
                'channelInd', 1, 'parseCluster', false, 'masterCompute', true);

            % the output keeps the input data type (uint16)
            out = readzarr([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            expected = max(0, single(chA) - 0.1 * imgaussfilt3(single(chB), [1, 1, 1]));
            testCase.verifyImagesClose(single(out), expected, 'absTol', 1, 'relTol', 1e-3);
        end
    end

    methods (Test, TestTags = {'Integration', 'KnownIssue'})
        function gaussianUnmixingThroughDataWrapper(testCase)
            % KNOWN ISSUE: XR_unmix_channels_data_wrapper.m:254 formats
            % unmixSigmas with mat2str_comma(unmixSigmas, 10), which only
            % works for row vectors, but the frame function indexes it per
            % channel as unmixSigmas(c, :) (XR_unmix_channels_frame.m), so a
            % per-channel sigma matrix can't be passed through the wrapper.
            % Fix: format it with strrep(mat2str(unmixSigmas), ' ', ',') as
            % XR_unmix_channels_zarr does.
            [dataPath, chA, chB] = testCase.twoChannelData(false);
            sigmas = [0, 0, 0; 1, 1, 1];

            XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [1, -0.1], 'mode', 'gaussian', ...
                'unmixSigmas', sigmas, 'channelPatterns', {'CamA', 'CamB'}, 'save16bit', false, ...
                'parseCluster', false, 'masterCompute', true);

            out = readtiff([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.tif']);
            expected = max(0, single(chA) - 0.1 * imgaussfilt3(single(chB), [1, 1, 1]));
            testCase.verifyImagesClose(out, expected, 'relTol', 1e-5);
        end
    end

    methods (Access = private)
        function [dataPath, chA, chB] = twoChannelData(testCase, zarrFile)
            % CamB bleeds into CamA
            dataPath = testCase.createTempFolder();
            chB = uint16(pktest.texture([48, 40, 24], 2, [100, 3000]));
            chA = uint16(pktest.texture([48, 40, 24], 2, [100, 1500]) + 0.1 * (single(chB) - 100));
            ims = {chA, chB};
            cams = {'CamA', 'CamB'};
            for c = 1 : 2
                fn = [dataPath, sprintf(testCase.Fsn, cams{c})];
                if zarrFile
                    writezarr(ims{c}, [fn, '.zarr'], blockSize=[16, 16, 16]);
                else
                    writetiff(ims{c}, [fn, '.tif']);
                end
            end
        end
    end
end
