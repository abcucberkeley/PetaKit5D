classdef ImarisConversionTest < pktest.TestCase
    % XR_imaris_conversion_data_wrapper ("Imaris converter" in
    % demo_useful_tools.m). It shells out to the bundled
    % Parallel_Imaris_Writer executable for this platform (via system()),
    % so this also checks that the executable and its libraries run.
    %
    % The image data in the .ims file is LZ4-compressed, which MATLAB's
    % HDF5 library can't decode without the filter plugin, so the tests
    % check structure, attributes and histograms instead of voxels.

    properties (Constant)
        ImSize = [40, 32, 20]
        PixelSizes = [0.108, 0.108, 0.2]
        Fsn = 'Scan_Iter_%04d_0000_%s_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
    end

    methods (Test, TestTags = {'Integration'})
        function convertsChannelsAndTimePoints(testCase)
            % two channels (8c93f1b regression: with several channel
            % patterns no file was selected) and two time points
            [dataPath, maxVals] = testCase.writeTimeSeries();
            imsFn = testCase.convert(dataPath, 'tiff');

            for t = 0 : 1
                for c = 0 : 1
                    g = sprintf('/DataSet/ResolutionLevel 0/TimePoint %d/Channel %d', t, c);
                    info = h5info(imsFn, g);
                    testCase.verifyTrue(any(strcmp({info.Datasets.Name}, 'Data')), sprintf('%s has no Data.', g));
                    % the histogram range covers the channel's intensities,
                    % which tells the channels and time points apart
                    testCase.verifyEqual(str2double(char(h5readatt(imsFn, g, 'HistogramMax')')), ...
                        maxVals(t + 1, c + 1), 'RelTol', 1e-6);
                end
            end
            timeInfo = '/DataSetInfo/TimeInfo';
            testCase.verifyEqual(str2double(char(h5readatt(imsFn, timeInfo, 'DatasetTimePoints')')), 2);
            testCase.verifyTrue(exist([dataPath, 'imaris/parameters.json'], 'file') == 2);
        end

        function converterArgumentsHaveNoSpaces(testCase)
            % c9f2f9c regression: num2str pads numbers of different widths
            % with spaces ([8, 8, 16] became '8, 8,16,'), which split the
            % converter's -b/-v arguments. The command is captured with a
            % framework spy instead of running the converter.
            dataPath = testCase.writeTimeSeries();
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);

            XR_imaris_conversion_data_wrapper({dataPath}, 'channelPatterns', {'CamA', 'CamB'}, ...
                'pixelSizes', [0.108, 0.108, 10.5], 'blockSize', [8, 8, 128], 'parseCluster', false, ...
                'masterCompute', true, 'cpusPerTask', 2);

            calls = spy.getCalls();
            cmd = calls{1}.funcStrs{1};
            testCase.verifySubstring(cmd, ' -b 8,8,128, ');
            voxelArg = regexp(cmd, ' -v (\S+) ', 'tokens', 'once');
            testCase.assertNotEmpty(voxelArg, 'No -v argument in the converter command.');
            testCase.verifyEqual(str2double(strsplit(strip(voxelArg{1}, ','), ',')), [0.108, 0.108, 10.5], 'RelTol', 1e-12);
            testCase.verifySubstring(cmd, ' -P "CamA,CamB" ');
        end

        function voxelSizesAreWritten(testCase)
            dataPath = testCase.writeTimeSeries();
            imsFn = testCase.convert(dataPath, 'tiff');
            g = '/DataSetInfo/Image';
            dims = cellfun(@(x) str2double(char(h5readatt(imsFn, g, x)')), {'X', 'Y', 'Z'});
            extMax = cellfun(@(x) str2double(char(h5readatt(imsFn, g, x)')), {'ExtMax0', 'ExtMax1', 'ExtMax2'});
            extMin = cellfun(@(x) str2double(char(h5readatt(imsFn, g, x)')), {'ExtMin0', 'ExtMin1', 'ExtMin2'});
            testCase.verifyEqual((extMax - extMin) ./ dims, testCase.PixelSizes, 'RelTol', 1e-6);
        end

        function zarrInput(testCase)
            dataPath = testCase.createTempFolder();
            for c = {'CamA', 'CamB'}
                writezarr(uint16(100 + randi(100, testCase.ImSize)), [dataPath, sprintf(testCase.Fsn, 0, c{1}), '.zarr']);
            end
            imsFn = testCase.convert(dataPath, 'zarr');
            info = h5info(imsFn, '/DataSet/ResolutionLevel 0/TimePoint 0');
            testCase.verifyNumElements(info.Groups, 2);
        end
    end

    methods (Test, TestTags = {'Integration', 'KnownIssue'})
        function imageSizeIsNotPaddedToBlocks(testCase)
            % KNOWN ISSUE: Parallel_Imaris_Writer records the image size
            % padded to multiples of blockSize (48 x 32 x 24 for a 40 x 32 x
            % 20 image with blockSize 16), and its histograms count the
            % padding, so Imaris shows zero borders and a histogram minimum
            % of 0. The ImageSize attributes should be the true size.
            dataPath = testCase.writeTimeSeries();
            imsFn = testCase.convert(dataPath, 'tiff');
            g = '/DataSet/ResolutionLevel 0/TimePoint 0/Channel 0';
            sz = cellfun(@(x) str2double(char(h5readatt(imsFn, g, x)')), {'ImageSizeX', 'ImageSizeY', 'ImageSizeZ'});
            % don't assume which of the first two axes Imaris calls X
            testCase.verifyEqual(sort(sz(1 : 2)), sort(testCase.ImSize(1 : 2)));
            testCase.verifyEqual(sz(3), testCase.ImSize(3));
            hist = h5read(imsFn, [g, '/Histogram']);
            testCase.verifyEqual(sum(double(hist)), prod(testCase.ImSize), 'The histogram counts padded voxels.');
        end
    end

    methods (Access = private)
        function [dataPath, maxVals] = writeTimeSeries(testCase)
            % two time points x two channels, with distinct intensity ranges
            % (maxima above 255, below which the histogram range is 0-256)
            dataPath = testCase.createTempFolder();
            cams = {'CamA', 'CamB'};
            maxVals = zeros(2, 2);
            for t = 0 : 1
                for c = 1 : 2
                    offset = 1000 * t + 500 * (c - 1);
                    im = uint16(400 + offset + randi(100, testCase.ImSize));
                    im(1) = 500 + offset;
                    maxVals(t + 1, c) = 500 + offset;
                    writetiff(im, [dataPath, sprintf(testCase.Fsn, t, cams{c}), '.tif']);
                end
            end
        end

        function [imsFn] = convert(testCase, dataPath, inputType)
            XR_imaris_conversion_data_wrapper({dataPath}, 'channelPatterns', {'CamA', 'CamB'}, ...
                'pixelSizes', testCase.PixelSizes, 'zarrFile', strcmp(inputType, 'zarr'), 'blockSize', [16, 16, 16], ...
                'parseCluster', false, 'masterCompute', true, 'cpusPerTask', 2);
            imsFn = [dataPath, 'imaris/', sprintf(testCase.Fsn, 0, 'CamA'), '.ims'];
            testCase.assertTrue(exist(imsFn, 'file') == 2, 'The Imaris file was not written.');
            h5info(imsFn);
        end
    end
end
