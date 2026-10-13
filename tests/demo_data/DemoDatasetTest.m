classdef DemoDatasetTest < pktest.TestCase
    % checks on the published demo dataset (light sheet cell images, see
    % demos/demo_light_sheet_cell_data_downloader.m). They are skipped
    % unless the environment variable PETAKIT5D_DEMO_DATA points to the
    % PetaKit5D_demo_cell_image_dataset folder (or its parent), e.g. in
    % Jenkins. The dataset is only read; outputs go to temporary folders.
    % Needs ~16 GB RAM for 100-slice tiles (more for the 500-slice ones).

    properties (Constant)
        Fsn = 'Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t'
        SkewAngle = 32.45
        XYPixelSize = 0.108
        Dz = 0.3
    end

    properties (SetAccess = private)
        DataPath
    end

    methods (TestClassSetup)
        function findDataset(testCase)
            dataPath = getenv('PETAKIT5D_DEMO_DATA');
            testCase.assumeNotEmpty(dataPath, 'Set PETAKIT5D_DEMO_DATA to the demo dataset folder to run these tests.');
            dataPath = [strrep(dataPath, '\', '/'), '/'];
            if isfolder([dataPath, 'PetaKit5D_demo_cell_image_dataset'])
                dataPath = [dataPath, 'PetaKit5D_demo_cell_image_dataset/'];
            end
            testCase.assumeTrue(isfile([dataPath, testCase.Fsn, '.tif']), ...
                sprintf('%s does not contain the demo dataset.', dataPath));
            testCase.DataPath = dataPath;
        end
    end

    methods (Test, TestTags = {'DemoData', 'Slow'})
        function cppTiffReaderMatchesLibtiff(testCase)
            fn = [testCase.DataPath, testCase.Fsn, '.tif'];
            im = readtiff(fn);
            sz = getImageSize(fn);
            testCase.verifyEqual(size(im), sz);
            testCase.verifyEqual(sz(1 : 2), [1800, 512]);
            testCase.verifyEqual(im, readtiff_matlab(fn));
            zr = round(sz(3) / 2) + [0, 9];
            testCase.verifyEqual(readtiff(fn, 'range', zr), im(:, :, zr(1) : zr(2)));
        end

        function writersRoundTripRealData(testCase)
            im = readtiff([testCase.DataPath, testCase.Fsn, '.tif']);
            outPath = testCase.createTempFolder();
            writetiff(im, [outPath, 'a.tif']);
            testCase.verifyEqual(readtiff([outPath, 'a.tif']), im);
            writezarr(im, [outPath, 'a.zarr'], 'blockSize', [256, 256, 256]);
            testCase.verifyEqual(readzarr([outPath, 'a.zarr']), im);
        end

        function combinedDeskewRotateMatchesSeparate(testCase)
            % the comparison made in demos/demo_geometric_transformation.m
            outPath = testCase.copyTile();
            fn = [outPath, testCase.Fsn, '.tif'];
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'skewAngle', testCase.SkewAngle, ...
                'reverse', true, 'DSRCombined', false, 'rotate', true, 'DSDirName', 'DS_separated/', ...
                'DSRDirName', 'DSR_separated/');
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'skewAngle', testCase.SkewAngle, ...
                'reverse', true);
            separate = readtiff([outPath, 'DSR_separated/', testCase.Fsn, '.tif']);
            combined = readtiff([outPath, 'DSR/', testCase.Fsn, '.tif']);
            testCase.verifyEqual(size(combined), size(separate));
            % the separate path interpolates twice; on this data the
            % difference is about 2%
            testCase.verifyImagesClose(combined, separate, 'relTol', 0.04);
        end

        function omwDeconvolutionOnCPU(testCase)
            % OMW decon as in demos/demo_RL_deconvolution.m, on the CPU
            outPath = testCase.copyTile();
            XR_decon_data_wrapper({outPath}, 'resultDirName', 'matlab_decon_omw', 'xyPixelSize', testCase.XYPixelSize, ...
                'dz', testCase.Dz, 'reverse', true, 'channelPatterns', {'CamB_ch0'}, ...
                'psfFullpaths', {[testCase.DataPath, 'PSF/488_2_c.tif']}, 'dzPSF', 0.3, 'RLMethod', 'omw', ...
                'wienerAlpha', 0.005, 'OTFCumThresh', 0.87, 'skewed', true, 'background', 100, 'deconIter', 2, ...
                'save16bit', true, 'parseCluster', false, 'GPUJob', false, 'cpusPerTask', feature('numcores'));

            raw = readtiff([outPath, testCase.Fsn, '.tif']);
            decon = readtiff([outPath, 'matlab_decon_omw/', testCase.Fsn, '.tif']);
            testCase.verifyClass(decon, 'uint16');
            testCase.verifyEqual(size(decon), size(raw));
            % deconvolution sharpens: more energy in high spatial frequencies
            testCase.verifyGreaterThan(testCase.gradientEnergy(decon), 2 * testCase.gradientEnergy(max(single(raw) - 100, 0)));
        end
    end

    methods (Access = private)
        function [outPath] = copyTile(testCase)
            outPath = testCase.createTempFolder();
            copyfile([testCase.DataPath, testCase.Fsn, '.tif'], outPath);
        end

        function [e] = gradientEnergy(~, im)
            nz = size(im, 3);
            im = single(im(:, :, round(nz / 3) : round(2 * nz / 3)));
            e = mean(diff(im, 1, 1) .^ 2, 'all') + mean(diff(im, 1, 2) .^ 2, 'all');
        end
    end
end
