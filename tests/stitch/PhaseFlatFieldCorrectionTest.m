classdef PhaseFlatFieldCorrectionTest < pktest.TestCase
    % XR_phase_image_flat_field_correction (the preprocessing of 2D phase
    % tiles in demo_phase_and_2photon_stitching.m): subtract a Gaussian
    % blurred image (sigma 100), add a constant, estimate the flat field of
    % all tiles with BaSiC, and divide by it.

    methods (Test, TestTags = {'Integration'})
        function correctsTilesWithEstimatedFlatField(testCase)
            dataPath = testCase.createTempFolder();
            outPath = [dataPath, 'ff_corrected/'];
            mkdir(outPath);
            const = 1100;

            % vignetting shared by all tiles: 1 at the center, 0.6 at corners
            [x, y] = meshgrid(linspace(-1, 1, 128));
            shading = single(1 - 0.2 * (x .^ 2 + y .^ 2));
            nT = 8;
            fns = cell(nT, 1);
            raw = cell(nT, 1);
            for i = 1 : nT
                raw{i} = uint16((pktest.texture([128, 128], 1.5, [500, 3000])) .* shading);
                fns{i} = sprintf(['%sScan_Iter_0000_CamB_ch0_CAM1_stack0000_642nm_0000000msec_', ...
                    '0000006528774msecAbs_000x_%03dy_000z_0000t.tif'], dataPath, i - 1);
                writetiff(raw{i}, fns{i});
            end

            XR_phase_image_flat_field_correction(fns, outPath, const);

            [~, fsn1] = fileparts(fns{1});
            ffFn = sprintf('%sBaSiC_flatfield/%s_flatfield.mat', outPath, fsn1);
            testCase.assertTrue(exist(ffFn, 'file') == 2, 'The flat field file was not saved.');
            a = load(ffFn);
            testCase.verifyEqual(size(a.flatfield), [128, 128]);
            % the estimate follows the vignetting (r is ~0.9 here; a broken
            % estimate would be close to 0)
            r = corrcoef(double(a.flatfield(:)), double(shading(:)));
            testCase.verifyGreaterThan(r(1, 2), 0.8, 'The estimated flat field does not follow the shading.');

            for i = 1 : nT
                [~, fsn] = fileparts(fns{i});
                out = readtiff([outPath, fsn, '.tif']);
                imf = single(raw{i});
                expected = uint16((imf - imgaussfilt(imf, 100) + const) ./ a.flatfield);
                testCase.verifyClass(out, 'uint16');
                testCase.verifyImagesClose(out, expected, 'relTol', 1e-6, 'absTol', 1, 'name', fsn);
            end
        end
    end
end
