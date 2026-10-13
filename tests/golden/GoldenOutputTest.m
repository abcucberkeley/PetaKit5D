classdef GoldenOutputTest < pktest.TestCase
    % compare outputs of the main algorithms on fixed synthetic inputs with
    % reference outputs saved in tests/golden/data. This catches unintended
    % numerical changes that the property tests may not (e.g., different
    % interpolation weights, iteration updates or blending weights).
    %
    % The references are made on Linux with update_golden_data.m. Other
    % platforms must match within a relative L2 error of 1e-4, which allows
    % floating point differences but not algorithm changes (>= 1e-3).
    %
    % If an output change is intended, regenerate the references with
    %   update_golden_data('<scenario>')
    % and commit them together with the change.

    properties (TestParameter)
        scenario = {'deskew_rotate_sample_scan', 'deskew_rotate_skewed_interp', 'deskew_rotate_resampled', ...
            'decon_rl_simplified', 'decon_omw', 'psf_gen_masked', 'stitch_feather_xcorr', 'point_detection'}
    end

    methods (Test, TestTags = {'Golden'})
        function outputMatchesReference(testCase, scenario)
            refFn = fullfile(fileparts(mfilename('fullpath')), 'data', [scenario, '.mat']);
            testCase.assertTrue(isfile(refFn), sprintf(['Reference output %s is missing; generate it with ', ...
                'update_golden_data(''%s'').'], refFn, scenario));
            if strcmp(scenario, 'point_detection')
                testCase.assumeMexAvailable('compute_gradient_hessian_mex', 'adtest_mex', 'KDTreeBallQuery');
            end

            ref = load(refFn);
            actual = GoldenOutputTest.computeScenario(scenario, testCase.createTempFolder());

            fns = fieldnames(ref.expected);
            testCase.verifyEqual(sort(fieldnames(actual)), sort(fns));
            for i = 1 : numel(fns)
                if ~isfield(actual, fns{i})
                    continue;
                end
                testCase.verifyClass(actual.(fns{i}), class(ref.expected.(fns{i})), fns{i});
                testCase.verifyImagesClose(actual.(fns{i}), ref.expected.(fns{i}), 'relTol', 1e-4, ...
                    'name', sprintf('%s.%s (reference from %s, %s)', scenario, fns{i}, ref.info.platform, ref.info.version));
            end
        end
    end

    methods (Static)
        function [out] = computeScenario(name, dataPath)
            % outputs of one scenario as a struct of numeric arrays. Inputs
            % are generated here with a fixed seed.
            rng(1, 'twister');
            px = 0.108;
            angle = 32.45;
            switch name
                % the float path, rounded to keep the references small; the
                % uint16 path rounds differently on macOS (see
                % DeskewRotateFrame3DTest/save16bitRoundsToNearest)
                case 'deskew_rotate_sample_scan'
                    vol = single(GoldenOutputTest.skewedVolume());
                    out.dsr = uint16(round(deskewRotateFrame3D(vol, angle, 0.2, px, 'reverse', true)));

                case 'deskew_rotate_skewed_interp'
                    % dz = 0.3 interpolates in skewed space before deskew/rotation
                    vol = single(GoldenOutputTest.skewedVolume());
                    out.dsr = uint16(round(deskewRotateFrame3D(vol, angle, 0.3, px, 'reverse', true)));

                case 'deskew_rotate_resampled'
                    vol = single(GoldenOutputTest.skewedVolume());
                    out.dsr = deskewRotateFrame3D(vol, angle, 0.2, px, 'reverse', true, 'resampleFactor', [2, 2, 2]);

                case 'decon_rl_simplified'
                    [im, psf] = GoldenOutputTest.blurredBeads();
                    out.decon = decon_lucy_function(im, psf, 10, 'Background', 100, 'useGPU', false, 'save16bit', false);

                case 'decon_omw'
                    [im, psf] = GoldenOutputTest.blurredBeads();
                    b = omw_backprojector_generation(psf, 0.005, false, 'OTFCumThresh', 0.9);
                    out.backProjector = single(b);
                    out.decon = decon_lucy_omw_function(im, psf, b, 2, 'Background', 100, 'useGPU', false, 'save16bit', false);

                case 'psf_gen_masked'
                    % measured-like PSF: camera offset, noise and a tilted sheet
                    psf = pktest.beads([63, 63, 41], [32, 32, 21], [1.6, 2.2], 2000, 100);
                    psf = psf + pktest.beads([63, 63, 41], [32, 32, 21], [12, 1.2], 50, 0);
                    psf = single(round(psf + 3 * randn(size(psf))));
                    out.psf = single(psf_gen_new(psf, 0.2, 0.2, 1.5, 'masked'));

                case 'stitch_feather_xcorr'
                    out.stitched = GoldenOutputTest.stitchTwoTiles(dataPath);

                case 'point_detection'
                    pos = [16.3, 20.7, 18.2; 16, 48, 32; 48, 24, 40; 40, 40, 20; 32, 32, 48];
                    vol = double(pktest.beads([64, 64, 64], pos, 1.5, 200, 100)) + 5 * randn(64, 64, 64);
                    clear XR_pointSourceDetection3D
                    ps = XR_pointSourceDetection3D(vol, [1.5, 1.5], 'BackgroundCorrection', false);
                    [~, order] = sortrows([ps.z(:), ps.y(:), ps.x(:)]);
                    out.position = [ps.y(order(:)'); ps.x(order(:)'); ps.z(order(:)')]';
                    out.amplitude = reshape(ps.A(order), [], 1);
                    out.background = reshape(ps.c(order), [], 1);

                otherwise
                    error('Unknown golden output scenario %s', name);
            end
        end
    end

    methods (Static, Access = private)
        function [vol] = skewedVolume()
            % beads plus smooth texture in skewed space, first dimension a
            % multiple of 8 (see DeskewMexSourceTest)
            [xr, zr] = ndgrid([12, 24, 36], [8, 20, 32]);
            raw = [4 + 2.6 * mod((0 : 8)' * 4, 9), xr(:), zr(:)];
            vol = pktest.beads([32, 48, 40], raw, [1.5, 2], 1000, 100) + pktest.texture([32, 48, 40], 3, [0, 200]);
            vol = uint16(vol);
        end

        function [im, psf] = blurredBeads()
            psf = pktest.beads([21, 21, 21], [11, 11, 11], [1.5, 2.5], 1, 0);
            psf = psf / sum(psf(:));
            obj = pktest.beads([32, 32, 32], [10, 12, 9; 22, 20, 24; 16, 8, 20], 0.8, 5000, 0) + ...
                pktest.texture([32, 32, 32], 2, [0, 300]);
            im = imfilter(obj, double(psf), 'circular', 'conv') + 100;
            im = single(round(im + 2 * randn(size(im))));
        end

        function [stitched] = stitchTwoTiles(dataPath)
            % two overlapping tiles where tile 2 is 400 brighter, so the
            % feather weights show in the overlap; tile 2 is listed 2
            % pixels off so xcorr has to correct it
            gt = pktest.texture([32, 104, 32], 3, [1000, 5000]);
            fsn = 'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_%03dx_000y_000z_0000t.tif';
            x0 = [0, 40];
            fns = cell(2, 1);
            for i = 1 : 2
                fns{i} = sprintf(fsn, i - 1);
                tile = gt(:, x0(i) + 1 : x0(i) + 64, :) + 400 * (i - 1);
                writetiff(uint16(tile), [dataPath, fns{i}]);
            end
            px = 0.108;
            t = table(fns, [0; 42 * px], [0; 0], [0; 0], 'VariableNames', {'Filename', 'StageX_um_', 'StageY_um_', 'StageZ_um_'});
            writetable(t, [dataPath, 'ImageList.csv']);
            XR_matlab_stitching_wrapper(dataPath, [dataPath, 'ImageList.csv'], 'DSR', true, 'resampleType', 'isotropic', ...
                'xyPixelSize', px, 'dz', px, 'axisOrder', 'x,y,z', 'xcorrShift', true, 'xcorrDownsample', [1, 1, 1], ...
                'xyMaxOffset', 10, 'zMaxOffset', 4, 'blendMethod', 'feather', 'channelPatterns', {'CamA_ch0'}, ...
                'primaryCh', 'CamA_ch0', 'parseCluster', false);
            stitched = readzarr([dataPath, 'matlab_stitch/Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs.zarr']);
        end
    end
end
