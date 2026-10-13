classdef DeconFrameTest < pktest.TestCase
    % XR_RLdeconFrame3D (single file decon) and
    % XR_visualize_OTF_mask_segmentation on small synthetic data

    properties (Constant)
        ImSize = [48, 48, 32]
        XYPixelSize = 0.108
        % dz * sind(32.45) < dzPSF, so RLdecon doesn't resample the PSF
        Dz = 0.2
        Background = 100
        BlurSigma = sqrt(0.8 ^ 2 + [1.5, 2] .^ 2)
    end

    methods (Test, TestTags = {'Integration'})
        function writesResultsToMatlabDeconByDefault(testCase)
            [fn, psfFn] = testCase.writeFrame();
            testCase.runFrame(fn, '', psfFn, 'DeconIter', 3);

            [dataPath, fsn] = fileparts(fn);
            deconPath = [dataPath, '/matlab_decon/'];
            out = readtiff([deconPath, fsn, '.tif']);
            testCase.verifyClass(out, 'uint16');
            testCase.verifyEqual(readtiff([deconPath, 'MIPs/', fsn, '_MIP_z.tif']), max(out, [], 3));
            [~, psfFsn] = fileparts(psfFn);
            psf = single(readtiff([deconPath, 'psfgen/', psfFsn, '_simplified.tif']));
            expected = decon_lucy_function(readtiff(fn), psf ./ sum(psf(:)), 3, 'Background', testCase.Background, ...
                'useGPU', false, 'save16bit', true);
            testCase.verifyEqual(out, expected);
        end

        function existingResultIsSkipped(testCase)
            [fn, psfFn] = testCase.writeFrame();
            deconPath = [fileparts(fn), '/decon/'];
            testCase.runFrame(fn, deconPath, psfFn, 'DeconIter', 3);
            [~, fsn] = fileparts(fn);
            before = readtiff([deconPath, fsn, '.tif']);
            % new input content: a recomputed result would differ
            testCase.writeFrame(fn);
            testCase.runFrame(fn, deconPath, psfFn, 'DeconIter', 3);
            testCase.verifyEqual(readtiff([deconPath, fsn, '.tif']), before, 'An existing result was recomputed.');
        end

        function psfGenFalseUsesTheRawPSF(testCase)
            % with psfGen=false the raw PSF is used even if a generated PSF
            % from an earlier run is cached in psfgen/ (regression test for
            % d280863 and a3eb8bd)
            bg = 100;
            psfRaw = uint16(pktest.deconPSF([31, 31, 31], [1.5, 2]) + bg);
            [fn1, psfFn] = testCase.writeFrame([], psfRaw);
            [dataPath, fsn1] = fileparts(fn1);
            fn2 = [dataPath, '/second_', fsn1, '.tif'];
            copyfile(fn1, fn2);
            deconPath = [dataPath, '/decon/'];

            testCase.runFrame(fn1, deconPath, psfFn, 'DeconIter', 3, 'psfGen', true);
            testCase.runFrame(fn2, deconPath, psfFn, 'DeconIter', 3, 'psfGen', false);

            psf = single(psfRaw);
            expected = decon_lucy_function(readtiff(fn2), psf ./ sum(psf(:)), 3, 'Background', testCase.Background, ...
                'useGPU', false, 'save16bit', true);
            out2 = readtiff([deconPath, 'second_', fsn1, '.tif']);
            testCase.verifyEqual(out2, expected);
            testCase.verifyNotEqual(out2, readtiff([deconPath, fsn1, '.tif']), ...
                'The raw and generated PSFs should give different results here.');
        end

        function savedErosionMaskCanBeReplaced(testCase)
            % rerunning the first time point replaces its saved erosion
            % mask (regression test for d654c22, rmdirs -> rmdir)
            [fn, psfFn] = testCase.writeFrame();
            [dataPath, fsn] = fileparts(fn);
            deconPath = [dataPath, '/decon/'];
            args = {'DeconIter', 2, 'EdgeErosion', 2, 'SaveMaskfile', true};
            testCase.runFrame(fn, deconPath, psfFn, args{:});
            maskFn = [deconPath, 'Masks/', fsn, '_eroded.zarr'];
            testCase.assertTrue(exist(maskFn, 'dir') == 7);
            delete([deconPath, fsn, '.tif']);

            testCase.runFrame(fn, deconPath, psfFn, args{:});
            testCase.verifyTrue(exist([deconPath, fsn, '.tif'], 'file') == 2);
            testCase.verifyEqual(readzarr(maskFn) > 0, decon_mask_edge_erosion(readtiff(fn) > 0, 2));
        end

        function otfMaskVisualizationSavesTheFigure(testCase)
            % demo step: check the OTF mask threshold before decon
            [~, psfFn] = testCase.writeFrame();
            evalc('XR_visualize_OTF_mask_segmentation(psfFn, 0.88, false, ''visible'', false, ''saveFig'', true)');
            [psfPath, psfFsn] = fileparts(psfFn);
            testCase.verifyTrue(exist([psfPath, '/figures/', psfFsn, '.png'], 'file') == 2);
            % automatic skewed space detection and a minimum intensity threshold
            evalc('XR_visualize_OTF_mask_segmentation(psfFn, 0.65, [], ''minIntThrsh'', 1e-3, ''visible'', false)');
        end
    end

    methods (Test, TestTags = {'Integration', 'KnownIssue'})
        function overwriteRecomputesTheResult(testCase)
            % KNOWN ISSUE: XR_RLdeconFrame3D.m:213 honors Overwrite, but
            % RLdecon.m:195-198 returns early whenever the output exists, so
            % Overwrite=true leaves the old result (non-large path). Fix:
            % pass Overwrite to RLdecon or delete the old output first.
            [fn, psfFn] = testCase.writeFrame();
            deconPath = [fileparts(fn), '/decon/'];
            testCase.runFrame(fn, deconPath, psfFn, 'DeconIter', 3);
            testCase.writeFrame(fn);
            testCase.runFrame(fn, deconPath, psfFn, 'DeconIter', 3, 'Overwrite', true);

            [~, fsn] = fileparts(fn);
            [~, psfFsn] = fileparts(psfFn);
            psf = single(readtiff([deconPath, 'psfgen/', psfFsn, '_simplified.tif']));
            expected = decon_lucy_function(readtiff(fn), psf ./ sum(psf(:)), 3, 'Background', testCase.Background, ...
                'useGPU', false, 'save16bit', true);
            testCase.verifyEqual(readtiff([deconPath, fsn, '.tif']), expected, 'Overwrite=true kept the old result.');
        end

        function flipZstackFlipsTheResult(testCase)
            % KNOWN ISSUE: XR_RLdeconFrame3D.m:260-269 doesn't pass
            % flipZstack to RLdecon, so the regular path never flips the
            % stack (the saved erosion mask and the large-file path do).
            % Fix: add 'flipZstack', flipZstack to the RLdecon call.
            [fn, psfFn] = testCase.writeFrame();
            dataPath = fileparts(fn);
            testCase.runFrame(fn, [dataPath, '/noflip/'], psfFn, 'DeconIter', 3, 'flipZstack', false);
            testCase.runFrame(fn, [dataPath, '/flip/'], psfFn, 'DeconIter', 3, 'flipZstack', true);
            [~, fsn] = fileparts(fn);
            % the PSF is symmetric in z, so decon commutes with the flip
            testCase.verifyImagesClose(readtiff([dataPath, '/flip/', fsn, '.tif']), ...
                flip(readtiff([dataPath, '/noflip/', fsn, '.tif']), 3), 'relTol', 1e-3, 'name', 'flipped result');
        end

        function simplifiedInfoFilesGoToTheDeconInfoFolder(testCase)
            % KNOWN ISSUE: RLdecon.m:537 builds the folder with
            % sprintf('%s/%s/decon_info/', deconPath) (two %s, one value),
            % so the info .mat/.txt files land next to the results instead
            % of in decon_info/. Fix: sprintf('%s/decon_info/', deconPath).
            [fn, psfFn] = testCase.writeFrame();
            deconPath = [fileparts(fn), '/decon/'];
            testCase.runFrame(fn, deconPath, psfFn, 'DeconIter', 3);
            [~, fsn] = fileparts(fn);
            testCase.verifyTrue(exist([deconPath, 'decon_info/', fsn, '_info.mat'], 'file') == 2, ...
                sprintf('decon_info/%s_info.mat is missing (found in the result folder: %d).', fsn, ...
                exist([deconPath, fsn, '_info.mat'], 'file') == 2));
        end

        function allZeroFrameIsDeconvolved(testCase)
            % KNOWN ISSUE: for an all-zero image RLdecon.m:368 skips the
            % decon, but RLdecon.m:536-544 still saves err_mat and iter_run
            % for 'simplified', which are undefined, so empty tiles (e.g.,
            % outside the sample) make the whole run error. Fix: only save
            % the info if the decon ran.
            [fn, psfFn] = testCase.writeFrame([], [], zeros(testCase.ImSize, 'uint16'));
            deconPath = [fileparts(fn), '/decon/'];
            try
                testCase.runFrame(fn, deconPath, psfFn, 'DeconIter', 3);
            catch ME
                testCase.verifyFail(sprintf('Deconvolving an all-zero frame errors: %s', ME.message));
                return;
            end
            [~, fsn] = fileparts(fn);
            testCase.verifyEqual(readtiff([deconPath, fsn, '.tif']), zeros(testCase.ImSize, 'uint16'));
        end
    end

    methods (Access = private)
        function [fn, psfFn] = writeFrame(testCase, fn, psf, im)
            % write a frame of blurred beads (and a PSF in <root>/psf/);
            % with a given fn, only overwrite that frame with new beads
            if nargin < 4 || isempty(im)
                pos = [8 + 32 * rand(6, 2), 6 + 20 * rand(6, 1)];
                im = uint16(pktest.beads(testCase.ImSize, pos, testCase.BlurSigma, 1000, testCase.Background));
            end
            psfFn = '';
            if nargin < 2 || isempty(fn)
                root = testCase.createTempFolder();
                mkdir([root, 'data']);
                mkdir([root, 'psf']);
                fn = [root, 'data/Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_000x_000y_000z_0000t.tif'];
                psfFn = [root, 'psf/PSF_CamA.tif'];
                if nargin < 3 || isempty(psf)
                    psf = uint16(pktest.deconPSF([31, 31, 31], [1.5, 2]));
                end
                writetiff(psf, psfFn);
            end
            writetiff(im, fn);
        end

        function [] = runFrame(testCase, fn, deconPath, psfFn, varargin)
            args = [{'PSFfile', psfFn, 'dzPSF', testCase.Dz, 'Background', testCase.Background, ...
                'RLMethod', 'simplified', 'EdgeErosion', 0, 'parseCluster', false, 'GPUJob', false}, varargin];
            evalc('XR_RLdeconFrame3D(fn, testCase.XYPixelSize, testCase.Dz, deconPath, args{:})');
        end
    end
end
