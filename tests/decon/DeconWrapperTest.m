classdef DeconWrapperTest < pktest.TestCase
    % XR_decon_data_wrapper on small synthetic light-sheet datasets: output
    % files, equality with the decon functions, channel to PSF mapping,
    % erosion by the first time point, large-file paths, and the workflows
    % in demo_RL_deconvolution, demo_widefield_and_confocal_deconvolution
    % and demo_large_scale_processing.

    properties (Constant)
        ImSize = [48, 48, 32]
        XYPixelSize = 0.108
        % dz == dzPSF, so neither the wrapper nor RLdecon resamples the PSF
        Dz = 0.2
        Background = 100
        BlurSigma = sqrt(0.8 ^ 2 + [1.5, 2] .^ 2)
    end

    methods (Test, TestTags = {'Integration'})
        function simplifiedWritesOutputsThatMatchTheFunction(testCase)
            d = testCase.createDataset({'A'}, 1);
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 5);

            deconPath = [d.dataPath, 'matlab_decon/'];
            out = readtiff([deconPath, d.fsns{1}, '.tif']);
            testCase.verifyClass(out, 'uint16');
            testCase.verifyEqual(size(out), testCase.ImSize);
            testCase.verifyEqual(readtiff([deconPath, 'MIPs/', d.fsns{1}, '_MIP_z.tif']), max(out, [], 3));

            psfgenFn = [deconPath, 'psfgen/', d.psfFsns{1}, '_simplified.tif'];
            testCase.verifyTrue(exist(psfgenFn, 'file') == 2, 'The generated PSF was not saved.');
            testCase.verifyTrue(exist([deconPath, 'parameters.mat'], 'file') == 2);
            params = jsondecode(fileread([deconPath, 'parameters.json']));
            testCase.verifyEqual(params.RLMethod, 'simplified');

            psf = single(readtiff(psfgenFn));
            expected = decon_lucy_function(readtiff(d.rawFns{1}), psf ./ sum(psf(:)), 5, ...
                'Background', testCase.Background, 'useGPU', false, 'save16bit', true);
            testCase.verifyEqual(out, expected);
        end

        function omwWritesBackProjectorFilesAndMatchesTheFunction(testCase)
            d = testCase.createDataset({'A'}, 1);
            testCase.runWrapper(d, 'RLMethod', 'omw', 'deconIter', 2, 'wienerAlpha', 0.005, 'OTFCumThresh', 0.9, 'skewed', false);

            psfgenPath = [d.dataPath, 'matlab_decon/psfgen/', d.psfFsns{1}];
            bpFn = [psfgenPath, '_back_projector_alpha_0.005000_otf_cum_thresh_0.900000.tif'];
            files = {[psfgenPath, '_omw.tif'], bpFn, [psfgenPath, '_OTF.tif'], ...
                [psfgenPath, '_OTF_mask_otf_cum_thresh_0.900000.tif'], ...
                [psfgenPath, '_back_projector_alpha_0.005000_otf_cum_thresh_0.900000_figure.png']};
            for i = 1 : numel(files)
                testCase.verifyTrue(exist(files{i}, 'file') == 2, sprintf('Missing %s', files{i}));
            end

            psf = single(readtiff([psfgenPath, '_omw.tif']));
            expected = decon_lucy_omw_function(readtiff(d.rawFns{1}), psf ./ sum(psf(:)), readtiff(bpFn), 2, ...
                'Background', testCase.Background, 'useGPU', false, 'save16bit', true);
            testCase.verifyEqual(readtiff([d.dataPath, 'matlab_decon/', d.fsns{1}, '.tif']), expected);
        end

        function save16bitFalseWritesSingle(testCase)
            d = testCase.createDataset({'A'}, 1);
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 3, 'save16bit', false);
            out = readtiff([d.dataPath, 'matlab_decon/', d.fsns{1}, '.tif']);
            testCase.verifyClass(out, 'single');

            psf = single(readtiff([d.dataPath, 'matlab_decon/psfgen/', d.psfFsns{1}, '_simplified.tif']));
            expected = decon_lucy_function(readtiff(d.rawFns{1}), psf ./ sum(psf(:)), 3, ...
                'Background', testCase.Background, 'useGPU', false, 'save16bit', false);
            testCase.verifyEqual(out, expected);
        end

        function zarrOutputMatchesTiffOutput(testCase)
            d = testCase.createDataset({'A'}, 1);
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 3, 'resultDirName', 'decon_tif');
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 3, 'resultDirName', 'decon_zarr', 'saveZarr', true);
            testCase.verifyEqual(readzarr([d.dataPath, 'decon_zarr/', d.fsns{1}, '.zarr']), ...
                readtiff([d.dataPath, 'decon_tif/', d.fsns{1}, '.tif']));
        end

        function eachChannelUsesItsOwnPSFAndBackground(testCase)
            % PSFs and per channel parameters map to files in the order of
            % channelPatterns
            d = testCase.createDataset({'A', 'B'}, 1, {[1.5, 2], [2.5, 3]});
            bg = [100, 80];
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 3, 'background', bg);

            for c = 1 : 2
                psf = single(readtiff(sprintf('%smatlab_decon/psfgen/%s_simplified.tif', d.dataPath, d.psfFsns{c})));
                expected = decon_lucy_function(readtiff(d.rawFns{c}), psf ./ sum(psf(:)), 3, ...
                    'Background', bg(c), 'useGPU', false, 'save16bit', true);
                testCase.verifyEqual(readtiff([d.dataPath, 'matlab_decon/', d.fsns{c}, '.tif']), expected, ...
                    sprintf('Channel %d was not deconvolved with its own PSF and background.', c));
            end
        end

        function omwParametersArePerChannel(testCase)
            % demo_large_scale_processing passes wienerAlpha and
            % OTFCumThresh per channel
            d = testCase.createDataset({'A', 'B'}, 1);
            testCase.runWrapper(d, 'RLMethod', 'omw', 'deconIter', 1, 'wienerAlpha', [0.05, 0.005], ...
                'OTFCumThresh', [0.88, 0.9], 'skewed', false);
            psfgenPath = [d.dataPath, 'matlab_decon/psfgen/'];
            testCase.verifyTrue(exist(sprintf('%s%s_back_projector_alpha_0.050000_otf_cum_thresh_0.880000.tif', ...
                psfgenPath, d.psfFsns{1}), 'file') == 2, 'Channel A did not use its own OMW parameters.');
            testCase.verifyTrue(exist(sprintf('%s%s_back_projector_alpha_0.005000_otf_cum_thresh_0.900000.tif', ...
                psfgenPath, d.psfFsns{2}), 'file') == 2, 'Channel B did not use its own OMW parameters.');
        end

        function erodeByFirstTimePointAppliesItsMaskToAllFiles(testCase)
            % the eroded mask of the first time point is saved and applied
            % to every file of the dataset
            d = testCase.createDataset({'A'}, 2);
            raw0 = readtiff(d.rawFns{1});
            raw0(:, 1 : 10, :) = 0;
            writetiff(raw0, d.rawFns{1});
            e = 2;
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 3, 'edgeErosion', e, 'erodeByFTP', true);

            mask = readzarr(sprintf('%smatlab_decon/Masks/%s_eroded.zarr', d.dataPath, d.fsns{1})) > 0;
            testCase.verifyEqual(mask, decon_mask_edge_erosion(raw0 > 0, e));
            for t = 1 : 2
                out = readtiff([d.dataPath, 'matlab_decon/', d.fsns{t}, '.tif']);
                testCase.verifyEqual(nnz(out(~mask)), 0, sprintf('Time point %d is not zero outside the first time point mask.', t));
                testCase.verifyGreaterThan(nnz(out(mask)), 0);
            end
        end

        function largeFileSingleBatchMatchesRegular(testCase)
            % with one batch the in-memory large-file path runs the same
            % decon on the same data, including the edge erosion
            d = testCase.createDataset({'A'}, 1);
            args = {'RLMethod', 'simplified', 'deconIter', 5, 'edgeErosion', 3, 'erodeByFTP', false};
            testCase.runWrapper(d, args{:}, 'resultDirName', 'regular');
            testCase.runWrapper(d, args{:}, 'resultDirName', 'large', 'largeFile', true, 'batchSize', [64, 64, 64]);
            large = readtiff([d.dataPath, 'large/', d.fsns{1}, '.tif']);
            testCase.verifyEqual(large, readtiff([d.dataPath, 'regular/', d.fsns{1}, '.tif']));
            testCase.verifyEqual(readtiff([d.dataPath, 'large/MIPs/', d.fsns{1}, '_MIP_z.tif']), max(large, [], 3));
        end

        function largeFileMultipleBatchesMatchRegular(testCase)
            % batches are deconvolved with a border of half the PSF; with 2
            % iterations (no acceleration) the stitched result matches the
            % whole volume decon away from the image border, where the
            % periodic FFT boxes differ
            d = testCase.createDataset({'A'}, 1);
            args = {'RLMethod', 'simplified', 'deconIter', 2, 'edgeErosion', 0};
            testCase.runWrapper(d, args{:}, 'resultDirName', 'regular');
            testCase.runWrapper(d, args{:}, 'resultDirName', 'large', 'largeFile', true, 'batchSize', [24, 24, 16]);
            regular = single(readtiff([d.dataPath, 'regular/', d.fsns{1}, '.tif']));
            large = single(readtiff([d.dataPath, 'large/', d.fsns{1}, '.tif']));
            b = 8;
            testCase.verifyImagesClose(large(b + 1 : end - b, b + 1 : end - b, b + 1 : end - b), ...
                regular(b + 1 : end - b, b + 1 : end - b, b + 1 : end - b), 'relTol', 0.01, 'absTol', 2);
        end

        function largeFileInPlaceMatchesRegular(testCase)
            % zarr in, zarr out (demo_large_scale_processing); one batch
            d = testCase.createDataset({'A'}, 1, {}, true);
            args = {'RLMethod', 'simplified', 'deconIter', 3, 'zarrFile', true, 'saveZarr', true};
            testCase.runWrapper(d, args{:}, 'resultDirName', 'regular');
            testCase.runWrapper(d, args{:}, 'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', ...
                'batchSize', [64, 64, 64], 'blockSize', [16, 16, 16]);
            large = readzarr([d.dataPath, 'large/', d.fsns{1}, '.zarr']);
            testCase.verifyEqual(large, readzarr([d.dataPath, 'regular/', d.fsns{1}, '.zarr']));
            testCase.verifyEqual(readtiff([d.dataPath, 'large/MIPs/', d.fsns{1}, '_MIP_z.tif']), max(large, [], 3));
        end

        function largeFileInPlaceMultipleBatchesMatchRegular(testCase)
            % the demo_large_scale_processing setup: zarr chunks deconvolved
            % batch by batch into a zarr, with edge erosion
            d = testCase.createDataset({'A'}, 1, {}, true);
            e = 3;
            args = {'RLMethod', 'simplified', 'deconIter', 2, 'edgeErosion', e, 'erodeByFTP', false, ...
                'zarrFile', true, 'saveZarr', true};
            testCase.runWrapper(d, args{:}, 'resultDirName', 'regular');
            testCase.runWrapper(d, args{:}, 'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', ...
                'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16]);
            regular = single(readzarr([d.dataPath, 'regular/', d.fsns{1}, '.zarr']));
            large = single(readzarr([d.dataPath, 'large/', d.fsns{1}, '.zarr']));

            shell = true(testCase.ImSize);
            shell(e + 1 : end - e, e + 1 : end - e, e + 1 : end - e) = false;
            testCase.verifyEqual(nnz(large(shell)), 0, 'The image border is not eroded.');
            b = 8;
            testCase.verifyImagesClose(large(b + 1 : end - b, b + 1 : end - b, b + 1 : end - b), ...
                regular(b + 1 : end - b, b + 1 : end - b, b + 1 : end - b), 'relTol', 0.01, 'absTol', 2);
        end

        function largeFileWithErodeByFirstTimePointCompletes(testCase)
            % regression test for d83468d: with largeFile=true and
            % erodeByFTP=true the first time point's mask was never saved,
            % so the other files waited for it forever. A watchdog writes
            % the mask after 60 s so the old code finishes (and fails here)
            % instead of hanging.
            d = testCase.createDataset({'A'}, 2);
            e = 3;
            watchdogFn = [d.root, 'watchdog_fired'];
            maskFn = sprintf('%smatlab_decon/Masks/%s_eroded.zarr', d.dataPath, d.fsns{1});
            sz = testCase.ImSize;
            watchdog = timer('StartDelay', 60, 'TimerFcn', @(~, ~) testCase.releaseWaitingWrapper(maskFn, sz, watchdogFn));
            testCase.addTeardown(@() delete(watchdog));
            testCase.addTeardown(@() stop(watchdog));
            start(watchdog);

            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'edgeErosion', e, 'erodeByFTP', true, ...
                'largeFile', true, 'batchSize', [64, 64, 64]);

            testCase.verifyFalse(exist(watchdogFn, 'file') == 2, ...
                'The wrapper waited for an erosion mask that the large-file path never writes.');
            testCase.verifyFalse(exist([d.dataPath, 'matlab_decon/Masks'], 'dir') == 7);
            shell = true(sz);
            shell(e + 1 : end - e, e + 1 : end - e, e + 1 : end - e) = false;
            for t = 1 : 2
                out = readtiff([d.dataPath, 'matlab_decon/', d.fsns{t}, '.tif']);
                testCase.verifyEqual(nnz(out(shell)), 0, sprintf('Time point %d is not eroded.', t));
                testCase.verifyGreaterThan(nnz(out(~shell)), 0);
            end
        end

        function widefieldDemoSettingsWithPlainFilenames(testCase)
            % demo_widefield_and_confocal_deconvolution: plain file names, PSF
            % next to the data, OMW with skewed=true and hannWinBounds [0.8, 1.2]
            root = testCase.createTempFolder();
            dataPath = [root, 'Widefield/'];
            mkdir(dataPath);
            pos = [24, 24, 16];
            writetiff(uint16(pktest.beads(testCase.ImSize, pos, testCase.BlurSigma, 1000, testCase.Background)), [dataPath, 'WF_raw.tif']);
            writetiff(uint16(pktest.deconPSF([31, 31, 31], [1.5, 2])), [dataPath, 'WF_PSF.tif']);

            evalc(['XR_decon_data_wrapper({dataPath}, ''resultDirName'', ''matlab_decon_omw'', ''xyPixelSize'', 0.157, ', ...
                '''dz'', 0.1, ''Reverse'', true, ''ChannelPatterns'', {''WF_raw''}, ''PSFFullpaths'', {[dataPath, ''WF_PSF.tif'']}, ', ...
                '''dzPSF'', 0.1, ''parseSettingFile'', false, ''RLmethod'', ''omw'', ''wienerAlpha'', 0.02, ', ...
                '''OTFCumThresh'', 0.65, ''skewed'', true, ''hannWinBounds'', [0.8, 1.2], ''Background'', 100, ', ...
                '''DeconIter'', 3, ''fixIter'', true, ''EdgeErosion'', 0, ''Save16bit'', true, ''zarrFile'', false, ', ...
                '''saveZarr'', false, ''parseCluster'', false, ''largeFile'', false, ''GPUJob'', false, ''debug'', false, ', ...
                '''cpusPerTask'', 4, ''ConfigFile'', '''', ''GPUConfigFile'', '''', ''mccMode'', false)']);

            out = readtiff([dataPath, 'matlab_decon_omw/WF_raw.tif']);
            testCase.verifyEqual(size(out), testCase.ImSize);
            testCase.verifyGreaterThan(max(out(:)), uint16(1.5 * 1000), 'The bead was not sharpened.');
            testCase.verifyEqual(pktest.centroid(single(out), pos, 6), pos, 'AbsTol', 0.1);
        end

        function deconvolvedDataCanBeDeskewedAndRotated(testCase)
            % demo_RL_deconvolution step 5: deskew/rotate raw and decon results
            d = testCase.createDataset({'A'}, 1);
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 3);
            evalc(['XR_deskew_rotate_data_wrapper({d.dataPath, [d.dataPath, ''matlab_decon/'']}, ''xyPixelSize'', 0.108, ', ...
                '''dz'', 0.2, ''Reverse'', true, ''ChannelPatterns'', {''CamA_ch0''}, ''largeFile'', false, ''zarrFile'', false, ', ...
                '''saveZarr'', false, ''Save16bit'', true, ''parseCluster'', false, ''masterCompute'', true)']);
            raw = readtiff([d.dataPath, 'DSR/', d.fsns{1}, '.tif']);
            decon = readtiff([d.dataPath, 'matlab_decon/DSR/', d.fsns{1}, '.tif']);
            testCase.verifyEqual(size(decon), size(raw));
            testCase.verifyGreaterThan(max(decon(:)), max(raw(:)));
        end
    end

    methods (Test, TestTags = {'Integration', 'KnownIssue'})
        function originalMethodWorksWithMultipleBatches(testCase)
            % KNOWN ISSUE: RLdecon.m:387-388 runs deconvlucy for
            % RLMethod='original' but never applies deconBbox (and it is not
            % applied for DeconIter=0 with a non-empty bbox either), so in
            % the in-memory large-file path every batch returns its padded
            % size and RLdecon_large_in_memory fails to place it. Fix: crop
            % with deconBbox after deconvlucy as the other methods do.
            d = testCase.createDataset({'A'}, 1);
            try
                testCase.runWrapper(d, 'RLMethod', 'original', 'deconIter', 3, 'edgeErosion', 0, ...
                    'largeFile', true, 'batchSize', [24, 24, 16]);
            catch ME
                testCase.verifyFail(sprintf('Multi-batch large-file decon with RLMethod=original errors: %s', ME.message));
                return;
            end
            testCase.verifyTrue(exist([d.dataPath, 'matlab_decon/', d.fsns{1}, '.tif'], 'file') == 2, ...
                'Multi-batch large-file decon with RLMethod=original did not produce a result.');
        end
    end

    methods (Access = private)
        function [d] = createDataset(testCase, cams, nT, psfSigmas, zarrFile)
            % light sheet style tiles (uint16, background 100, blurred
            % beads, no zero voxels) in <root>/data/ and one PSF per camera
            % in <root>/psf/
            if nargin < 4 || isempty(psfSigmas)
                psfSigmas = repmat({[1.5, 2]}, 1, numel(cams));
            end
            if nargin < 5
                zarrFile = false;
            end
            d.root = testCase.createTempFolder();
            d.dataPath = [d.root, 'data/'];
            psfPath = [d.root, 'psf/'];
            mkdir(d.dataPath);
            mkdir(psfPath);

            d.fsns = {};
            d.rawFns = {};
            d.psfFns = {};
            d.psfFsns = {};
            d.channelPatterns = {};
            for c = 1 : numel(cams)
                d.channelPatterns{c} = sprintf('Cam%s_ch0', cams{c});
                d.psfFsns{c} = sprintf('PSF_Cam%s', cams{c});
                d.psfFns{c} = [psfPath, d.psfFsns{c}, '.tif'];
                writetiff(uint16(pktest.deconPSF([31, 31, 31], psfSigmas{c})), d.psfFns{c});
                for t = 1 : nT
                    fsn = sprintf('Scan_Iter_0000_0000_Cam%s_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_000x_000y_000z_%04dt', ...
                        cams{c}, t - 1);
                    pos = [8 + 32 * rand(6, 2), 6 + 20 * rand(6, 1)];
                    im = uint16(pktest.beads(testCase.ImSize, pos, testCase.BlurSigma, 1000, testCase.Background));
                    if zarrFile
                        d.rawFns{end + 1} = [d.dataPath, fsn, '.zarr'];
                        writezarr(im, d.rawFns{end}, 'blockSize', [16, 16, 16]);
                    else
                        d.rawFns{end + 1} = [d.dataPath, fsn, '.tif'];
                        writetiff(im, d.rawFns{end});
                    end
                    d.fsns{end + 1} = fsn;
                end
            end
        end

        function [] = runWrapper(testCase, d, varargin)
            % XR_decon_data_wrapper with the dataset settings; varargin
            % adds or overrides parameters
            args = [{'channelPatterns', d.channelPatterns, 'psfFullpaths', d.psfFns, ...
                'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, 'dzPSF', testCase.Dz, ...
                'background', testCase.Background, 'edgeErosion', 0, 'save16bit', true, ...
                'parseCluster', false, 'GPUJob', false}, varargin];
            evalc('XR_decon_data_wrapper({d.dataPath}, args{:})');
        end
    end

    methods (Static, Access = private)
        function [] = releaseWaitingWrapper(maskFn, sz, flagFn)
            % watchdog for the d83468d regression test: write the mask the
            % old code waits for, and record that it was needed
            fclose(fopen(flagFn, 'w'));
            writezarr(uint8(ones(sz)), maskFn);
        end
    end
end
