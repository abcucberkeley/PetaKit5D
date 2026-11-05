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
            expected = decon_lucy_function(readtiff(d.rawFns{1}), psf ./ sum(psf(:), 'double'), 5, ...
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
            expected = decon_lucy_omw_function(readtiff(d.rawFns{1}), psf ./ sum(psf(:), 'double'), readtiff(bpFn), 2, ...
                'Background', testCase.Background, 'useGPU', false, 'save16bit', true);
            testCase.verifyEqual(readtiff([d.dataPath, 'matlab_decon/', d.fsns{1}, '.tif']), expected);
        end

        function save16bitFalseWritesSingle(testCase)
            d = testCase.createDataset({'A'}, 1);
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 3, 'save16bit', false);
            out = readtiff([d.dataPath, 'matlab_decon/', d.fsns{1}, '.tif']);
            testCase.verifyClass(out, 'single');

            psf = single(readtiff([d.dataPath, 'matlab_decon/psfgen/', d.psfFsns{1}, '_simplified.tif']));
            expected = decon_lucy_function(readtiff(d.rawFns{1}), psf ./ sum(psf(:), 'double'), 3, ...
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
                expected = decon_lucy_function(readtiff(d.rawFns{c}), psf ./ sum(psf(:), 'double'), 3, ...
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

        function largeFileMakesMIPBeforePublishing(testCase)
            % the wrapper counts a file as done once its result exists; the
            % in-memory large-file path moved the result into place and then
            % wrote its MIP
            d = testCase.createDataset({'A'}, 1);
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'resultDirName', 'large', ...
                'largeFile', true, 'batchSize', [64, 64, 64]);
            deconFn = [d.dataPath, 'large/', d.fsns{1}, '.tif'];
            testCase.assertTrue(moves.wasMoved(deconFn), 'The result was not moved into place.');
            published = moves.filesWhenMoved(deconFn);
            testCase.verifyTrue(ismember(['MIPs/', d.fsns{1}, '_MIP_z.tif'], published), ...
                'The MIP was written after the result was moved into place.');
            added = setdiff(pktest.MovefileSpyFixture.listFiles([d.dataPath, 'large']), [published; {[d.fsns{1}, '.tif']}]);
            testCase.verifyEmpty(added, 'Files were written after the result was moved into place.');
        end

        function largeFileErosionMaskMatchesDistanceMap(testCase)
            % the in-memory large-file path erodes the edges with a distance
            % map of the data; without zeros that is the volume without
            % EdgeErosion voxels at each face, made without bwdist. Data with
            % zeros (e.g., outside a deskewed region) keeps bwdist.
            d = testCase.createDataset({'A'}, 1);
            im = readtiff(d.rawFns{1});
            % a slanted corner without data, where the distance map differs
            % from the 2D disk and z line erosion used for other masks
            [y, x, z] = ndgrid(1 : size(im, 1), 1 : size(im, 2), 1 : size(im, 3));
            im(y + 2 * x + 3 * z < 60) = 0;
            fsnZeros = strrep(d.fsns{1}, '0000t', '0001t');
            writetiff(im, [d.dataPath, fsnZeros, '.tif']);
            e = 3;
            deconPath = [d.root, 'large/'];
            for fsn = {d.fsns{1}, fsnZeros}
                % SaveMaskfile isn't passed on by the wrapper
                evalc(['RLdecon_large_in_memory([d.dataPath, fsn{1}, ''.tif''], d.psfFns{1}, [deconPath, fsn{1}, ''.tif''], ', ...
                    'testCase.XYPixelSize, testCase.Dz, ''Background'', testCase.Background, ''dzPSF'', testCase.Dz, ', ...
                    '''DeconIter'', 1, ''EdgeErosion'', e, ''SaveMaskfile'', true, ''batchSize'', [64, 64, 64], ', ...
                    '''useGPU'', false, ''parseCluster'', false)']);
                bw = readtiff([d.dataPath, fsn{1}, '.tif']) > 0;
                bw([1, end], :, :) = false;
                bw(:, [1, end], :) = false;
                bw(:, :, [1, end]) = false;
                expected = bwdist(~bw) > e - 1;
                mask = readtiff(sprintf('%sMasks/%s_eroded.tif', deconPath, fsn{1})) > 0;
                testCase.verifyEqual(mask, expected, fsn{1});
            end
        end

        function largeFileInMemoryUsesBlockSize(testCase)
            % the in-memory large-file path wrote its zarr result with
            % writezarr's default chunks (500^3) instead of blockSize
            d = testCase.createDataset({'A'}, 1, {}, true);
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 1, 'zarrFile', true, 'saveZarr', true, ...
                'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inmemory', 'batchSize', [64, 64, 64], ...
                'blockSize', [16, 16, 16]);
            info = getZarrInfo([d.dataPath, 'large/', d.fsns{1}, '.zarr']);
            testCase.verifyEqual(info.blockSize, [16, 16, 16]);
        end

        function largeFileInPlaceMakesMIPsBeforePublishing(testCase)
            % as above for the in-place large-file path, which made its MIPs
            % with XR_MIP_zarr after moving the result into place
            d = testCase.createDataset({'A'}, 1, {}, true);
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'zarrFile', true, 'saveZarr', true, ...
                'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', 'batchSize', [64, 64, 64], ...
                'blockSize', [16, 16, 16]);
            deconFn = [d.dataPath, 'large/', d.fsns{1}, '.zarr'];
            testCase.assertTrue(moves.wasMoved(deconFn), 'The result was not moved into place.');
            published = moves.filesWhenMoved(deconFn);
            mips = strcat('MIPs/', d.fsns{1}, {'_MIP_y.tif'; '_MIP_x.tif'; '_MIP_z.tif'});
            testCase.verifyEqual(ismember(mips, published), true(3, 1), 'MIPs were made after the result was moved into place.');
            added = setdiff(pktest.MovefileSpyFixture.listFiles([d.dataPath, 'large']), [published; {[d.fsns{1}, '.zarr']}]);
            testCase.verifyEmpty(added, 'Files were written after the result was moved into place.');
        end

        function largeFileInPlaceMIPsAreMergedFromTheBlocks(testCase)
            % the decon tasks save the MIPs of their blocks in their flags,
            % and XR_MIP_zarr merges them into the MIP slabs instead of
            % reading the result again in MIP_block tasks. The blocks are
            % single and the result uint16: the MIPs are rounded as written.
            d = testCase.createDataset({'A'}, 1, {}, true);
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            % MIP_block tasks would read the result; the spy errors at once, as the
            % framework would retry a task without output for minutes
            mipTasks = testCase.applyFixture(pktest.FunctionSpyFixture('MIP_block', [], 0));
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'zarrFile', true, 'saveZarr', true, ...
                'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', 'batchSize', [32, 32, 16], ...
                'blockSize', [16, 16, 16]);
            testCase.verifyEmpty(mipTasks.getCalls(), 'The MIPs were computed from the result.');
            calls = moves.getCalls();
            testCase.verifyTrue(any(contains({calls.dst}, '/zarr_flag/') & endsWith({calls.dst}, '.mat')), ...
                'The block flags were not moved into place.');
            pktest.verifyLargeFileMIPs(testCase, [d.dataPath, 'large/', d.fsns{1}, '.zarr']);
        end

        function largeFileInPlaceMIPsIncludeMaskedBlocks(testCase)
            % blocks outside the 2D masks are written as zeros without decon,
            % and their MIPs are saved too (here the blocks below row 32)
            d = testCase.createDataset({'A'}, 1, {}, true);
            sz = testCase.ImSize;
            maskFns = strcat(d.root, {'mask_xy.zarr', 'mask_xz.zarr', 'mask_yz.zarr'});
            xy = ones(sz(1 : 2), 'uint8');
            xy(33 : end, :) = 0;
            writezarr(xy, maskFns{1});
            writezarr(ones(sz(2 : 3), 'uint8'), maskFns{2});
            writezarr(ones(sz([1, 3]), 'uint8'), maskFns{3});
            mipTasks = testCase.applyFixture(pktest.FunctionSpyFixture('MIP_block', [], 0));
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'zarrFile', true, 'saveZarr', true, ...
                'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', 'batchSize', [32, 32, 16], ...
                'blockSize', [16, 16, 16], 'maskFullpaths', maskFns);
            testCase.verifyEmpty(mipTasks.getCalls(), 'The MIPs were computed from the result.');
            resultFn = [d.dataPath, 'large/', d.fsns{1}, '.zarr'];
            out = readzarr(resultFn);
            testCase.assertEqual(nnz(out(33 : end, :, :)), 0, 'The blocks outside the masks were deconvolved.');
            testCase.assertGreaterThan(nnz(out(1 : 32, :, :)), 0);
            pktest.verifyLargeFileMIPs(testCase, resultFn);
        end

        function largeFileInPlacePassesJobSettingsToMIPs(testCase)
            % RLdecon_large_in_place didn't pass parseCluster, masterCompute
            % or jobLogDir to XR_MIP_zarr, so on a machine with Slurm the
            % MIPs of a local run were made in Slurm jobs
            d = testCase.createDataset({'A'}, 1, {}, true);
            spy = testCase.applyFixture(pktest.FunctionSpyFixture('XR_MIP_zarr'));
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'zarrFile', true, 'saveZarr', true, ...
                'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', 'batchSize', [64, 64, 64], ...
                'blockSize', [16, 16, 16]);
            testCase.assertNotEmpty(spy.getCalls(), 'XR_MIP_zarr was not called.');
            testCase.verifyEqual(spy.getOption(1, 'parseCluster'), false);
            testCase.verifyEqual(spy.getOption(1, 'masterCompute'), true);
            testCase.verifyNotEmpty(spy.getOption(1, 'jobLogDir'));
        end

        function largeFileInPlaceGPUJobMIPsKeepMasterCompute(testCase)
            % GPU jobs turn masterCompute off for the decon blocks only; the
            % MIPs, made on CPUs, keep the caller's setting
            d = testCase.createDataset({'A'}, 1, {}, true);
            % the (GPU) decon blocks don't run
            testCase.applyFixture(pktest.FrameworkSpyFixture);
            spy = testCase.applyFixture(pktest.FunctionSpyFixture('XR_MIP_zarr'));
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'zarrFile', true, 'saveZarr', true, ...
                'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', 'batchSize', [64, 64, 64], ...
                'blockSize', [16, 16, 16], 'GPUJob', true);
            testCase.assertNotEmpty(spy.getCalls(), 'XR_MIP_zarr was not called.');
            testCase.verifyEqual(spy.getOption(1, 'masterCompute'), true);
        end

        function largeFileInPlaceKeepsAResultAnotherWriterPublished(testCase)
            % two job sets on the same data (a controller started again, or
            % another call): the one that finished last deleted the result
            % the other had published, which the next step (deskew/rotate)
            % was reading, and moved its own into place. Now it keeps the
            % other's result and removes its own. Here the other result
            % comes while the blocks run (the framework is a spy that
            % reports them done), so this one's result is incomplete.
            d = testCase.createDataset({'A'}, 1, {}, true);
            resultFn = [d.dataPath, 'large/', d.fsns{1}, '.zarr'];
            other = uint16(7 * ones(testCase.ImSize));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture(@() pktest.atCall([d.root, 'calls.txt'], 1, ...
                @() writezarr(other, resultFn, 'blockSize', [16, 16, 16]))));
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'zarrFile', true, 'saveZarr', true, ...
                'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', 'batchSize', [64, 64, 64], ...
                'blockSize', [16, 16, 16]);
            testCase.verifyEqual(readzarr(resultFn), other, 'The other writer''s result was replaced.');
            testCase.verifyEmpty(dir([resultFn, '_*']), 'The temporary result was not removed.');
            testCase.verifyEmpty(dir([d.dataPath, 'large/zarr_flag/', d.fsns{1}, '_*']), 'The flags were not removed.');
            % (no MIPs of this one's incomplete result, named after the result)
            testCase.verifyEmpty(dir([d.dataPath, 'large/MIPs/*.tif']), 'MIPs were made from the incomplete result.');
            % (and the blocks stop once the result is there; on Windows the
            % decon folder ends with filesep)
            testCase.verifyEqual(regexprep(strrep(spy.getOption(1, 'finalOutFullpath'), '\', '/'), '/+', '/'), resultFn);
        end

        function largeFileInPlaceBlocksOfAnEndedControllerAreSkipped(testCase)
            % the block jobs of a controller whose job was cancelled ran on
            % (the queued ones when they started), for a result nothing
            % would publish; a block now skips its work when its
            % controller's job has ended (its writer record; a fake Slurm
            % that no longer lists the job). The block task is the one the
            % controller made (a framework spy that reports it failed).
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 1, {}, true);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture([], false));
            testCase.verifyError(@() testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'zarrFile', true, ...
                'saveZarr', true, 'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', ...
                'batchSize', [64, 64, 64], 'blockSize', [16, 16, 16]), ?MException);
            calls = spy.getCalls();
            flagFn = calls{1}.outputFullpaths{1};
            recordFn = [fileparts(flagFn), '/writer.txt'];
            testCase.assertTrue(isfile(recordFn), 'The controller did not record itself.');
            % (it ended by an error, and says so: its process may stay, e.g.
            % an interactive session, and no one may wait for it)
            testCase.verifySubstring(fileread(recordFn), sprintf('\nended=1'), ...
                'The controller that failed did not mark its record ended.');
            % a controller whose job was cancelled can't say so: its job is gone
            txt = regexprep(fileread(recordFn), {'owner=[^\n]*', 'ended=[^\n]*\n'}, {'owner=777_3:othernode:5', ''});
            fid = fopen(recordFn, 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
            testCase.applyFixture(pktest.FakeSlurmFixture);
            % (squeue's answers are kept per job for a minute; fake job ids repeat)
            clear('result_writers');
            funcStr = calls{1}.funcStrs{1};
            evalc('feval(str2func([''@()'', funcStr]))');
            testCase.verifyFalse(isfile(flagFn), 'The block of an ended controller was deconvolved.');
        end

        function largeFileInPlaceKeepsAResultPublishedDuringItsMIPs(testCase)
            % as above, with the other result coming while the MIPs are made
            % (a spy of XR_MIP_zarr makes it)
            d = testCase.createDataset({'A'}, 1, {}, true);
            resultFn = [d.dataPath, 'large/', d.fsns{1}, '.zarr'];
            other = uint16(7 * ones(testCase.ImSize));
            testCase.applyFixture(pktest.FunctionSpyFixture('XR_MIP_zarr', @() writezarr(other, resultFn, ...
                'blockSize', [16, 16, 16])));
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'zarrFile', true, 'saveZarr', true, ...
                'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', 'batchSize', [32, 32, 16], ...
                'blockSize', [16, 16, 16]);
            testCase.verifyEqual(readzarr(resultFn), other, 'The other writer''s result was replaced.');
            testCase.verifyEmpty(dir([resultFn, '_*']), 'The temporary result was not removed.');
            testCase.verifyEmpty(dir([resultFn, '/', d.fsns{1}, '.zarr_*']), 'The result was moved into the other''s.');
        end

        function largeFileInPlaceWaitsForAControllerOfAnotherCall(testCase)
            % a controller of another call (another uuid) that may still
            % make the volume, e.g. of a call cancelled and run again: the
            % second controller made it as well (and replaced the first
            % one's). Now it waits for it. Here the other controller is a
            % process on this node, and its result comes during the wait (a
            % spy of wait_for_files); the blocks must not run.
            testCase.assumeTrue(isunix && ~ismac, 'Writers outside Slurm are checked on Linux only.');
            d = testCase.createDataset({'A'}, 1, {}, true);
            resultFn = [d.dataPath, 'large/', d.fsns{1}, '.zarr'];
            pktest.writeWriterRecord(testCase, resultFn, 'otheruuid');
            other = uint16(7 * ones(testCase.ImSize));
            testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', @() pktest.atCall([d.root, 'calls.txt'], 1, ...
                @() writezarr(other, resultFn, 'blockSize', [16, 16, 16])), 10));
            blocks = testCase.applyFixture(pktest.FunctionSpyFixture('RLdecon_for_zarr_block', [], 0));
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'zarrFile', true, 'saveZarr', true, ...
                'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', 'batchSize', [64, 64, 64], ...
                'blockSize', [16, 16, 16]);
            testCase.verifyEmpty(blocks.getCalls(), 'The volume was deconvolved a second time.');
            testCase.verifyEqual(readzarr(resultFn), other, 'The other controller''s result was replaced.');
            testCase.verifyEmpty(dir([resultFn, '_*']), 'A temporary result was left.');
        end

        function clusterModeWaitsForAControllerOfAnotherCall(testCase)
            % the wrapper submitted a controller for a volume that a
            % controller of another call (a job in the queue) still made;
            % it now waits for that job, without a try. Here the other job's
            % result comes at the second wait.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 1, {}, true);
            resultFn = [d.dataPath, 'large/', d.fsns{1}, '.zarr'];
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('777_1', {'R'});
            pktest.writeWriterRecord(testCase, resultFn, 'otheruuid', '777_1:othernode:5');
            testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', @() pktest.atCall([d.root, 'calls.txt'], 2, ...
                @() mkdir(resultFn)), 20));
            testCase.runWrapper(d, 'RLMethod', 'simplified', 'zarrFile', true, 'saveZarr', true, 'resultDirName', 'large', ...
                'largeFile', true, 'largeMethod', 'inplace', 'batchSize', [64, 64, 64], 'blockSize', [16, 16, 16], ...
                'parseCluster', true, 'maxTrialNum', 1, 'jobLogDir', [d.root, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 0, 'A controller was submitted while another call''s made the volume.');
        end

        function largeFileInPlaceControllersAskForWhatTheyHold(testCase)
            % the job of a volume (RLdecon_large_in_place) runs its block
            % jobs. It asked for a block-based estimate it never holds:
            % prod(batchSize * 2) x 40 B on the CPU (16 CPUs here) and
            % prod(batchSize) x 80 B with GPUs (4), though GPU controllers
            % compute no block (1.5-1.7 GB of 86 GB were used), and
            % with masterCompute off (the config's, as in a user's, or the
            % wrapper's) no CPU controller does: 4 of them held 4 nodes
            % while their blocks pended. Now a controller asks for MATLAB
            % and its MIP merge with one CPU, and for the CPUs of a block's
            % job if it deconvolves blocks itself. The volume is a 34
            % G-voxel header (the jobs don't run).
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 1, {}, true);
            rmdir(d.rawFns{1}, 's');
            createzarr(d.rawFns{1}, dataSize=[4096, 4096, 2048], blockSize=[256, 256, 256], dtype='uint16');
            offConfig = [d.root, 'cpu_config.json'];
            fid = fopen(offConfig, 'w');
            fprintf(fid, '{"masterCompute": false}');
            fclose(fid);
            % GPUJob, masterCompute, configFile, cpusPerTask, controller CPUs
            cases = {false, true, '', 4, 4; false, true, offConfig, 4, 1; false, false, '', 4, 1; true, true, '', 2, 1};
            testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', [], 40));
            % (the controllers wait for their block jobs: these go first)
            testCase.applyFixture(pktest.FunctionStubFixture('coordinator_nice', 1000));
            for i = 1 : size(cases, 1)
                slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
                testCase.runWrapper(d, 'RLMethod', 'simplified', 'zarrFile', true, 'saveZarr', true, ...
                    'resultDirName', sprintf('decon_%d', i), 'largeFile', true, 'largeMethod', 'inplace', ...
                    'batchSize', [1024, 1024, 1024], 'blockSize', [256, 256, 256], 'GPUJob', cases{i, 1}, ...
                    'masterCompute', cases{i, 2}, 'configFile', cases{i, 3}, 'cpusPerTask', cases{i, 4}, ...
                    'parseCluster', true, 'maxTrialNum', 1, 'jobLogDir', [d.root, 'logs']);
                lines = testCase.sbatchLines(slurm);
                testCase.assertNumElements(lines, 1, sprintf('Case %d.', i));
                cpus = str2double(regexp(lines{1}, '--cpus-per-task=(\d+)', 'tokens', 'once'));
                testCase.verifyEqual(cpus, cases{i, 5}, sprintf('Case %d: the controller''s CPUs.', i));
                % (the wrapper's masterCompute reaches the controller)
                testCase.verifyEqual(contains(lines{1}, '''masterCompute'',false'), ~cases{i, 2}, sprintf('Case %d.', i));
                testCase.verifySubstring(lines{1}, ' --nice=1000', sprintf('Case %d: the controller went before its blocks.', i));
            end
        end

        function largeFileInPlaceBlocksLogToTheCallersJobLogDir(testCase)
            % the controller of a volume got no job log folder, so its block
            % jobs logged to its default (../job_logs, or the data folder's
            % job_logs where that doesn't exist) instead of the caller's
            % jobLogDir, where the wrapper's own jobs log. Here the
            % controller runs here and its block framework is a spy.
            d = testCase.createDataset({'A'}, 1, {}, true);
            logDir = [d.root, 'caller_logs'];
            mkdir(logDir);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture([], false));
            testCase.verifyError(@() testCase.runWrapper(d, 'RLMethod', 'simplified', 'deconIter', 2, 'zarrFile', true, ...
                'saveZarr', true, 'resultDirName', 'large', 'largeFile', true, 'largeMethod', 'inplace', ...
                'batchSize', [64, 64, 64], 'blockSize', [16, 16, 16], 'jobLogDir', logDir), ?MException);
            testCase.assertNotEmpty(spy.getCalls(), 'The controller ran no block framework.');
            testCase.verifyEqual(spy.getOption(1, 'jobLogDir'), logDir);
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

        function originalMethodWorksWithMultipleBatches(testCase)
            % 'original' used to skip the deconBbox crop, so the batches of the
            % in-memory large-file path couldn't be placed
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

        function originalMethodBatchesAreCopiedByTheMex(testCase)
            % RLMethod 'original' returns single batches, which the in-memory
            % large-file path placed by assignment (converting them to uint16)
            % instead of with indexing3d_mex, which takes the class of the
            % result: they are converted first now, with the same result
            d = testCase.createDataset({'A'}, 1);
            deconPath = [d.root, 'large/'];
            outFn = [deconPath, d.fsns{1}, '.tif'];
            % batches of 24 x 24 x 16 of the 48 x 48 x 32 image: 8
            args = {d.rawFns{1}, d.psfFns{1}, outFn, testCase.XYPixelSize, testCase.Dz, 'Background', testCase.Background, ...
                'dzPSF', testCase.Dz, 'DeconIter', 3, 'EdgeErosion', 0, 'RLMethod', 'original', 'batchSize', [24, 24, 16], ...
                'useGPU', false, 'parseCluster', false};
            evalc('RLdecon_large_in_memory(args{:})');
            out = readtiff(outFn);
            delete(outFn);
            testCase.verifyClass(out, 'uint16');
            testCase.verifyGreaterThan(max(out(:)), 0);

            % the same image as the assignment (the mex hidden by a stub that errors)
            testCase.applyFixture(pktest.FunctionStubFixture('indexing3d_mex'));
            testCase.expectMexFallback('indexing3d_mex');
            evalc('RLdecon_large_in_memory(args{:})');
            testCase.verifyEqual(readtiff(outFn), out);
            delete(outFn);

            % one copy per batch, of the output class (a spy, first on the path)
            spy = testCase.applyFixture(pktest.FunctionSpyFixture('indexing3d_mex'));
            evalc('RLdecon_large_in_memory(args{:})');
            calls = spy.getCalls();
            testCase.verifyNumElements(calls, 8);
            testCase.verifyTrue(all(cellfun(@(c) isa(c{3}, 'uint16'), calls)), 'A batch was not converted to uint16.');
        end

        function clusterModeChecksFewJobsSoon(testCase)
            % with jobs on the cluster the wrapper always waited 30 s between
            % rounds, so the last results were seen up to 30 s late; it now
            % waits 0.3 s per unfinished file, at least 3 s and at most 30 s
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 1);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            waitSpy = testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', [], 20));
            % the fake jobs never run: one round submits, the next resubmits,
            % the last checks that the second job is gone
            testCase.runWrapper(d, 'parseCluster', true, 'maxTrialNum', 2, 'jobLogDir', [d.root, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 2);
            waits = cellfun(@(x) x{2}, waitSpy.getCalls());
            testCase.verifyEqual(waits, [3, 3]);
        end

        function clusterModeEndsSoonAfterTheLastResult(testCase)
            % the wrapper waited the whole wait between its rounds (0.3 s per
            % unfinished file, at least 3 s), so it saw the last result up
            % to that late; it now looks for the results still to come
            % each second and goes on once they are all there. The fake
            % job doesn't run; its result comes during the second second.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 1);
            outFn = [d.dataPath, 'matlab_decon/', d.fsns{1}, '.tif'];
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [d.root, 'calls.txt'], 2, @() fclose(fopen(outFn, 'w'))), 20));
            testCase.runWrapper(d, 'parseCluster', true, 'maxTrialNum', 2, 'jobLogDir', [d.root, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 1, 'The job was submitted again.');
            testCase.verifyEqual(cellfun(@(x) x{1}, pauses.getCalls()), [1, 1], ...
                'The wrapper did not go on within a second of the result.');
        end

        function clusterModeSubmitsFilesSoonAfterTheirMask(testCase)
            % (erodeByFTP) the other files of a dataset wait for the erosion
            % mask of its first time point; the wrapper now looks for it
            % each second and submits them once it is there. The first
            % file's job makes its result when first checked (its mask
            % comes during the second second); the second file's job
            % doesn't run.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 2);
            maskFn = sprintf('%smatlab_decon/Masks/%s_eroded.zarr', d.dataPath, d.fsns{1});
            mkdir(fileparts(maskFn));
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('1000_1', {'R', 'done'}, {[d.dataPath, 'matlab_decon/', d.fsns{1}, '.tif']});
            sz = testCase.ImSize;
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [d.root, 'calls.txt'], 2, @() writezarr(uint8(ones(sz)), maskFn)), 20));
            testCase.runWrapper(d, 'edgeErosion', 2, 'parseCluster', true, 'maxTrialNum', 1, ...
                'jobLogDir', [d.root, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 2);
            % two seconds until the mask, then a whole wait for the second
            % file's job
            testCase.verifyEqual(cellfun(@(x) x{1}, pauses.getCalls()), ones(1, 5), ...
                'The wrapper did not submit the second file within a second of the mask.');
        end

        function clusterModeSeesMasksOtherNodesJustMade(testCase)
            % (erodeByFTP) on NFS, exist() can miss a mask another node just
            % made for 30 s or more, while the wait sees it at once
            % (batch_file_exist), so the wrapper checks masks as the wait
            % does: a wait that ends at a mask is followed by the
            % submission of the files that waited for it, not by rounds a
            % second apart. The stub reports the first file's mask made
            % (only batch_file_exist sees it): the second file is submitted
            % in the first round, and the wrapper waits once.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 2);
            stubFolder = testCase.createTempFolder();
            fid = fopen([stubFolder, 'batch_file_exist.m'], 'w');
            fprintf(fid, ['function [tf] = batch_file_exist(paths, varargin)\n', ...
                'tf = endsWith(paths(:), ''_eroded.zarr'') | cellfun(@(x) isfile(x) || isfolder(x), paths(:));\n', ...
                'end\n']);
            fclose(fid);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(stubFolder));
            testCase.applyFixture(pktest.FunctionStubFixture('getImageSize', testCase.ImSize));
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 10));
            testCase.runWrapper(d, 'edgeErosion', 2, 'parseCluster', true, 'maxTrialNum', 1, ...
                'jobLogDir', [d.root, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 2);
            testCase.verifyEqual(cellfun(@(x) x{1}, pauses.getCalls()), ones(1, 3), ...
                'The second file was not submitted in the first round.');
        end

        function clusterModeWaitsForTheLastTry(testCase)
            % after submitting the last try, the wrapper returned without
            % waiting for its job; it now checks the job until it ends,
            % without submitting again
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 1);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('1000_1', {'PD', 'R', 'R', 'done'});
            waitSpy = testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', [], 20));
            testCase.runWrapper(d, 'parseCluster', true, 'maxTrialNum', 1, 'jobLogDir', [d.root, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 1, 'The last try was submitted again.');
            % one check per round until the job is gone (none before
            % submitting: the task had no job)
            testCase.verifyEqual(slurm.callCount('squeue'), 4);
            testCase.verifyNumElements(waitSpy.getCalls(), 4);
        end

        function clusterModeWaitsForTheLastTryOfAJobOfSeveralFiles(testCase)
            % files that share a job also wait for it on their last try,
            % without submitting it again
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 2);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            % the job of both files has the task id of the first
            slurm.setJobStates('1000_1', {'PD', 'R', 'R', 'done'});
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            testCase.runWrapper(d, 'parseCluster', true, 'maxTrialNum', 1, 'taskBatchNum', 2, ...
                'jobLogDir', [d.root, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 1, 'The last try was submitted again.');
            checks = slurm.history();
            checks = checks(startsWith(checks, 'squeue') & contains(checks, '1000_1'));
            testCase.verifyNumElements(checks, 4, 'The job was not checked until it ended.');
        end

        function clusterModeWaitGrowsWithUnfinishedFiles(testCase)
            % the wait between rounds is 0.3 s per unfinished file: 6 s for
            % 20 files (between the 3 s minimum and the 30 s maximum)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 20);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            waitSpy = testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', [], 20));
            % the fake jobs never run: one round submits, the next resubmits
            testCase.runWrapper(d, 'parseCluster', true, 'maxTrialNum', 2, 'jobLogDir', [d.root, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 40);
            waits = cellfun(@(x) x{2}, waitSpy.getCalls());
            % the second round also waits for the last tries' jobs to end
            testCase.verifyEqual(waits, [6, 6], 'AbsTol', 1e-9, 'The wait for 20 unfinished files should be 6 s.');
        end

        function lastTryWaitsForATaskRunnerThatHasTheTask(testCase)
            % (task runners) file 2's job makes its result and then runs
            % file 3's task (a task runner; file 1, the first of its kind,
            % runs as before), so file 3's own job ends at once: the
            % wrapper counts file 3 as running while file 2's job runs,
            % instead of as finished, and does not submit it again
            testCase.assumeTrue(isunix && ~ismac, 'Task runners run on Linux only.');
            d = testCase.createDataset({'A'}, 3);
            logPath = [d.root, 'logs/'];
            mkdir(logPath);
            outFns = strcat(d.dataPath, 'matlab_decon/', d.fsns, '.tif');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('1000_1', {'R', 'done'}, outFns(1));
            slurm.setJobStates('1001_2', {'R', 'R', 'R', 'R', 'done'}, outFns(2));
            % file 3's task is the second of the pool
            claim = testCase.claimCommand(logPath, 2, '1001_2');
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() system(claim), 40));
            testCase.runWrapper(d, 'parseCluster', true, 'maxTrialNum', 1, 'jobLogDir', logPath);
            testCase.verifyEqual(slurm.callCount('sbatch'), 3, 'File 3 was submitted again.');
            checks = slurm.history();
            checks = checks(startsWith(checks, 'squeue') & contains(checks, '1001_2'));
            testCase.verifyGreaterThanOrEqual(numel(checks), 5, 'The wrapper did not wait for the job that had file 3.');
            % (its jobs gone, the master removed the pool at the end)
            testCase.verifyEmpty(dir([logPath, 'task_runner_single_*']), 'The pool is still there.');
        end

        function filesOfAJobWaitForATaskRunnerThatHasThem(testCase)
            % (task runners) GPU jobs of two files each: the third job's
            % task is run by a task runner of the second job (the first
            % job, the first of its kind, runs as before), so the third job
            % ends at once; its files are submitted again (one job each, as
            % before, without task runners) only once the second job ended
            % without their results, not while it may still make them
            testCase.assumeTrue(isunix && ~ismac, 'Task runners run on Linux only.');
            d = testCase.createDataset({'A'}, 6);
            logPath = [d.root, 'logs/'];
            mkdir(logPath);
            outFns = strcat(d.dataPath, 'matlab_decon/', d.fsns, '.tif');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            % the jobs have the task ids of their first files: 1, 3, 5
            slurm.setJobStates('1000_1', {'R', 'done'}, outFns(1 : 2));
            slurm.setJobStates('1001_3', {'R', 'R', 'R', 'R', 'done'}, outFns(3 : 4));
            claim = testCase.claimCommand(logPath, 2, '1001_3');
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() system(claim), 40));
            testCase.runWrapper(d, 'parseCluster', true, 'GPUJob', true, 'taskBatchNum', 2, 'maxTrialNum', 2, ...
                'GPUConfigFile', testCase.writeGPUConfig(d.root, struct()), 'jobLogDir', logPath);
            cmds = slurm.history();
            submits = find(startsWith(cmds, 'sbatch'));
            testCase.assertNumElements(submits, 5);
            running = find(startsWith(cmds, 'squeue') & contains(cmds, '1001_3'));
            testCase.assertGreaterThanOrEqual(numel(running), 5);
            testCase.verifyGreaterThan(submits(4), running(5), 'Files were submitted again while the job that had them ran.');
            testCase.verifyFalse(contains(cmds{submits(1)}, 'run_tasks.sh'), 'The first job of its kind has task runners.');
            testCase.verifyTrue(all(contains(cmds(submits(2 : 3)), 'run_tasks.sh matlab')), 'The pool jobs have no task runners.');
            testCase.verifyFalse(any(contains(cmds(submits(4 : 5)), 'run_tasks.sh')), 'Files were submitted again with task runners.');
        end

        function gpuJobsTakeAtMostFourCPUs(testCase)
            % GPU jobs got the CPU-oriented cpusPerTask (24 in the demos),
            % capped only by the GPU config's maxCPUNum (16), so one GPU job
            % took a whole 16-core 4-GPU node and 3 of its GPUs idled. The GPU
            % config's minCPUNum still raises it (nodes with more cores per
            % GPU), and CPU jobs keep cpusPerTask.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 2);
            cases = {struct('minCPUNum', 1), true, 4; struct('minCPUNum', 8), true, 8; struct(), false, 24};
            for i = 1 : size(cases, 1)
                slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
                config = testCase.writeGPUConfig(d.root, cases{i, 1});
                testCase.runWrapper(d, 'resultDirName', sprintf('decon_%d', i), 'parseCluster', true, ...
                    'GPUJob', cases{i, 2}, 'GPUConfigFile', config, 'cpusPerTask', 24, 'maxTrialNum', 1, ...
                    'jobLogDir', [d.root, 'logs']);
                cpus = cellfun(@(x) str2double(regexp(x, '--cpus-per-task=(\d+)', 'tokens', 'once')), testCase.sbatchLines(slurm));
                testCase.verifyEqual(cpus, repmat(cases{i, 3}, 2, 1), sprintf('Case %d (GPUJob %d).', i, cases{i, 2}));
            end
        end

        function gpuJobsTakeSeveralFilesWithTheSameSettings(testCase)
            % a job per file starts MATLAB (or the MATLAB Runtime) and the
            % GPU for about as long as a GPU decon of a tile takes; with
            % taskBatchNum, files with the same settings (here the channel)
            % share a job, and the jobs give the same results as one file
            % at a time
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A', 'B'}, 3, {[1.5, 2], [2.5, 3]});
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            args = {'RLMethod', 'simplified', 'deconIter', 3, 'GPUJob', true};
            testCase.runWrapper(d, args{:}, 'resultDirName', 'decon_jobs', 'parseCluster', true, ...
                'GPUConfigFile', testCase.writeGPUConfig(d.root, struct()), 'taskBatchNum', 2, ...
                'maxTrialNum', 1, 'jobLogDir', [d.root, 'logs']);
            funcStrs = testCase.sbatchFuncStrs(slurm);
            testCase.assertEqual(numel(funcStrs), 4, 'Expected 2 jobs per channel (2 + 1 files).');
            jobFiles = cellfun(@(x) regexp(x, '/data/(\w+)\.tif''', 'tokens'), funcStrs, 'unif', 0);
            jobFiles = cellfun(@(x) [x{:}], jobFiles, 'unif', 0);
            testCase.verifyEqual(jobFiles, {d.fsns(1 : 2); d.fsns(3); d.fsns(4 : 5); d.fsns(6)});

            % run the jobs here and compare with one file at a time
            for i = 1 : numel(funcStrs)
                evalc('feval(str2func([''@()'', funcStrs{i}]))');
            end
            testCase.runWrapper(d, args{:}, 'resultDirName', 'decon_local');
            for f = 1 : numel(d.fsns)
                testCase.verifyEqual(readtiff(sprintf('%sdecon_jobs/%s.tif', d.dataPath, d.fsns{f})), ...
                    readtiff(sprintf('%sdecon_local/%s.tif', d.dataPath, d.fsns{f})), d.fsns{f});
            end
        end

        function runningJobIsNotSubmittedAgain(testCase)
            % the job of two files runs for two rounds and then ends without
            % results, so its files are submitted again, one per job
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 2);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.applyFixture(pktest.FunctionSpyFixture('pause'));
            % the first job is 1000, array task 1 (the task id of its first file)
            slurm.setJobStates('1000_1', {'R', 'R', 'done'});
            testCase.runWrapper(d, 'parseCluster', true, 'GPUJob', true, 'taskBatchNum', 2, 'maxTrialNum', 2, ...
                'GPUConfigFile', testCase.writeGPUConfig(d.root, struct()), 'jobLogDir', [d.root, 'logs']);
            jobFiles = cellfun(@(x) regexp(x, '/data/(\w+)\.tif''', 'tokens'), testCase.sbatchFuncStrs(slurm), 'unif', 0);
            jobFiles = cellfun(@(x) [x{:}], jobFiles, 'unif', 0);
            testCase.verifyEqual(jobFiles, {d.fsns; d.fsns(1); d.fsns(2)});
            % one status check per round for the shared job
            testCase.verifyEqual(nnz(contains(slurm.calls('squeue'), '1000_1')), 3);
        end

        function retriedFilesGetAJobEach(testCase)
            % a file that fails would make its whole job fail again and use
            % up the trials of the other files, so files are retried alone.
            % The fake jobs never run: the second round retries both files.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 2);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.applyFixture(pktest.FunctionSpyFixture('pause'));
            testCase.runWrapper(d, 'parseCluster', true, 'GPUJob', true, 'taskBatchNum', 2, 'maxTrialNum', 2, ...
                'GPUConfigFile', testCase.writeGPUConfig(d.root, struct()), 'jobLogDir', [d.root, 'logs']);
            jobFiles = cellfun(@(x) regexp(x, '/data/(\w+)\.tif''', 'tokens'), testCase.sbatchFuncStrs(slurm), 'unif', 0);
            jobFiles = cellfun(@(x) [x{:}], jobFiles, 'unif', 0);
            testCase.verifyEqual(jobFiles, {d.fsns; d.fsns(1); d.fsns(2)});
        end

        function filesOfDifferentSizesGoToDifferentJobs(testCase)
            % a job keeps the OTFs of its first file on the GPU, so a file
            % of another size would need memory for both
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 3);
            im = readtiff(d.rawFns{2});
            writetiff(im(:, :, 1 : end - 4), d.rawFns{2});
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.runWrapper(d, 'parseCluster', true, 'GPUJob', true, 'taskBatchNum', 3, 'maxTrialNum', 1, ...
                'GPUConfigFile', testCase.writeGPUConfig(d.root, struct()), 'jobLogDir', [d.root, 'logs']);
            jobFiles = cellfun(@(x) regexp(x, '/data/(\w+)\.tif''', 'tokens'), testCase.sbatchFuncStrs(slurm), 'unif', 0);
            jobFiles = cellfun(@(x) [x{:}], jobFiles, 'unif', 0);
            testCase.verifyEqual(jobFiles, {d.fsns([1, 3]); d.fsns(2)});
        end

        function gpuJobsTakeUpToEightFilesWhileThereAreSixteenJobs(testCase)
            % without taskBatchNum, GPU jobs take floor(files / 16) files,
            % at most 8, so small datasets keep one file per job
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 40);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.runWrapper(d, 'parseCluster', true, 'GPUJob', true, 'maxTrialNum', 1, ...
                'GPUConfigFile', testCase.writeGPUConfig(d.root, struct()), 'jobLogDir', [d.root, 'logs']);
            funcStrs = testCase.sbatchFuncStrs(slurm);
            nFiles = cellfun(@(x) numel(regexp(x, '\.tif''', 'match')) - 1, funcStrs);
            testCase.verifyEqual(nFiles, 2 * ones(20, 1), 'Expected 20 jobs of 2 files for 40 files.');
        end

        function omwBackProjectorIsMadeOnceBeforeTheJobs(testCase)
            % on the cluster every first job of a PSF made the back projector
            % (with the job's threads, which change it in the last bits); the
            % wrapper now makes it once, as RLdecon would from the PSF
            % with two datasets, made once and copied to the second one
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            d = testCase.createDataset({'A'}, 2);
            d2 = testCase.createDataset({'A'}, 1);
            d.dataPaths = {d.dataPath, d2.dataPath};
            testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.runWrapper(d, 'RLMethod', 'omw', 'deconIter', 2, ...
                'wienerAlpha', 0.005, 'OTFCumThresh', 0.9, 'skewed', false, 'parseCluster', true, 'GPUJob', true, ...
                'maxTrialNum', 1, 'GPUConfigFile', testCase.writeGPUConfig(d.root, struct()), 'jobLogDir', [d.root, 'logs']);
            for dataPath = {d.dataPath, d2.dataPath}
                psfgenPath = [dataPath{1}, 'matlab_decon/psfgen/', d.psfFsns{1}];
                bpFn = [psfgenPath, '_back_projector_alpha_0.005000_otf_cum_thresh_0.900000.tif'];
                testCase.assertTrue(exist(bpFn, 'file') == 2, 'The wrapper did not make the back projector.');
                psf = single(readtiff([psfgenPath, '_omw.tif']));
                [~, expected] = evalc('omw_backprojector_generation(psf, 0.005, false, ''OTFCumThresh'', 0.9, ''hannWinBounds'', [0.8, 1])');
                testCase.verifyEqual(readtiff(bpFn), expected);
                testCase.verifyTrue(exist([psfgenPath, '_OTF.tif'], 'file') == 2);
                testCase.verifyTrue(exist([psfgenPath, '_OTF_mask_otf_cum_thresh_0.900000.tif'], 'file') == 2);
            end
        end

        function compiledJobGetsAllFilesOfItsBatch(testCase)
            % in mcc mode the job's command line goes through the MATLAB
            % Runtime parser, which must turn the file list back into a cell
            testCase.assumeFalse(ispc, 'The fake Slurm and the argument split need a POSIX shell.');
            d = testCase.createDataset({'A'}, 2);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.runWrapper(d, 'parseCluster', true, 'GPUJob', true, 'mccMode', true, 'taskBatchNum', 2, ...
                'maxTrialNum', 1, 'GPUConfigFile', testCase.writeGPUConfig(d.root, struct()), 'jobLogDir', [d.root, 'logs']);
            log = fileread([slurm.Folder, '/sbatch.log']);
            argStr = regexp(log, 'XR_RLdeconFrame3D ([^\n]*)', 'tokens', 'once');
            testCase.assertNotEmpty(argStr, 'No compiled XR_RLdeconFrame3D command was submitted.');
            % split the arguments as the job's shell does, in a script run
            % with bash: system() uses the login shell, and zsh (the macOS
            % default) rejects unquoted [...] that match no file
            splitFn = [d.root, 'split_args.sh'];
            fid = fopen(splitFn, 'w');
            fprintf(fid, 'printf ''[%%s]\\n'' %s\n', argStr{1});
            fclose(fid);
            [status, out] = system(sprintf('bash %s', splitFn));
            testCase.assertEqual(status, 0);
            mccArgs = regexp(strip(out), '^\[(.*)\]$', 'tokens', 'lineanchors', 'dotexceptnewline');
            mccArgs = cellfun(@(x) x{1}, mccArgs, 'unif', 0);
            spy = testCase.applyFixture(pktest.FunctionSpyFixture('XR_RLdeconFrame3D'));
            XR_RLdeconFrame3D_parser(mccArgs{:});
            calls = spy.getCalls();
            testCase.assertNumElements(calls, 1);
            testCase.verifyEqual(calls{1}{1}, strcat(d.dataPath, d.fsns, '.tif'));
            testCase.verifyEqual(spy.getOption(1, 'GPUJob'), true);
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
            % XR_decon_data_wrapper with the dataset settings (on
            % d.dataPaths if set); varargin adds or overrides parameters
            args = [{'channelPatterns', d.channelPatterns, 'psfFullpaths', d.psfFns, ...
                'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, 'dzPSF', testCase.Dz, ...
                'background', testCase.Background, 'edgeErosion', 0, 'save16bit', true, ...
                'parseCluster', false, 'GPUJob', false}, varargin];
            dataPaths = {d.dataPath};
            if isfield(d, 'dataPaths')
                dataPaths = d.dataPaths;
            end
            evalc('XR_decon_data_wrapper(dataPaths, args{:})');
        end

        function [fn] = writeGPUConfig(~, folder, overrides)
            % a GPU job config like the production one for 16-core nodes
            % with 4 GPUs; overrides replace or add keys. Each config gets a
            % new file name, because the framework caches a config by name
            % and modification time (1 s resolution).
            config = struct('SlurmParam', '-p gpu_test -n1 --mem-per-cpu=32128M --gres=gpu:1', ...
                'MemPerCPU', 31.375, 'maxCPUNum', 16, 'minCPUNum', 4);
            names = fieldnames(overrides);
            for i = 1 : numel(names)
                config.(names{i}) = overrides.(names{i});
            end
            [~, uuid] = fileparts(tempname);
            fn = sprintf('%sgpu_config_%s.json', folder, uuid);
            fid = fopen(fn, 'w');
            fprintf(fid, '%s', jsonencode(config));
            fclose(fid);
        end

        function [cmd] = claimCommand(~, logPath, i, job)
            % a shell command that has job's task runner claim task i of the
            % single jobs' pool in logPath, once the pool exists
            cmd = sprintf(['for d in "%s"task_runner_single_*/claims; do [ -d "$d" ] && ', ...
                'ln -s "%s:othernode:5:1" "$d/t%d" 2> /dev/null; done; true'], logPath, job, i);
        end

        function [lines] = sbatchLines(~, slurm)
            % the arguments of each sbatch call of MATLAB jobs (one line each)
            lines = slurm.calls('sbatch');
        end

        function [funcStrs] = sbatchFuncStrs(testCase, slurm)
            % the function string of each submitted MATLAB job, in order
            funcStrs = regexp(testCase.sbatchLines(slurm), 'Matlab command:\s+"(XR_RLdeconFrame3D\(.+?\))";', 'tokens', 'once');
            funcStrs = cellfun(@(x) x{1}, funcStrs, 'unif', 0);
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
