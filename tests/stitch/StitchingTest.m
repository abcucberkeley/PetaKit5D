classdef StitchingTest < pktest.TestCase
    % XR_matlab_stitching_wrapper on tiles cut from synthetic ground truth
    % volumes at known positions. With exact positions (xcorrShift false)
    % the stitched volume must reproduce the ground truth voxel by voxel.
    %
    % Canvas convention: the minimum tile corner lands at voxel 1 and the
    % canvas has 2 extra (zero) voxels at the end of every axis. With
    % xcorrShift, 2 * #tiles (per axis) voxels of padding are added on both
    % sides. 0 means "no data" in stitching, so ground truths are nonzero.

    properties (Constant)
        XYPixelSize = 0.108
        StitchName = 'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs'
    end

    properties (TestParameter)
        % minutes, -1 for no marker
        markerAge = struct('none', -1, 'old', 11)
        blendMethod = {'none', 'feather', 'max', 'mean', 'median'}
        tileBlendMethod = {'none', 'max', 'mean', 'median'}
        axisOrder = struct('x', 'x,y,z', 'minusX', '-x,y,z')
        thinTileDepth = struct('nz8', 8, 'nz20', 20)
        % DSR, or skewed space (where flipped tiles are also flipped in z)
        flipDSR = struct('DSR', true, 'skewedSpace', false)
        % volumes the wrapper sizes jobs for (see the tests of the same layouts)
        sizeSetup = struct( ...
            'grid', struct('gtSize', [100, 120, 32], 'starts', [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], ...
                'tileSize', [64, 72, 32], 'params', {{}}), ...
            'gridPadded', struct('gtSize', [100, 120, 32], 'starts', [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], ...
                'tileSize', [64, 72, 32], 'params', {{'xcorrShift', true}}), ...
            'gridOutBbox', struct('gtSize', [100, 120, 32], 'starts', [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], ...
                'tileSize', [64, 72, 32], 'params', {{'outBbox', [3, 5, 2, 90, 100, 30]}}), ...
            'tiles2D', struct('gtSize', [100, 120], 'starts', [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], ...
                'tileSize', [64, 72, 1], 'params', {{'blendMethod', 'feather'}}), ...
            'ioScan', struct('gtSize', [64, 104, 40], 'starts', [0, 0, 0; 0, 40, 0; 0, 0, 16; 0, 40, 16], ...
                'tileSize', [64, 64, 24], 'pixelSizes', [0.108, 0.108, 0.5], ...
                'params', {{'DSR', false, 'IOScan', true, 'dz', 0.5}}), ...
            'skewedSpace', struct('gtSize', [100, 64, 50], 'starts', [0, 0, 0; 36, 0, 0; 0, 0, 10; 36, 0, 10], ...
                'tileSize', [64, 64, 40], 'stageXYZ', [0, 0, 0; 0, 36 * 0.108, 0; 3, 0, 0; 3, 36 * 0.108, 0], ...
                'tileIdx', [0, 0, 0; 0, 1, 0; 1, 0, 0; 1, 1, 0], ...
                'params', {{'DSR', false, 'dz', 0.3, 'skewAngle', 32.45, 'reverse', true, 'axisOrder', '-x,y,z'}}), ...
            'resampledDSR', struct('gtSize', [80, 108, 32], 'starts', [0, 0, 0; 0, 36, 0; 30, 0, 0; 30, 36, 0], ...
                'tileSize', [50, 72, 32], 'pixelSizes', [0.216, 0.108, 0.108], ...
                'params', {{'resampleType', 'given', 'resampleFactor', [2, 1, 1]}}))
    end

    methods (Test, TestTags = {'Integration'})
        function exactPositionsReproduceGroundTruth(testCase, blendMethod)
            % 2 x 2 tiles in DSR space, as in demo_zarr_stitching.m
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);

            testCase.stitch(dataPath, csvFn, 'blendMethod', blendMethod);

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(size(out), [102, 122, 34]);
            testCase.verifyEqual(out(1 : 100, 1 : 120, 1 : 32), gt, 'The stitched volume differs from the ground truth.');
            testCase.verifyEqual(nnz(out) - nnz(out(1 : 100, 1 : 120, 1 : 32)), 0, 'Voxels outside the tiles are nonzero.');

            mip = readtiff(sprintf('%smatlab_stitch/MIPs/%s_MIP_z.tif', dataPath, testCase.StitchName));
            testCase.verifyEqual(mip, max(out, [], 3), 'The z MIP does not match the stitched volume.');

            info = load(sprintf('%smatlab_stitch/stitchInfo/%s.mat', dataPath, testCase.StitchName));
            testCase.verifyEqual(info.pImSz, [102, 122, 34]);
            testCase.verifyEqual(size(info.pTileSizes, 1), 4);
        end

        function featherBlendingOfUint8Tiles(testCase)
            % the feather blending mex takes uint16 and single tiles; uint8
            % tiles (kept as uint8 without save16bit) are blended by the
            % MATLAB code, and its NaNs replaced in replace_nan_with_value
            % (the MATLAB code where it isn't built)
            dataPath = testCase.createTempFolder();
            gt = uint8(testCase.randomVolume([100, 120, 32], [1, 255]));
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);

            testCase.stitch(dataPath, csvFn, 'blendMethod', 'feather', 'save16bit', false);

            out = testCase.readResult(dataPath);
            testCase.verifyClass(out, 'uint8');
            testCase.verifyEqual(out(1 : 100, 1 : 120, 1 : 32), gt, 'The stitched volume differs from the ground truth.');
            testCase.verifyEqual(nnz(out) - nnz(out(1 : 100, 1 : 120, 1 : 32)), 0, 'Voxels outside the tiles are nonzero.');
        end

        function featherBlendingOfUint8TilesInABatchThatIsNotFull(testCase)
            % two tiles with a narrow overlap: a batch that isn't full, whose
            % overlap is blended into the block of the tiles. uint8 tiles
            % kept the second tile's values there; they must blend as the
            % same tiles in uint16 (the mex, which may round a voxel the
            % other way where the compiler fuses its multiply-adds). Tile 2
            % is offset by 100, so blended voxels have values of neither tile.
            gt = testCase.randomVolume([64, 104, 32], [1, 120]);
            ims = {gt(:, 1 : 64, :), gt(:, 41 : 104, :) + 100};
            dtypes = {'uint8', 'uint16'};
            outs = cell(1, 2);
            for k = 1 : 2
                dataPath = testCase.createTempFolder();
                tiles = struct('im', cellfun(@(x) cast(x, dtypes{k}), ims, 'unif', 0), 'idx', {[0, 0, 0], [1, 0, 0]}, ...
                    'xyz', {[0, 0, 0], [40 * testCase.XYPixelSize, 0, 0]});
                [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);
                testCase.stitch(dataPath, csvFn, 'blendMethod', 'feather', 'save16bit', false);
                outs{k} = testCase.readResult(dataPath);
            end
            testCase.verifyClass(outs{1}, 'uint8');
            testCase.verifyEqual(double(outs{1}), double(outs{2}), 'AbsTol', 1, ...
                'uint8 tiles blend unlike the same tiles in uint16.');
            overlap = double(outs{1}(1 : 64, 41 : 64, 1 : 32));
            testCase.verifyTrue(any(overlap ~= double(ims{1}(:, 41 : 64, :)) & overlap ~= double(ims{2}(:, 1 : 24, :)), 'all'), ...
                'The overlap is not blended.');
        end

        function mipIsMadeBeforePublishing(testCase)
            % the stitching wrapper and the pipeline count a volume as
            % stitched once its result exists; the stitched MIP was written
            % after the result was moved into place
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            testCase.stitch(dataPath, csvFn);
            resultFn = sprintf('%smatlab_stitch/%s.zarr', dataPath, testCase.StitchName);
            testCase.assertTrue(moves.wasMoved(resultFn), 'The result was not moved into place.');
            published = moves.filesWhenMoved(resultFn);
            testCase.verifyTrue(ismember(sprintf('MIPs/%s_MIP_z.tif', testCase.StitchName), published), ...
                'The MIP was written after the result was moved into place.');
            added = setdiff(pktest.MovefileSpyFixture.listFiles([dataPath, 'matlab_stitch']), ...
                [published; {[testCase.StitchName, '.zarr']}]);
            testCase.verifyEmpty(added, 'Files were written after the result was moved into place.');
        end

        function mipIsMergedFromBlendingTasks(testCase)
            % the blending tasks save the z MIPs of their batches, so the
            % result isn't read again for its MIP. With mean blending of
            % tiles that differ by 1 the overlap ends in .5, which the uint16
            % result rounds as the MIP must. Small batches, so that a task
            % merges the MIPs of batches at the same y and x.
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 72, 32]);
            tiles = struct('im', {gt(1 : 64, :, :), gt(37 : 100, :, :) + 1}, 'idx', {[0, 0, 0], [0, 1, 0]}, ...
                'xyz', {[0, 0, 0], [0, 36 * testCase.XYPixelSize, 0]});
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);
            % memory for an in-memory MIP whatever the job's limit
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 1024));
            spy = testCase.applyFixture(pktest.FunctionSpyFixture('saveMIP_zarr'));
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            testCase.stitch(dataPath, csvFn, 'blendMethod', 'mean', 'batchSize', [64, 64, 16], 'blockSize', [16, 16, 16]);
            % a task's flag, which holds its MIPs, marks it done: it is
            % written under a temporary name first
            calls = moves.getCalls();
            testCase.verifyTrue(any(contains({calls.dst}, '/zarr_flags/') & endsWith({calls.dst}, '.mat')), ...
                'The blending flags were not moved into place.');
            out = testCase.readResult(dataPath);
            testCase.assertEqual(out(40, 10, 5), uint16(single(gt(40, 10, 5)) + 0.5), ...
                'The overlap should be the rounded mean of the tiles.');
            mip = readtiff(sprintf('%smatlab_stitch/MIPs/%s_MIP_z.tif', dataPath, testCase.StitchName));
            testCase.verifyEqual(mip, max(out, [], 3), 'The z MIP does not match the stitched volume.');
            testCase.verifyEmpty(spy.getCalls(), 'The MIP was computed from the stitched volume.');
        end

        function resultOfAnotherRunIsKept(testCase)
            % a run that lost the race to another run of the same volume (a
            % rerun while the earlier jobs still run) wrote the MIP of its
            % own unfinished result over the published one
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);
            testCase.stitch(dataPath, csvFn);
            resultFn = sprintf('%smatlab_stitch/%s.zarr', dataPath, testCase.StitchName);
            mipFn = sprintf('%smatlab_stitch/MIPs/%s_MIP_z.tif', dataPath, testCase.StitchName);
            otherPath = testCase.createTempFolder();
            movefile(resultFn, [otherPath, 'result.zarr']);
            movefile(mipFn, [otherPath, 'mip.tif']);
            % the other run publishes its result and MIP while this run blends
            testCase.applyFixture(pktest.FunctionSpyFixture('processStitchBlock', ...
                @() cellfun(@copyfile, {[otherPath, 'result.zarr'], [otherPath, 'mip.tif']}, {resultFn, mipFn})));
            testCase.stitch(dataPath, csvFn);
            testCase.verifyEqual(readtiff(mipFn), readtiff([otherPath, 'mip.tif']), 'The published MIP was overwritten.');
            testCase.verifyEqual(readzarr(resultFn), readzarr([otherPath, 'result.zarr']));
            testCase.verifyEmpty(dir(sprintf('%smatlab_stitch/%s.zarr_*', dataPath, testCase.StitchName)), ...
                'The temporary result was left behind.');
        end

        function largeVolumeMIPsAreMadeBeforePublishing(testCase)
            % volumes over half the memory get their MIPs from XR_MIP_zarr,
            % which names them after its input: the temporary result must
            % be named after the final one
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 0));
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            testCase.stitch(dataPath, csvFn);
            resultFn = sprintf('%smatlab_stitch/%s.zarr', dataPath, testCase.StitchName);
            testCase.assertTrue(moves.wasMoved(resultFn), 'The result was not moved into place.');
            mips = strcat('MIPs/', testCase.StitchName, {'_MIP_y.tif'; '_MIP_x.tif'; '_MIP_z.tif'});
            testCase.verifyEqual(ismember(mips, moves.filesWhenMoved(resultFn)), true(3, 1), ...
                'MIPs were missing when the result was moved into place.');
            out = testCase.readResult(dataPath);
            mip = readtiff(sprintf('%smatlab_stitch/MIPs/%s_MIP_z.tif', dataPath, testCase.StitchName));
            testCase.verifyEqual(mip, max(out, [], 3), 'The z MIP does not match the stitched volume.');
        end

        function separateBlendWritesEachTileOnTheCanvas(testCase)
            % blendMethod 'separate' writes each tile, padded to the stitched
            % canvas, as its own result (cropped with XR_crop_zarr, with MIPs)
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            starts = [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0];
            tileSize = [64, 72, 32];
            csvFn = testCase.writeTiles(dataPath, gt, starts, tileSize);

            testCase.stitch(dataPath, csvFn, 'blendMethod', 'separate');

            resultPath = [dataPath, 'matlab_stitch/'];
            testCase.verifyTrue(isfolder([resultPath, testCase.StitchName, '.zarr']), 'The stitched result is missing.');
            outs = dir([resultPath, '*.zarr']);
            outs = setdiff({outs.name}, {[testCase.StitchName, '.zarr']});
            testCase.assertNumElements(outs, 4);
            placed = zeros(0, 3);
            for i = 1 : numel(outs)
                out = readzarr([resultPath, outs{i}]);
                testCase.verifyEqual(size(out), [102, 122, 34], outs{i});
                % the tile sits at its position, with zeros around it
                [y, x, z] = ind2sub(size(out), find(out));
                s = [min(y), min(x), min(z)] - 1;
                placed(end + 1, :) = s;
                tile = zeros(size(out), 'like', out);
                tile(s(1) + 1 : s(1) + tileSize(1), s(2) + 1 : s(2) + tileSize(2), s(3) + 1 : s(3) + tileSize(3)) = ...
                    gt(s(1) + 1 : s(1) + tileSize(1), s(2) + 1 : s(2) + tileSize(2), s(3) + 1 : s(3) + tileSize(3));
                testCase.verifyEqual(out, tile, outs{i});
                mip = readtiff(sprintf('%sMIPs/%s_MIP_z.tif', resultPath, outs{i}(1 : end - 5)));
                testCase.verifyEqual(mip, max(out, [], 3), outs{i});
            end
            testCase.verifyEqual(sortrows(placed), sortrows(starts), 'The tiles are not at their positions.');
        end

        function processFunctionFileHasANameOfItsOwn(testCase)
            % with EdgeArtifacts, a stitching run saves its tile processing
            % function to tmp/. The file was named by the minute and the
            % hundredths of a second ('yyyyMMddHHmmSS'), so jobs of other
            % volumes starting in the same minute could write the same file;
            % it now gets a uuid of its own.
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);
            testCase.applyFixture(pktest.FunctionStubFixture('get_uuid', 'a1b2c3d4'));
            testCase.stitch(dataPath, csvFn, 'EdgeArtifacts', 2);
            fns = dir(fullfile(dataPath, '**', 'processFunction_*.mat'));
            testCase.verifyEqual({fns.name}, {'processFunction_a1b2c3d4.mat'});
        end

        function masterDoesNotStitchAVolumeWhoseJobRuns(testCase)
            % with masterCompute, the master stitches the last unfinished
            % volume itself. Once that one was done, the last unfinished
            % volume could be one whose job already ran, which the master
            % then stitched a second time with the same temporary files (a
            % corrupt stitch info file in a real run); it now waits for the
            % job. The stitching function is a spy that makes the result of
            % the second volume, so each call is one the master made.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [dataPath, csvFn, resultFns] = testCase.writeTwoVolumes();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            % the job of the first volume runs and makes its result when first checked
            slurm.setJobStates('1000_1', {'R', 'done'}, resultFns(1));
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            frames = testCase.applyFixture(pktest.FunctionSpyFixture('XR_stitching_frame_zarr_dev_v1', ...
                @() fclose(fopen(resultFns{2}, 'w'))));
            testCase.stitch(dataPath, csvFn, 'parseCluster', true, 'masterCompute', true, 'xcorrMode', 'all', ...
                'jobLogDir', [dataPath, 'logs']);
            testCase.verifyNumElements(frames.getCalls(), 1, 'The master stitched a volume whose job was running.');
            testCase.verifySubstring(frames.getOption(1, 'resultPath'), 'Scan_Iter_0001', ...
                'The master did not stitch the second volume.');
        end

        function masterInAJobStitchesNoVolumeThatDoesNotFitIt(testCase)
            % with masterCompute the master stitches the last unfinished
            % volume itself, also in a Slurm job without the memory for it
            % (the pipeline's stitching job, sized for one tile); now it
            % submits it then. The stitching function is a spy, so each
            % call is one the master made.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [dataPath, csvFn, resultFns] = testCase.writeTwoVolumes();
            orig = getenv('SLURM_JOB_ID');
            if isempty(orig)
                testCase.addTeardown(@unsetenv, 'SLURM_JOB_ID');
            else
                testCase.addTeardown(@setenv, 'SLURM_JOB_ID', orig);
            end
            setenv('SLURM_JOB_ID', '4242');
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 1e-6));
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('1000_1', {'R', 'done'}, resultFns(1));
            slurm.setJobStates('1001_2', {'R', 'done'}, resultFns(2));
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            frames = testCase.applyFixture(pktest.FunctionSpyFixture('XR_stitching_frame_zarr_dev_v1'));
            testCase.stitch(dataPath, csvFn, 'parseCluster', true, 'masterCompute', true, 'xcorrMode', 'all', ...
                'jobLogDir', [dataPath, 'logs']);
            testCase.verifyEmpty(frames.getCalls(), 'The master stitched a volume that does not fit its job.');
            testCase.verifyEqual(slurm.callCount('sbatch'), 2);
        end

        function clusterJobsAskForWhatTheyHold(testCase)
            % the stitching job of a volume asked for the memory of all its
            % tiles (8 B per voxel) on at least cpusPerTask (8) CPUs, which
            % it held while it waited for its block jobs (80 such jobs held
            % 840 CPUs, and none of their 624 block jobs could start). Now
            % it asks for what it holds, on one CPU, or with masterCompute
            % (here the config's) for a block's job (2 CPUs), unless it
            % stitches the volume itself (fewer than 30 batches: as
            % before); the wrapper's masterCompute reaches it, and its jobs
            % go first (nice). The stitching function is a spy; the fake
            % jobs make the results.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.applyFixture(pktest.FunctionStubFixture('coordinator_nice', 1000));
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 60));
            testCase.applyFixture(pktest.FunctionSpyFixture('XR_stitching_frame_zarr_dev_v1'));
            configFolder = testCase.createTempFolder();
            onConfig = [configFolder, 'cpu_config.json'];
            fid = fopen(onConfig, 'w');
            fprintf(fid, '{"masterCompute": true}');
            fclose(fid);
            % batch size, config, CPUs of the jobs, more parameters (168
            % batches; 1 batch; with registration, whose max offsets could
            % make it 1 batch: stitched in its job as before, or small
            % offsets: still 168)
            cases = {[16, 16, 16], '', 1, {}; [16, 16, 16], onConfig, 2, {}; [512, 512, 512], '', 8, {}; ...
                [16, 16, 16], '', 8, {'xcorrShift', true}; ...
                [16, 16, 16], '', 1, {'xcorrShift', true, 'xyMaxOffset', 2, 'zMaxOffset', 1}};
            for i = 1 : size(cases, 1)
                [dataPath, csvFn, resultFns] = testCase.writeTwoVolumes();
                slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
                slurm.setJobStates('1000_1', {'R', 'done'}, resultFns(1));
                slurm.setJobStates('1001_2', {'R', 'done'}, resultFns(2));
                testCase.stitch(dataPath, csvFn, 'parseCluster', true, 'masterCompute', false, 'xcorrMode', 'all', ...
                    'batchSize', cases{i, 1}, 'blockSize', [16, 16, 16], 'configFile', cases{i, 2}, ...
                    'jobLogDir', [dataPath, 'logs'], cases{i, 4}{:});
                lines = slurm.calls('sbatch');
                testCase.assertNumElements(lines, 2, sprintf('Case %d.', i));
                cpus = cellfun(@(x) str2double(regexp(x, '--cpus-per-task=(\d+)', 'tokens', 'once')), lines);
                testCase.verifyEqual(cpus(:)', [1, 1] * cases{i, 3}, sprintf('Case %d: the jobs'' CPUs.', i));
                testCase.verifyTrue(all(contains(lines, '--nice=1000')), sprintf('Case %d: no nice.', i));
                testCase.verifyTrue(all(contains(lines, '''masterCompute'',false')), ...
                    sprintf('Case %d: the wrapper''s masterCompute did not reach the jobs.', i));
            end
        end

        function volumeSizeEstimateMatchesTheStitchedVolume(testCase, sizeSetup)
            % the stitching wrapper sizes a volume's job by the size of the
            % volume (stitch_volume_size_estimate, which mirrors the
            % stitching function's geometry): it is the size of the
            % stitched volume (exact tiles: no registration shifts)
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume(sizeSetup.gtSize);
            pixelSizes = testCase.XYPixelSize * [1, 1, 1];
            if isfield(sizeSetup, 'pixelSizes')
                pixelSizes = sizeSetup.pixelSizes;
            end
            nT = size(sizeSetup.starts, 1);
            stageXYZ = sizeSetup.starts(:, [2, 1, 3]) .* pixelSizes([2, 1, 3]);
            if isfield(sizeSetup, 'stageXYZ')
                stageXYZ = sizeSetup.stageXYZ;
            end
            tileIdx = [];
            if isfield(sizeSetup, 'tileIdx')
                tileIdx = sizeSetup.tileIdx;
            end
            csvFn = testCase.writeTiles(dataPath, gt, sizeSetup.starts, sizeSetup.tileSize, pixelSizes, stageXYZ, tileIdx);
            testCase.stitch(dataPath, csvFn, sizeSetup.params{:});
            out = testCase.readResult(dataPath);
            actual = [size(out), ones(1, 3 - ndims(out))];

            % the stitching defaults of the test, then the setup's
            params = [{'DSR', true, 'resampleType', 'isotropic', 'xyPixelSize', testCase.XYPixelSize, ...
                'axisOrder', 'x,y,z', 'xcorrShift', false}, sizeSetup.params];
            known = {'DSR', 'resampleType', 'resampleFactor', 'xyPixelSize', 'axisOrder', 'dataOrder', 'xcorrShift', ...
                'dz', 'skewAngle', 'IOScan', 'objectiveScan', 'outBbox', 'batchSize', 'blockSize'};
            keep = ismember(lower(params(1 : 2 : end)), lower(known));
            params = reshape(params, 2, []);
            params = params(:, keep);
            if isempty(tileIdx)
                xyzStarts = sizeSetup.starts(:, [2, 1, 3]);
                tileIdx = zeros(nT, 3);
                for d = 1 : 3
                    [~, ~, tileIdx(:, d)] = unique(xyzStarts(:, d));
                end
                tileIdx = tileIdx - 1;
            end
            [estimate, ~, range] = stitch_volume_size_estimate(stageXYZ, sizeSetup.tileSize, ...
                'tileIdx', [tileIdx, zeros(nT, 1)], params{:});
            testCase.verifyEqual(estimate, actual);
            % (with registration, smaller or larger by up to its max
            % offsets, 300 in x and y and 50 in z by default, at both ends)
            k = find(strcmpi(params(1, :), 'xcorrShift'), 1, 'last');
            margin = 2 * [300, 300, 50] * (params{2, k} && nT > 1 && ~any(strcmpi(params(1, :), 'outBbox')));
            margin(estimate == 1) = 0;
            testCase.verifyEqual(range, [max(estimate - margin, 1); estimate + margin]);
        end

        function masterWaitsForTheLastTryOfAJob(testCase)
            % with masterCompute, once the last try of every unfinished
            % volume was submitted, the wrapper errored while waiting for
            % those jobs (no last volume for the master to take over)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [dataPath, csvFn, resultFns] = testCase.writeTwoVolumes();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('1000_1', {'PD', 'R', 'done'});
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            testCase.applyFixture(pktest.FunctionSpyFixture('XR_stitching_frame_zarr_dev_v1', ...
                @() fclose(fopen(resultFns{2}, 'w'))));
            testCase.stitch(dataPath, csvFn, 'parseCluster', true, 'masterCompute', true, 'xcorrMode', 'all', ...
                'maxTrialNum', 1, 'jobLogDir', [dataPath, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 1, 'The last try was submitted again.');
            checks = slurm.history();
            checks = checks(startsWith(checks, 'squeue') & contains(checks, '1000_1'));
            testCase.verifyNumElements(checks, 3, 'The job was not checked until it ended.');
        end

        function lastTryWaitsForATaskRunnerThatHasTheTask(testCase)
            % (task runners) the second volume's job makes its result and
            % then runs the third volume's task (a task runner; the first
            % volume, the first of its kind, runs as before), so the third
            % volume's own job ends at once: the wrapper counts it as
            % running while the second volume's job runs, instead of as
            % finished, and does not submit it again
            testCase.assumeTrue(isunix && ~ismac, 'Task runners run on Linux only.');
            [dataPath, csvFn, resultFns] = testCase.writeTwoVolumes(3);
            logPath = [dataPath, 'logs/'];
            mkdir(logPath);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('1000_1', {'R', 'done'}, resultFns(1));
            slurm.setJobStates('1001_2', {'R', 'R', 'R', 'R', 'done'}, resultFns(2));
            % after the first round, the second volume's job (1001, array
            % task 2) has the third volume's task (the second of the pool)
            claim = sprintf(['for d in "%s"task_runner_single_*/claims; do [ -d "$d" ] && ', ...
                'ln -s "1001_2:othernode:5:1" "$d/t2" 2> /dev/null; done; true'], logPath);
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() system(claim), 40));
            testCase.stitch(dataPath, csvFn, 'parseCluster', true, 'masterCompute', false, 'xcorrMode', 'all', ...
                'maxTrialNum', 1, 'jobLogDir', logPath);
            testCase.verifyEqual(slurm.callCount('sbatch'), 3, 'The third volume was submitted again.');
            checks = slurm.history();
            checks = checks(startsWith(checks, 'squeue') & contains(checks, '1001_2'));
            testCase.verifyGreaterThanOrEqual(numel(checks), 5, 'The wrapper did not wait for the job that had the third volume.');
            % (its jobs gone, the master removed the pool at the end)
            testCase.verifyEmpty(dir([logPath, 'task_runner_single_*']), 'The pool is still there.');
        end

        function clusterModeSeesResultsOtherNodesJustMade(testCase)
            % on NFS, exist() can miss a result another node just made for 30
            % s or more, so the wrapper saw finished volumes (and secondary
            % channels the primary's stitch info) a round late; it now checks
            % them with batch_file_exist, which sees them at once. The stub
            % reports the result made, so no job is submitted.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.applyFixture(pktest.FunctionStubFixture('batch_file_exist', true));
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            testCase.stitch(dataPath, csvFn, 'parseCluster', true, 'masterCompute', false, 'maxTrialNum', 1, ...
                'jobLogDir', [dataPath, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 0, 'A stitching job was submitted for a finished result.');
        end

        function processedTilesOtherNodesJustMadeAreSeen(testCase)
            % on NFS, exist() can miss tiles another node just made for 30 s
            % or more, so the wrapper ran out of tries for a time point whose
            % deskewed tiles were done and left it to a second stitching job
            % (the pipeline demo then took 110-140 s instead of 65 s); it now
            % checks them with batch_file_exist. The stub reports the tiles
            % missing although they are on disk, so the wrapper must wait for
            % them instead of stitching.
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32], ...
                [], [], [], 'tileDir', 'DSR', 'saveFormat', 'zarr');
            testCase.applyFixture(pktest.FunctionStubFixture('batch_file_exist', false));
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            frames = testCase.applyFixture(pktest.FunctionSpyFixture('XR_stitching_frame_zarr_dev_v1'));
            testCase.stitch(dataPath, csvFn, 'zarrFile', true, 'processedDirStr', 'DSR', 'maxTrialNum', 1);
            testCase.verifyEmpty(frames.getCalls(), 'The wrapper stitched tiles that batch_file_exist did not see.');
        end

        function clusterModeWaitsForTheLastTry(testCase)
            % after submitting the last try, the wrapper returned without
            % waiting for its job; it now checks the job until it ends,
            % without submitting again
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('1000_1', {'PD', 'R', 'done'});
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            testCase.stitch(dataPath, csvFn, 'parseCluster', true, 'masterCompute', false, 'maxTrialNum', 1, ...
                'jobLogDir', [dataPath, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 1, 'The last try was submitted again.');
            % one check per round until the job is gone (none before
            % submitting: the task had no job)
            testCase.verifyEqual(slurm.callCount('squeue'), 3);
        end

        function clusterModeEndsSoonAfterTheLastResult(testCase)
            % the wrapper looked for its jobs' results in rounds 5 s apart,
            % so it ended up to 5 s after the last one; it now looks for the
            % results still to come each second and goes on once they are
            % all there. The job runs; its result comes during the second
            % second.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);
            resultFn = sprintf('%smatlab_stitch/%s.zarr', dataPath, testCase.StitchName);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('1000_1', {'R'});
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [dataPath, 'calls.txt'], 2, @() mkdir(resultFn)), 20));
            testCase.stitch(dataPath, csvFn, 'parseCluster', true, 'masterCompute', false, 'maxTrialNum', 1, ...
                'jobLogDir', [dataPath, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 1);
            pktest.verifyWaitsOfASecond(testCase, pauses, 2, ...
                'The wrapper did not end within a second of the result.');
        end

        function clusterModeWaitsForTilesJobsStillMake(testCase)
            % in cluster mode the wrapper gave up a volume whose processed
            % tiles were missing after maxTrialNum checks 5 s apart, also
            % while the pipeline's deskew jobs were still making them (the
            % pipeline then needed another stitching job for it, or lost
            % it). While the pipeline's marker in their folder is recent,
            % it now waits for them without counting tries. The tiles come
            % at the third wait; the master stitches the volume (the spy
            % makes its result).
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [dataPath, csvFn, stagePath] = testCase.writeTilesToCome();
            fclose(fopen([dataPath, 'DSR/deskew_in_progress.tmp'], 'w'));
            resultFn = sprintf('%smatlab_stitch/%s.zarr', dataPath, testCase.StitchName);
            testCase.applyFixture(pktest.FakeSlurmFixture);
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', ...
                @() testCase.atCall([stagePath, 'calls.txt'], 3, @() movefile([stagePath, 'DSR/*'], [dataPath, 'DSR/'])), 20));
            frames = testCase.applyFixture(pktest.FunctionSpyFixture('XR_stitching_frame_zarr_dev_v1', ...
                @() fclose(fopen(resultFn, 'w'))));
            testCase.stitch(dataPath, csvFn, 'processedDirStr', 'DSR', 'parseCluster', true, 'masterCompute', true, ...
                'maxTrialNum', 1, 'jobLogDir', [dataPath, 'logs']);
            testCase.verifyNumElements(frames.getCalls(), 1, 'The wrapper gave up tiles that jobs still made.');
            % it looks for them each second (and not after 5 s rounds), so
            % that it stitches soon after they come
            pktest.verifyWaitsOfASecond(testCase, pauses, 3, 'The wrapper did not look for the tiles each second.');
        end

        function clusterModeLooksForThePrimaryStitchInfoEachSecond(testCase)
            % a secondary channel waits for the stitch info of the primary
            % channel, which the wrapper looked for in rounds 5 s apart; in
            % cluster mode it now looks for it each second, so that the
            % secondary channel is submitted soon after it comes, also
            % while another volume waits for tiles that jobs still make.
            % Time point 0: the primary channel is done, its stitch info
            % comes at the second wait, and the job of the secondary
            % channel makes its result; time point 1 (primary channel)
            % waits for its tiles until the marker goes (third wait).
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([64, 104, 32]);
            stageXYZ = [0, 0, 0; 40 * testCase.XYPixelSize, 0, 0];
            tileIdx = [0, 0, 0; 1, 0, 0];
            tiles = {[0, 0, 0; 0, 40, 0], [64, 64, 32], [], stageXYZ, tileIdx};
            testCase.writeTiles(dataPath, gt, tiles{:}, 'cam', 'A');
            testCase.writeTiles(dataPath, gt, tiles{:}, 'cam', 'B');
            csvFn = testCase.writeTiles(dataPath, gt, tiles{:}, 'cam', 'A', 'iter', 1);
            for cam = 'AB'
                testCase.writeTiles(dataPath, gt, tiles{:}, 'cam', cam, 'tileDir', 'DSR', 'imageList', 'ImageList_DSR.csv');
            end
            markerFn = [dataPath, 'DSR/deskew_in_progress.tmp'];
            fclose(fopen(markerFn, 'w'));
            stitchPath = [dataPath, 'matlab_stitch/'];
            mkdir([stitchPath, 'stitchInfo']);
            fclose(fopen([stitchPath, testCase.StitchName, '.zarr'], 'w'));
            infoFn = [stitchPath, 'stitchInfo/', testCase.StitchName, '.mat'];
            resultFn = [stitchPath, strrep(testCase.StitchName, 'CamA', 'CamB'), '.zarr'];
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            % the secondary channel of time point 0 is the third volume
            slurm.setJobStates('1000_3', {'R', 'done'}, {resultFn});
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() testCase.atCall( ...
                [dataPath, 'calls.txt'], [2, 3], {@() fclose(fopen(infoFn, 'w')), @() delete(markerFn)}), 20));
            testCase.stitch(dataPath, csvFn, 'channelPatterns', {'CamA_ch0', 'CamB_ch0'}, 'processedDirStr', 'DSR', ...
                'parseCluster', true, 'masterCompute', false, 'maxTrialNum', 1, 'jobLogDir', [dataPath, 'logs']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 1);
            % waits of a second: two for the stitch info, five (a round) for
            % the tiles, which no longer come
            pktest.verifyWaitsOfASecond(testCase, pauses, 7, ...
                'The wrapper did not look for the stitch info and tiles each second.');
        end

        function clusterModeGivesUpTilesWithoutARecentMarker(testCase, markerAge)
            % without a pipeline that keeps its marker recent (none, or a
            % marker left from one that ended), missing processed tiles
            % count as tries, as before
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [dataPath, csvFn] = testCase.writeTilesToCome();
            if markerAge >= 0
                markerFn = [dataPath, 'DSR/deskew_in_progress.tmp'];
                fclose(fopen(markerFn, 'w'));
                testCase.assumeTrue(usejava('jvm'), 'Setting file times needs Java.');
                java.io.File(markerFn).setLastModified(java.lang.System.currentTimeMillis() - markerAge * 60 * 1000);
            end
            testCase.applyFixture(pktest.FakeSlurmFixture);
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            frames = testCase.applyFixture(pktest.FunctionSpyFixture('XR_stitching_frame_zarr_dev_v1'));
            testCase.stitch(dataPath, csvFn, 'processedDirStr', 'DSR', 'parseCluster', true, 'masterCompute', true, ...
                'maxTrialNum', 2, 'jobLogDir', [dataPath, 'logs']);
            testCase.verifyEmpty(frames.getCalls(), 'The wrapper stitched missing tiles.');
            testCase.verifyNumElements(pauses.getCalls(), 1, 'The tiles did not count as tries.');
        end

        function streamingSeesTilesOtherMachinesJustMade(testCase)
            % in streaming mode the tiles were checked with exist(), which
            % on NFS can miss files another machine just made for 30 s or
            % more; they are now checked with batch_file_exist. The stub
            % reports all files made (the result too, so nothing is
            % stitched), while there are no raw tiles on disk: the wrapper
            % must not wait for them.
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32], ...
                [], [], [], 'tileDir', 'DSR');
            testCase.applyFixture(pktest.FunctionStubFixture('batch_file_exist', true));
            % (each wait notes the batch_file_exist in effect, for the
            % diagnostic)
            noteFn = [dataPath, 'pauses.txt'];
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', ...
                @() writelines(string(which('batch_file_exist')), noteFn, 'WriteMode', 'append'), 20));
            testCase.stitch(dataPath, csvFn, 'processedDirStr', 'DSR', 'streaming', true);
            notes = '';
            if isfile(noteFn)
                notes = strjoin(readlines(noteFn), '; ');
            end
            testCase.verifyEmpty(pauses.getCalls(), ['The wrapper waited for tiles batch_file_exist reported made ', ...
                '(batch_file_exist at the waits: ', notes, ').']);
        end

        function secondaryChannelUsesTheTaskSizeOfThePrimary(testCase)
            % the secondary channel reads the per-task block info of the
            % primary channel, whose task size depends on the local
            % workers it was stitched with. The stub leaves no memory for
            % the workers, so all tasks run in this MATLAB.
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 0));
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            gtA = testCase.randomVolume([64, 104, 32]);
            gtB = testCase.randomVolume([64, 104, 32]);
            stageXYZ = [0, 0, 0; 40 * px, 0, 0];
            tileIdx = [0, 0, 0; 1, 0, 0];
            testCase.writeTiles(dataPath, gtA, [0, 0, 0; 0, 40, 0], [64, 64, 32], [], stageXYZ, tileIdx, 'cam', 'A');
            csvFn = testCase.writeTiles(dataPath, gtB, [0, 0, 0; 0, 40, 0], [64, 64, 32], [], stageXYZ, tileIdx, 'cam', 'B');
            opts = {'primaryCh', 'CamA_ch0', 'blendMethod', 'feather', 'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16]};

            testCase.stitch(dataPath, csvFn, opts{:}, 'channelPatterns', {'CamA_ch0'}, 'localWorkerNum', 4);
            testCase.stitch(dataPath, csvFn, opts{:}, 'channelPatterns', {'CamA_ch0', 'CamB_ch0'}, 'localWorkerNum', 1);

            outB = testCase.readResult(dataPath, 'matlab_stitch', strrep(testCase.StitchName, 'CamA', 'CamB'));
            testCase.verifyEqual(outB(1 : 64, 1 : 104, 1 : 32), gtB);
        end

        function secondaryWithAShortTileStitchesWhatItHas(testCase)
            % a secondary (here time point 1) whose tile 4 is shorter in z
            % than the primary's (its acquisition stopped early) used the
            % primary's per-block info, made for the primary's tile sizes,
            % so its blocks read past the end of that tile (readzarr: "Upper
            % bound is invalid"). Now it makes its own, with the primary's
            % tile positions (its own tile 2 is one voxel off them) and
            % shifts: its result is that of its tiles stitched on their own
            % at those positions
            px = testCase.XYPixelSize;
            starts = [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0];
            dataPath = testCase.createTempFolder();
            alonePath = testCase.createTempFolder();
            testCase.writeTiles(dataPath, testCase.randomVolume([100, 120, 32]), starts, [64, 72, 32]);
            tiles = testCase.cutTiles(testCase.randomVolume([100, 120, 32]), starts, [64, 72, 32]);
            tiles(4).im = tiles(4).im(:, :, 1 : 20);
            pktest.writeStitchTiles(alonePath, tiles, 'iter', 1);
            tiles(2).xyz = tiles(2).xyz + [px, 0, 0];
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles, 'iter', 1);

            testCase.verifyWarning(@() testCase.stitch(dataPath, csvFn), 'XR_stitching_frame_zarr_dev_v1:tileSizes');
            testCase.stitch(alonePath, [alonePath, 'ImageList.csv']);

            name = strrep(testCase.StitchName, 'Iter_0000', 'Iter_0001');
            testCase.verifyEqual(testCase.readResult(dataPath, 'matlab_stitch', name), ...
                testCase.readResult(alonePath, 'matlab_stitch', name), ...
                'The time point with a short tile is not its tiles stitched on their own.');
            infoPath = [dataPath, 'matlab_stitch/stitchInfo/'];
            testCase.verifyNotEmpty(dir([infoPath, name, '_task_size_*']), 'The secondary has no block info of its own.');
            testCase.verifyNotEmpty(dir([infoPath, testCase.StitchName, '_task_size_*']));
        end

        function secondaryOfThePrimarySizesUsesItsBlockInfo(testCase)
            % a secondary whose tiles and stitched size are the primary's
            % still stitches with the primary's per-block info (none of its
            % own), as before: the same data, the same result
            [dataPath, csvFn] = testCase.writeTwoVolumes();
            testCase.stitch(dataPath, csvFn);

            name = strrep(testCase.StitchName, 'Iter_0000', 'Iter_0001');
            infoPath = [dataPath, 'matlab_stitch/stitchInfo/'];
            testCase.verifyNotEmpty(dir([infoPath, testCase.StitchName, '_task_size_*']));
            testCase.verifyEmpty(dir([infoPath, name, '_task_size_*']), 'The secondary made block info of its own.');
            testCase.verifyEqual(testCase.readResult(dataPath, 'matlab_stitch', name), testCase.readResult(dataPath));
        end

        function secondaryOfAnotherStitchedSizeGetsItsOwnBlockInfo(testCase)
            % a secondary stitched later with another outBbox than the
            % primary's (so of another stitched size) used the primary's
            % per-block info, made for the primary's volume. Now it makes its
            % own: its result is its box of the volume
            gt = testCase.randomVolume([100, 120, 32]);
            starts = [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0];
            dataPath = testCase.createTempFolder();
            csvFn = testCase.writeTiles(dataPath, gt, starts, [64, 72, 32]);
            testCase.stitch(dataPath, csvFn, 'outBbox', [1, 1, 1, 100, 120, 32]);
            testCase.writeTiles(dataPath, gt, starts, [64, 72, 32], [], [], [], 'iter', 1);

            testCase.verifyWarning(@() testCase.stitch(dataPath, csvFn, 'outBbox', [11, 21, 1, 90, 100, 32]), ...
                'XR_stitching_frame_zarr_dev_v1:stitchedSize');

            name = strrep(testCase.StitchName, 'Iter_0000', 'Iter_0001');
            testCase.verifyEqual(testCase.readResult(dataPath, 'matlab_stitch', name), gt(11 : 90, 21 : 100, :));
            testCase.verifyNotEmpty(dir([dataPath, 'matlab_stitch/stitchInfo/', name, '_task_size_*']), ...
                'The secondary has no block info of its own.');
        end

        function largeDatasetsUseLocalWorkers(testCase)
            % without localWorkerNum, the volumes of datasets with 32 or more
            % tiles per volume to register are stitched with local workers,
            % smaller ones one task at a time in this MATLAB
            testCase.assumeTrue(license('test', 'Distrib_Computing_Toolbox') && ~isempty(which('gcp')), ...
                'Needs the Parallel Computing Toolbox.');
            spy = testCase.applyFixture(pktest.FunctionSpyFixture('XR_stitching_frame_zarr_dev_v1'));
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            for nT = [4, 32]
                dataPath = testCase.createTempFolder();
                [x, y] = ndgrid(0 : 3, 0 : nT / 4 - 1);
                tiles = struct('im', repmat({ones(8, 8, 2, 'uint16')}, nT, 1), 'idx', num2cell([x(:), y(:), zeros(nT, 1)], 2), ...
                    'xyz', num2cell([x(:), y(:), zeros(nT, 1)] * 6 * testCase.XYPixelSize, 2));
                [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);
                testCase.stitch(dataPath, csvFn, 'xcorrShift', true, 'maxTrialNum', 1);
            end
            testCase.assertNumElements(spy.getCalls(), 2);
            testCase.verifyEmpty(spy.getOption(1, 'localWorkerNum'), 'The small dataset used local workers.');
            workerNum = max(1, min(8, floor(maxNumCompThreads / 3)));
            if workerNum > 1
                testCase.verifyEqual(spy.getOption(2, 'localWorkerNum'), workerNum);
            else
                testCase.verifyEmpty(spy.getOption(2, 'localWorkerNum'));
            end
        end

        function tiles2D(testCase)
            % single page tiles, as for the 2D phase data in
            % demo_phase_and_2photon_stitching.m
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 1]);

            testCase.stitch(dataPath, csvFn, 'blendMethod', 'feather');

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(size(out), [102, 122]);
            testCase.verifyEqual(out(1 : 100, 1 : 120), gt);
        end

        function ioScanUsesStagePositionsAsIs(testCase)
            % objective/IO scan tiles in x and z (z offsets in units of dz),
            % as for the 3D 2-photon data in demo_phase_and_2photon_stitching.m
            dataPath = testCase.createTempFolder();
            dz = 0.5;
            gt = testCase.randomVolume([64, 104, 40]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 0, 40, 0; 0, 0, 16; 0, 40, 16], [64, 64, 24], ...
                [testCase.XYPixelSize, testCase.XYPixelSize, dz]);

            testCase.stitch(dataPath, csvFn, 'DSR', false, 'IOScan', true, 'dz', dz, 'blendMethod', 'none');

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(size(out), [66, 106, 42]);
            testCase.verifyEqual(out(1 : 64, 1 : 104, 1 : 40), gt);
        end

        function skewedSpacePlacement(testCase, axisOrder)
            % stitching in skewed space (demo_skewed_space_stitching.m): a
            % stage x step of n * dz is n raw frames in z; stage y is y. With
            % '-x,y,z' the stage x direction is flipped.
            dataPath = testCase.createTempFolder();
            dz = 0.3;
            n = 10;
            gt = testCase.randomVolume([100, 64, 50]);
            % [y0, x0, z0] in the ground truth; stage X of the tiles is set below
            xFlipped = startsWith(axisOrder, '-');
            zStart = [n, 0];
            if xFlipped
                zStart = [0, n];
            end
            starts = [0, 0, zStart(1); 36, 0, zStart(1); 0, 0, zStart(2); 36, 0, zStart(2)];
            stageXYZ = [0, 0, 0; 0, 36 * testCase.XYPixelSize, 0; n * dz, 0, 0; n * dz, 36 * testCase.XYPixelSize, 0];
            tileIdx = [0, 0, 0; 0, 1, 0; 1, 0, 0; 1, 1, 0];
            csvFn = testCase.writeTiles(dataPath, gt, starts, [64, 64, 40], [], stageXYZ, tileIdx);

            testCase.stitch(dataPath, csvFn, 'DSR', false, 'dz', dz, 'skewAngle', 32.45, 'reverse', true, ...
                'axisOrder', axisOrder, 'blendMethod', 'feather');

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(size(out), [102, 66, 52]);
            testCase.verifyEqual(out(1 : 100, 1 : 64, 1 : 50), gt);
        end

        function skewedSpaceRegistration(testCase)
            % registration in skewed space with the demo settings ('-x,y,z',
            % feather): the y = 1 tiles are listed 4 voxels off in y
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            dz = 0.3;
            n = 10;
            gt = testCase.randomVolume([100, 64, 50]);
            starts = [0, 0, 0; 36, 0, 0; 0, 0, n; 36, 0, n];
            stageXYZ = [0, 0, 0; 0, 40 * px, 0; n * dz, 0, 0; n * dz, 40 * px, 0];
            tileIdx = [0, 0, 0; 0, 1, 0; 1, 0, 0; 1, 1, 0];
            csvFn = testCase.writeTiles(dataPath, gt, starts, [64, 64, 40], [], stageXYZ, tileIdx);

            testCase.stitch(dataPath, csvFn, 'DSR', false, 'dz', dz, 'skewAngle', 32.45, 'reverse', true, ...
                'axisOrder', '-x,y,z', 'blendMethod', 'feather', 'xcorrShift', true, 'xcorrDownsample', [1, 1, 1], ...
                'xyMaxOffset', 20, 'zMaxOffset', 10);

            testCase.verifyEqual(testCase.cropToData(testCase.readResult(dataPath)), gt, ...
                'Registration in skewed space did not recover the tile positions.');
        end

        function registration2DTiles(testCase)
            % 2D registration with xcorrMode 'primary' (all time points of the
            % primary channel), as for the phase data in
            % demo_phase_and_2photon_stitching.m
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            gt = testCase.randomVolume([64, 104]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 0, 40, 0], [64, 64, 1], [], ...
                [0, 0, 0; 44 * px, 0, 0], [0, 0, 0; 1, 0, 0]);

            testCase.stitch(dataPath, csvFn, 'DSR', false, 'IOScan', true, 'xcorrShift', true, 'xcorrMode', 'primary', ...
                'xcorrDownsample', [1, 1, 1], 'xyMaxOffset', 20, 'zMaxOffset', 10, 'blendMethod', 'feather', ...
                'blockSize', [64, 64, 1], 'batchSize', [256, 256, 1]);

            testCase.verifyEqual(testCase.cropToData(testCase.readResult(dataPath)), gt, ...
                '2D registration did not recover the tile position.');
        end

        function processedZarrTiles(testCase)
            % DSR tiles saved as zarr in a subfolder with the image list
            % refering to the raw tif names (ProcessedDirStr='DSR',
            % zarrFile=true), as in demo_zarr_stitching.m
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32], ...
                [], [], [], 'tileDir', 'DSR', 'saveFormat', 'zarr');

            testCase.stitch(dataPath, csvFn, 'zarrFile', true, 'processedDirStr', 'DSR', 'blendMethod', 'feather');

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(size(out), [102, 122, 34]);
            testCase.verifyEqual(out(1 : 100, 1 : 120, 1 : 32), gt);
        end

        function processFunctionIsAppliedToTiles(testCase)
            % user process function saved in a .mat file, as the flat field
            % correction in demo_phase_and_2photon_stitching.m. With a flat
            % flat-field image, background 100 and offset 1, x -> x - 99.
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);

            ffFn = [dataPath, 'flatfield/flatfield.tif'];
            bgFn = [dataPath, 'flatfield/background.tif'];
            writetiff(1000 * ones(64, 72, 'uint16'), ffFn);
            writetiff(100 * ones(64, 72, 'uint16'), bgFn);
            usrFun = sprintf("@(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'%s','%s',%d,%d,%s),%0.10d)", ...
                ffFn, bgFn, 1, 0.4, string(false), 0);
            funFn = [dataPath, 'processFunction.mat'];
            save('-v7.3', funFn, 'usrFun');

            testCase.stitch(dataPath, csvFn, 'blendMethod', 'none', 'processFunPath', {funFn});

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(out(1 : 100, 1 : 120, 1 : 32), gt - 99);
        end

        function xcorrRecoversPerturbedTilePosition(testCase)
            % the second tile is listed 4 voxels off in x; registration must
            % recover the true position. Even shifts avoid the +/-1 rounding
            % of the global position solver for two tiles.
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            gt = testCase.randomVolume([64, 104, 32]);
            stageXYZ = [0, 0, 0; 44 * px, 0, 0];
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 0, 40, 0], [64, 64, 32], [], stageXYZ, [0, 0, 0; 1, 0, 0]);

            testCase.stitch(dataPath, csvFn, 'xcorrShift', true, 'xcorrDownsample', [1, 1, 1], ...
                'xyMaxOffset', 20, 'zMaxOffset', 10, 'blendMethod', 'none');

            xcorrInfo = load(sprintf('%smatlab_stitch/xcorr/%s.mat', dataPath, testCase.StitchName));
            testCase.verifyEqual(round(xcorrInfo.dxyz_shift / px), [0, 0, 0; -4, 0, 0], ...
                'Registration did not recover the true tile offset.');

            out = testCase.readResult(dataPath);
            % padding of 2 * #tiles voxels per axis on both sides
            testCase.verifyEqual(size(out), [64 + 2 * 2 + 2, 104 + 2 * 4 + 2, 32 + 2 * 2 + 2]);
            testCase.verifyEqual(out(3 : 66, 5 : 108, 3 : 34), gt);
        end

        function primaryChannelRegistrationAppliesToOtherChannels(testCase)
            % xcorrMode 'primaryFirst': the shift found in the primary channel
            % (CamA) is used for CamB, whose overlap is featureless so it
            % could not be registered on its own
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            gtA = testCase.randomVolume([64, 104, 32]);
            gtB = testCase.randomVolume([64, 104, 32]);
            gtB(:, 41 : 64, :) = 30000;
            stageXYZ = [0, 0, 0; 44 * px, 0, 0];
            tileIdx = [0, 0, 0; 1, 0, 0];
            testCase.writeTiles(dataPath, gtA, [0, 0, 0; 0, 40, 0], [64, 64, 32], [], stageXYZ, tileIdx, 'cam', 'A');
            csvFn = testCase.writeTiles(dataPath, gtB, [0, 0, 0; 0, 40, 0], [64, 64, 32], [], stageXYZ, tileIdx, 'cam', 'B');

            testCase.stitch(dataPath, csvFn, 'channelPatterns', {'CamA_ch0', 'CamB_ch0'}, 'primaryCh', 'CamA_ch0', ...
                'xcorrShift', true, 'xcorrDownsample', [1, 1, 1], 'xyMaxOffset', 20, 'zMaxOffset', 10, 'blendMethod', 'none');

            outA = testCase.readResult(dataPath);
            outB = testCase.readResult(dataPath, 'matlab_stitch', strrep(testCase.StitchName, 'CamA', 'CamB'));
            testCase.verifyEqual(outA(3 : 66, 5 : 108, 3 : 34), gtA);
            testCase.verifyEqual(outB(3 : 66, 5 : 108, 3 : 34), gtB, ...
                'The secondary channel was not stitched with the primary channel registration.');
        end

        function largeFileStitching(testCase)
            % large-file stitching used to fail making the MIP slabs (bim undefined
            % in XR_MIP_zarr) and registering them (three slips in
            % cross_correlation_registration_3d_mip_slabs from a86c36d)
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);

            testCase.stitch(dataPath, csvFn, 'blendMethod', 'feather', 'largeFile', true, 'poolSize', [2, 2, 2], ...
                'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16], 'xcorrShift', true, ...
                'xcorrDownsample', [1, 1, 1], 'xyMaxOffset', 20, 'zMaxOffset', 10);

            out = testCase.readResult(dataPath);
            % padding of 2 * #tiles voxels per axis on both sides
            testCase.verifyEqual(size(out), [100 + 2 * 4 + 2, 120 + 2 * 4 + 2, 32 + 2 * 2 + 2]);
            testCase.verifyEqual(out(5 : 104, 5 : 124, 3 : 34), gt);
        end

        function largeFileRegistrationRecoversPerturbedTilePosition(testCase)
            % large-file registration built the task of each tile pair with
            % the tile and box at the pair's index instead of its first tile,
            % so beyond two tiles it registered the wrong tiles. The third
            % tile is listed 2 voxels off in x.
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            gt = testCase.randomVolume([100, 120, 32]);
            stageXYZ = [0, 0, 0; 0, 36 * px, 0; 50 * px, 0, 0; 48 * px, 36 * px, 0];
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32], ...
                [], stageXYZ, [0, 0, 0; 0, 1, 0; 1, 0, 0; 1, 1, 0]);

            testCase.stitch(dataPath, csvFn, 'blendMethod', 'feather', 'largeFile', true, 'poolSize', [2, 2, 2], ...
                'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16], 'xcorrShift', true, ...
                'xcorrDownsample', [1, 1, 1], 'xyMaxOffset', 20, 'zMaxOffset', 10);

            xcorrInfo = load(sprintf('%smatlab_stitch/xcorr/%s.mat', dataPath, testCase.StitchName));
            % one tile moves 2 voxels back in x (the rows follow the tile
            % order of the stitching, x first)
            shifts = round(xcorrInfo.dxyz_shift / px);
            testCase.verifyEqual(sortrows(shifts), [-2, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0], ...
                sprintf('Registration did not recover the true tile offset: %s.', mat2str(shifts)));
            out = testCase.readResult(dataPath);
            testCase.verifyEqual(out(5 : 104, 5 : 124, 3 : 34), gt);
        end

        function largeFileWithoutRegistration(testCase)
            % the distance maps use MIP slabs, which used to be made only by the
            % registration (xcorrShift)
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32]);

            testCase.stitch(dataPath, csvFn, 'blendMethod', 'feather', 'largeFile', true, 'poolSize', [2, 2, 2], ...
                'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16]);

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(out(1 : 100, 1 : 120, 1 : 32), gt);
        end

        function noneBlendTakesHalfOverlapFromEachTile(testCase)
            % 'none' blending used to take the whole overlap from the lower index tile
            % (regressions from d9a099b and 4a78ade). The split axis is the first in
            % halfOrder (default [3, 2, 1], xyz) with more than one voxel of overlap,
            % z for these tiles, so this checks only that half of the overlap comes
            % from each tile.
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([64, 104, 32], [1000, 50000]);
            tiles = struct('im', {gt(:, 1 : 64, :), gt(:, 41 : 104, :) + 5000}, 'idx', {[0, 0, 0], [1, 0, 0]}, ...
                'xyz', {[0, 0, 0], [40 * testCase.XYPixelSize, 0, 0]});
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);

            testCase.stitch(dataPath, csvFn, 'blendMethod', 'none');

            out = testCase.readResult(dataPath);
            out = out(1 : 64, 1 : 104, 1 : 32);
            testCase.verifyEqual(out(:, 1 : 40, :), gt(:, 1 : 40, :), 'Tile 1 outside the overlap.');
            testCase.verifyEqual(out(:, 65 : 104, :), gt(:, 65 : 104, :) + 5000, 'Tile 2 outside the overlap.');
            % overlap columns 41-64: tile 2 is offset by 5000, so each voxel
            % shows which tile it came from
            overlap = out(:, 41 : 64, :);
            fromTile1 = overlap == gt(:, 41 : 64, :);
            fromTile2 = overlap == gt(:, 41 : 64, :) + 5000;
            testCase.verifyTrue(all(fromTile1 | fromTile2, 'all'), 'The overlap has voxels from neither tile.');
            testCase.verifyEqual(nnz(fromTile1) / numel(overlap), 0.5, 'AbsTol', 0.1, ...
                'The overlap is not split between the two tiles.');
        end

        function blendingSplitsTheOverlapInQuarters(testCase)
            % 'max', 'mean' and 'median' blend the middle half of an overlap
            % and take one tile in each outer quarter (overlap type 'half';
            % 'none' takes half from each tile). The split axis is z for these
            % tiles (see noneBlendTakesHalfOverlapFromEachTile). Tile 2 is
            % offset by 5000, so each voxel shows which tiles it came from.
            gt = testCase.randomVolume([64, 104, 32], [1000, 50000]);
            tiles = struct('im', {gt(:, 1 : 64, :), gt(:, 41 : 104, :) + 5000}, 'idx', {[0, 0, 0], [1, 0, 0]}, ...
                'xyz', {[0, 0, 0], [40 * testCase.XYPixelSize, 0, 0]});
            for method = {'max', 'mean', 'median'}
                dataPath = testCase.createTempFolder();
                [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);
                testCase.stitch(dataPath, csvFn, 'blendMethod', method{1});
                out = testCase.readResult(dataPath);
                overlap = out(1 : 64, 41 : 64, 1 : 32);
                tile1 = gt(:, 41 : 64, :);
                fractions = [nnz(overlap == tile1), nnz(overlap == tile1 + 2500), nnz(overlap == tile1 + 5000)] / numel(overlap);
                if strcmp(method{1}, 'max')
                    expected = [0.25, 0, 0.75];
                else
                    expected = [0.25, 0.5, 0.25];
                end
                testCase.verifyEqual(fractions, expected, 'AbsTol', 0.01, sprintf(['%s: the fractions of the overlap ', ...
                    'from tile 1, both and tile 2 are wrong.'], method{1}));
            end
        end

        function blendingWithoutTheMexGivesTheSameResult(testCase, tileBlendMethod)
            % tile_blending_3d_mex blends the batches of 'none', 'max', 'mean'
            % and 'median' in one pass; without it (a stub that errors) the
            % MATLAB code must give the same volume and MIP. The tiles differ
            % in the overlaps (tile k is offset by 1000 * k) and have zero
            % regions, and the batches hold 1 to 4 tiles.
            testCase.assumeMexAvailable('tile_blending_3d_mex');
            gt = testCase.randomVolume([100, 120, 32], [1000, 50000]);
            starts = [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0];
            tiles = struct('im', cell(4, 1), 'idx', [], 'xyz', []);
            for i = 1 : 4
                s = starts(i, :);
                im = gt(s(1) + 1 : s(1) + 64, s(2) + 1 : s(2) + 72, :) + 1000 * i;
                im(1 : 16, 1 : 24, :) = 0;
                im(randperm(numel(im), 3000)) = 0;
                tiles(i).im = im;
                tiles(i).idx = [s(2) > 0, s(1) > 0, 0];
                tiles(i).xyz = s([2, 1, 3]) * testCase.XYPixelSize;
            end
            opts = {'blendMethod', tileBlendMethod, 'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16]};
            outs = cell(2, 1);
            mips = cell(2, 1);
            for k = 1 : 2
                dataPath = testCase.createTempFolder();
                [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);
                if k == 2
                    testCase.applyFixture(pktest.FunctionStubFixture('tile_blending_3d_mex'));
                    testCase.expectMexFallback('tile_blending_3d_mex');
                end
                evalc('testCase.stitch(dataPath, csvFn, opts{:})');
                outs{k} = testCase.readResult(dataPath);
                mips{k} = readtiff(sprintf('%smatlab_stitch/MIPs/%s_MIP_z.tif', dataPath, testCase.StitchName));
            end
            testCase.verifyEqual(outs{1}, outs{2}, 'The volume differs from that of the MATLAB code.');
            testCase.verifyEqual(mips{1}, mips{2}, 'The MIP differs from that of the MATLAB code.');
            testCase.verifyEqual(mips{1}, max(outs{1}, [], 3), 'The z MIP does not match the stitched volume.');
        end

        function featherBlendingWithThinTiles(testCase, thinTileDepth)
            % the feather z window used to give the first slice weight 0 for 10-29
            % slices and index 0 (error) for fewer than 10
            nz = thinTileDepth;
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, nz]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, nz]);

            testCase.stitch(dataPath, csvFn, 'blendMethod', 'feather');

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(out(1 : 100, 1 : 120, 1 : nz), gt, ...
                sprintf('Feather blending of %d-slice tiles does not reproduce the ground truth.', nz));
        end

        function tilesOfDifferentSizesGetTheirOwnDistanceMaps(testCase)
            % tiles of one location with different sizes get a distance map
            % each. Their overlap regions were cropped with a 3 x 2 bounding
            % box, which crop3d_mex rejected (the error was printed and a
            % MATLAB copy used instead).
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            starts = [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0];
            % tile 1 is narrower; tile 3 still covers the columns it lacks
            sizes = [64, 68, 32; 64, 72, 32; 64, 72, 32; 64, 72, 32];
            tiles = struct('im', cell(4, 1), 'idx', [], 'xyz', []);
            for i = 1 : 4
                s = starts(i, :);
                tiles(i).im = gt(s(1) + 1 : s(1) + sizes(i, 1), s(2) + 1 : s(2) + sizes(i, 2), s(3) + 1 : s(3) + sizes(i, 3));
                tiles(i).idx = [s(2) > 0, s(1) > 0, 0];
                tiles(i).xyz = s([2, 1, 3]) * testCase.XYPixelSize;
            end
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);

            out = evalc('testCase.stitch(dataPath, csvFn, ''blendMethod'', ''feather'')');

            testCase.verifyNumElements(dir(sprintf('%smatlab_stitch/imdist/**/*.zarr', dataPath)), 4, ...
                'Each tile should have its own distance map.');
            testCase.verifyEmpty(strfind(out, 'Input range for bbox'), 'Cropping the distance map regions failed.');
            result = testCase.readResult(dataPath);
            testCase.verifyEqual(result(1 : 100, 1 : 120, 1 : 32), gt);
        end

        function localShiftMethod(testCase)
            % 'local' used to read the shift table in its layout before 812679e
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            gt = testCase.randomVolume([64, 104, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 0, 40, 0], [64, 64, 32], [], ...
                [0, 0, 0; 44 * px, 0, 0], [0, 0, 0; 1, 0, 0]);

            testCase.stitch(dataPath, csvFn, 'xcorrShift', true, 'xcorrDownsample', [1, 1, 1], ...
                'xyMaxOffset', 20, 'zMaxOffset', 10, 'blendMethod', 'none', 'shiftMethod', 'local');

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(out(3 : 66, 5 : 108, 3 : 34), gt);
        end

        function flippedTilesStayWithTheirTilesWhenSorted(testCase, flipDSR)
            % the stitching sorts the tiles by their indices (zigzag
            % listings), but the flags of flipped tiles (flippedTile, from
            % the settings files) stayed in the listed order and so went to
            % other tiles: the stitched volume depended on the order of the
            % list. The frame is called directly, as the wrapper takes the
            % flags only from settings files; each run has tiles of its own.
            px = testCase.XYPixelSize;
            gt = testCase.randomVolume([64, 108, 32]);
            tiles = testCase.cutTiles(gt, [0, 0, 0; 0, 36, 0], [64, 72, 32]);
            idx = [cat(1, tiles.idx), zeros(2, 1)];
            runs = {[1, 2], [false, true]; [2, 1], [true, false]; [1, 2], [false, false]};
            out = cell(1, 3);
            for k = 1 : 3
                dataPath = testCase.createTempFolder();
                fsns = pktest.writeStitchTiles(dataPath, tiles);
                order = runs{k, 1};
                fns = strcat(dataPath, fsns(order), '.tif');
                xyz = cat(1, tiles(order).xyz);
                mkdir([dataPath, 'res']);
                evalc(['XR_stitching_frame_zarr_dev_v1(fns, xyz, ''ResultPath'', [dataPath, ''res/out.zarr''], ', ...
                    '''tileIdx'', idx(order, :), ''flippedTile'', runs{k, 2}, ''DSR'', flipDSR, ''resampleType'', ''isotropic'', ', ...
                    '''xyPixelSize'', px, ''dz'', 0.3, ''axisOrder'', ''x,y,z'', ''xcorrShift'', false, ', ...
                    '''parseCluster'', false, ''BlendMethod'', ''none'')']);
                out{k} = readzarr([dataPath, 'res/out.zarr']);
            end
            testCase.verifyEqual(out{2}, out{1}, 'The order of the list changed the stitched volume.');
            testCase.verifyNotEqual(out{3}, out{1}, 'The flag had no effect.');
        end

        function zxyTilesOfOneXPlaneAreNot2D(testCase)
            % the first 2D check read the size at z's place in the data order
            % (imSize(data_order_mat == 3)) instead of z's dimension, which
            % differ for 'y,z,x' and 'z,x,y': 3D tiles in 'z,x,y' order with
            % one x plane passed for 2D tiles, and only their first z plane
            % was stitched
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            gt = permute(testCase.randomVolume([100, 1, 8]), [3, 2, 1]);
            tiles = struct('im', {gt(:, :, 1 : 60), gt(:, :, 41 : 100)}, 'idx', {[0, 0, 0], [0, 1, 0]}, ...
                'xyz', {[0, 0, 0], [0, 40 * px, 0]});
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);

            testCase.stitch(dataPath, csvFn, 'dataOrder', 'z,x,y');

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(size(out), [10, 3, 102]);
            testCase.verifyEqual(out(1 : 8, 1, 1 : 100), gt);
        end

        function resampleFactorIsInTheTilesAxisOrder(testCase)
            % resampleFactor is in yxz order, the order in which tiles are
            % resampled (and DSR tiles were deskewed with it); the stitching
            % took it as x, y, z for the voxel sizes, so tiles whose y and x
            % factors differ were placed with the x and y voxel sizes
            % swapped. DSR tiles with y voxels of 2 px, as a factor of
            % [2, 1, 1] makes them.
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            gt = testCase.randomVolume([80, 108, 32]);
            csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 0, 36, 0; 30, 0, 0; 30, 36, 0], [50, 72, 32], [2 * px, px, px]);

            testCase.stitch(dataPath, csvFn, 'resampleType', 'given', 'resampleFactor', [2, 1, 1]);

            out = testCase.readResult(dataPath);
            testCase.verifyEqual(size(out), [82, 110, 34]);
            testCase.verifyEqual(out(1 : 80, 1 : 108, 1 : 32), gt);
        end
    end

    methods (Test, TestTags = {'Integration', 'Slow'})
        function localWorkersGiveTheSameResult(testCase)
            % with local workers, the tiles are registered and blended (in
            % several tasks) on a local process pool; the result is that of
            % one task at a time in this MATLAB
            testCase.assumeTrue(license('test', 'Distrib_Computing_Toolbox') && ~isempty(which('gcp')), ...
                'Needs the Parallel Computing Toolbox.');
            if isempty(gcp('nocreate'))
                testCase.addTeardown(@() delete(gcp('nocreate')));
            end
            % memory for the workers whatever the machine
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 1024));
            px = testCase.XYPixelSize;
            gt = testCase.randomVolume([100, 120, 32]);
            % the third tile is listed 2 voxels off in x
            stageXYZ = [0, 0, 0; 0, 36 * px, 0; 50 * px, 0, 0; 48 * px, 36 * px, 0];
            opts = {'blendMethod', 'feather', 'xcorrShift', true, 'xcorrDownsample', [1, 1, 1], 'xyMaxOffset', 10, ...
                'zMaxOffset', 5, 'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16]};
            outs = cell(2, 1);
            mips = cell(2, 1);
            for k = 1 : 2
                dataPath = testCase.createTempFolder();
                csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], ...
                    [64, 72, 32], [], stageXYZ, [0, 0, 0; 0, 1, 0; 1, 0, 0; 1, 1, 0]);
                workerNum = 3 - k;
                out = evalc('testCase.stitch(dataPath, csvFn, opts{:}, ''localWorkerNum'', workerNum)');
                outs{k} = testCase.readResult(dataPath);
                mips{k} = readtiff(sprintf('%smatlab_stitch/MIPs/%s_MIP_z.tif', dataPath, testCase.StitchName));
                info = load(sprintf('%smatlab_stitch/stitchInfo/%s.mat', dataPath, testCase.StitchName), 'PerBlockInfoFullpaths');
                if k == 1
                    testCase.verifyGreaterThan(numel(info.PerBlockInfoFullpaths), 1, 'The blending ran in one task.');
                    testCase.verifyNumElements(strfind(out, 'on the local pool'), 2, ...
                        'The registration and the blending tasks did not run on the local pool.');
                else
                    testCase.verifyNumElements(info.PerBlockInfoFullpaths, 1);
                end
            end
            testCase.verifyEqual(outs{1}, outs{2}, 'The result with local workers differs.');
            testCase.verifyEqual(mips{1}, mips{2}, 'The MIP with local workers differs.');
        end
    end

    methods (Access = private)
        function [gt] = randomVolume(~, sz, intRange)
            % nonzero random uint16 values (0 means no data in stitching)
            if nargin < 3
                intRange = [1000, 60000];
            end
            gt = uint16(randi(intRange, sz));
        end

        function [dataPath, csvFn, resultFns] = writeTwoVolumes(testCase, n)
            % two (or n) time points of the 2 x 2 tile layout, and the paths
            % of their stitched results
            if nargin < 2
                n = 2;
            end
            dataPath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            for iter = 0 : n - 1
                csvFn = testCase.writeTiles(dataPath, gt, [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0], [64, 72, 32], ...
                    [], [], [], 'iter', iter);
            end
            resultFns = arrayfun(@(x) sprintf('%smatlab_stitch/%s.zarr', dataPath, ...
                strrep(testCase.StitchName, 'Iter_0000', sprintf('Iter_%04d', x))), 0 : n - 1, 'unif', 0);
        end

        function [dataPath, csvFn, stagePath] = writeTilesToCome(testCase)
            % the raw tiles of the 2 x 2 tile layout, and their processed
            % (DSR) tiles in stagePath/DSR/, to be moved to dataPath/DSR/
            dataPath = testCase.createTempFolder();
            stagePath = testCase.createTempFolder();
            gt = testCase.randomVolume([100, 120, 32]);
            starts = [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0];
            csvFn = testCase.writeTiles(dataPath, gt, starts, [64, 72, 32]);
            testCase.writeTiles(stagePath, gt, starts, [64, 72, 32], [], [], [], 'tileDir', 'DSR');
            mkdir([dataPath, 'DSR']);
        end

        function [] = atCall(~, counterFn, n, fcn)
            % count the calls in counterFn; at the n-th, run fcn (n a
            % vector: at the n(i)-th, run fcn{i})
            k = 1;
            if isfile(counterFn)
                k = str2double(fileread(counterFn)) + 1;
            end
            fid = fopen(counterFn, 'w');
            fprintf(fid, '%d', k);
            fclose(fid);
            if ~iscell(fcn)
                fcn = {fcn};
            end
            if any(k == n)
                fcn{k == n}();
            end
        end

        function [csvFn] = writeTiles(testCase, dataPath, gt, starts, tileSize, pixelSizes, stageXYZ, tileIdx, varargin)
            % write the tiles of cutTiles
            if nargin < 6
                pixelSizes = [];
            end
            if nargin < 7
                stageXYZ = [];
            end
            if nargin < 8
                tileIdx = [];
            end
            tiles = testCase.cutTiles(gt, starts, tileSize, pixelSizes, stageXYZ, tileIdx);
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles, varargin{:});
        end

        function [tiles] = cutTiles(testCase, gt, starts, tileSize, pixelSizes, stageXYZ, tileIdx)
            % cut tiles of tileSize ([y, x, z]) at 0-based voxel offsets
            % starts ([y0, x0, z0] per row) from the ground truth. By default
            % the stage coordinates are the offsets times the voxel size
            % and the tile indices follow the positions.
            if nargin < 5 || isempty(pixelSizes)
                pixelSizes = testCase.XYPixelSize * [1, 1, 1];
            end
            nT = size(starts, 1);
            if nargin < 6 || isempty(stageXYZ)
                stageXYZ = starts(:, [2, 1, 3]) .* pixelSizes([2, 1, 3]);
            end
            if nargin < 7 || isempty(tileIdx)
                xyzStarts = starts(:, [2, 1, 3]);
                tileIdx = zeros(nT, 3);
                for i = 1 : 3
                    [~, ~, tileIdx(:, i)] = unique(xyzStarts(:, i));
                end
                tileIdx = tileIdx - 1;
            end
            tiles = struct('im', cell(nT, 1), 'idx', [], 'xyz', []);
            for i = 1 : nT
                s = starts(i, :);
                if ismatrix(gt)
                    tiles(i).im = gt(s(1) + 1 : s(1) + tileSize(1), s(2) + 1 : s(2) + tileSize(2));
                else
                    tiles(i).im = gt(s(1) + 1 : s(1) + tileSize(1), s(2) + 1 : s(2) + tileSize(2), s(3) + 1 : s(3) + tileSize(3));
                end
                tiles(i).idx = tileIdx(i, :);
                tiles(i).xyz = stageXYZ(i, :);
            end
        end

        function [] = stitch(testCase, dataPath, csvFn, varargin)
            % defaults for exact DSR stitching; later name-value pairs win
            XR_matlab_stitching_wrapper(dataPath, csvFn, 'DSR', true, 'resampleType', 'isotropic', ...
                'xyPixelSize', testCase.XYPixelSize, 'axisOrder', 'x,y,z', 'xcorrShift', false, ...
                'channelPatterns', {'CamA_ch0'}, 'primaryCh', 'CamA_ch0', 'parseCluster', false, varargin{:});
        end

        function [im] = cropToData(~, im)
            % crop to the bounding box of the nonzero voxels
            ranges = cell(1, ndims(im));
            for d = 1 : ndims(im)
                inds = find(any(im ~= 0, setdiff(1 : ndims(im), d)));
                ranges{d} = inds(1) : inds(end);
            end
            im = im(ranges{:});
        end

        function [out] = readResult(testCase, dataPath, resultDirName, name)
            if nargin < 3
                resultDirName = 'matlab_stitch';
            end
            if nargin < 4
                name = testCase.StitchName;
            end
            out = readzarr(sprintf('%s%s/%s.zarr', dataPath, resultDirName, name));
        end
    end
end
