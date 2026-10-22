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
        blendMethod = {'none', 'feather', 'max', 'mean', 'median'}
        axisOrder = struct('x', 'x,y,z', 'minusX', '-x,y,z')
        thinTileDepth = struct('nz8', 8, 'nz20', 20)
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
            % a check before submitting, then one per round until the job is gone
            testCase.verifyEqual(slurm.callCount('squeue'), 4);
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
    end

    methods (Access = private)
        function [gt] = randomVolume(~, sz, intRange)
            % nonzero random uint16 values (0 means no data in stitching)
            if nargin < 3
                intRange = [1000, 60000];
            end
            gt = uint16(randi(intRange, sz));
        end

        function [csvFn] = writeTiles(testCase, dataPath, gt, starts, tileSize, pixelSizes, stageXYZ, tileIdx, varargin)
            % cut tiles of tileSize ([y, x, z]) at 0-based voxel offsets
            % starts ([y0, x0, z0] per row) from the ground truth. By default
            % the stage coordinates are the offsets times the voxel size
            % and the tile indices follow the positions.
            if nargin < 6 || isempty(pixelSizes)
                pixelSizes = testCase.XYPixelSize * [1, 1, 1];
            end
            nT = size(starts, 1);
            if nargin < 7 || isempty(stageXYZ)
                stageXYZ = starts(:, [2, 1, 3]) .* pixelSizes([2, 1, 3]);
            end
            if nargin < 8 || isempty(tileIdx)
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
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles, varargin{:});
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
