classdef ImageListGenerationTest < pktest.TestCase
    % XR_generate_image_list_wrapper: image list csv files (stage
    % coordinates for stitching) from tile lists, tile positions, and
    % encoder slice lists. The tile_position cases follow
    % demo_phase_and_2photon_stitching.m and check that the generated list
    % stitches synthetic tiles back to the ground truth.

    properties (Constant)
        XYPixelSize = 0.108
    end

    methods (Test, TestTags = {'Unit'})
        function tileListCoordinates(testCase)
            dataPath = testCase.createTempFolder();
            fns = {'tile_a.tif'; 'tile_b.tif'; 'tile_c.tif'; 'tile_d.tif'};
            % [t, c, x, y, z]
            tileIndices = [0, 0, 1, 0, 2; 0, 0, 2, 0, 2; 0, 0, 1, 1, 2; 0, 1, 2, 1, 3];
            tileInterval = [50, 40, 10];

            XR_generate_image_list_wrapper(dataPath, 'tile_list', 'tileFilenames', fns, ...
                'tileIndices', tileIndices, 'tileInterval', tileInterval);

            tab = readtable([dataPath, 'ImageList_from_tile_list.csv'], 'Delimiter', 'comma');
            testCase.verifyEqual(tab.Properties.VariableNames, {'Filepath', 'Filename', 'mappedFilename', ...
                'StageX_um_', 'StageY_um_', 'StageZ_um_', 'ObjectiveX_um_', 'ObjectiveY_um_', 'ObjectiveZ_um_'});
            testCase.verifyEqual(tab.Filename, fns);
            % positions relative to the minimum tile index times the interval
            expected = (tileIndices(:, 3 : 5) - min(tileIndices(:, 3 : 5))) .* tileInterval;
            testCase.verifyEqual([tab.StageX_um_, tab.StageY_um_, tab.StageZ_um_], expected, 'AbsTol', 1e-9);
            testCase.verifyEqual([tab.ObjectiveX_um_, tab.ObjectiveY_um_, tab.ObjectiveZ_um_], expected, 'AbsTol', 1e-9);
            testCase.verifyEqual(tab.mappedFilename{4}, ...
                'Scan_Iter_0000_CamA_ch1_CAM1_stack0000_488nm_000000msec_000000msecAbs_002x_001y_003z_0000t.tif');
        end

        function encoderSliceListCoordinates(testCase)
            % stage coordinates come from the first MCS2 position of each
            % tile in the SliceList csv files
            dataPath = testCase.createTempFolder();
            mcs2 = [490.05, -8052.12, -5473.28; 490.04, -7862.72, -5473.27];
            fns = testCase.writeSliceLists(dataPath, mcs2);

            XR_generate_image_list_wrapper(dataPath, 'encoder', 'dz', 0.2, 'channelPatterns', {'CamA_ch0'});

            testCase.verifyEncoderImageList(dataPath, fns, mcs2);
            testCase.verifyTrue(exist([dataPath, 'encoder_info/coordinates_info.mat'], 'file') == 2);
        end

        function encoderListsForMultipleDatasets(testCase)
            % one image list per dataset folder (regression for ef90ba3:
            % a cell of data paths was passed on as is and errored)
            dataPaths = {testCase.createTempFolder(), testCase.createTempFolder()};
            mcs2 = {[10, 20, 30; 10, 220, 30], [-5, 0, 2; -5, 150, 2]};
            fns = cell(1, 2);
            for d = 1 : 2
                fns{d} = testCase.writeSliceLists(dataPaths{d}, mcs2{d});
            end

            XR_generate_image_list_wrapper(dataPaths, 'encoder', 'dz', 0.2, 'channelPatterns', {'CamA_ch0'});

            for d = 1 : 2
                testCase.verifyEncoderImageList(dataPaths{d}, fns{d}, mcs2{d});
            end
        end
    end

    methods (Test, TestTags = {'Integration'})
        function tilePositions2DPhase(testCase)
            % 2D tiles with standard names, objective/IO scan and overlap in
            % pixels; the generated list stitches the tiles exactly
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            gt = uint16(randi([1000, 60000], 120, 120));
            step = 56;
            tiles = struct('im', {}, 'idx', {}, 'xyz', {});
            for xi = 0 : 1
                for yi = 0 : 1
                    tiles(end + 1) = struct('im', gt(yi * step + (1 : 64), xi * step + (1 : 64)), ...
                        'idx', [xi, yi, 0], 'xyz', [0, 0, 0]); %#ok<AGROW>
                end
            end
            pktest.writeStitchTiles(dataPath, tiles, 'cam', 'B');

            XR_generate_image_list_wrapper(dataPath, 'tile_position', 'xyPixelSize', px, ...
                'channelPatterns', {'Scan'}, 'overlapSize', [8, 8, 1], 'overlapSizeType', 'pixel', ...
                'objectiveScan', false, 'IOScan', true);

            csvFn = [dataPath, 'ImageList_from_tile_positions.csv'];
            tab = readtable(csvFn, 'Delimiter', 'comma');
            idx = [cellfun(@(x) str2double(regexp(x, '(\d+)x_', 'tokens', 'once')), tab.Filename), ...
                cellfun(@(x) str2double(regexp(x, '(\d+)y_', 'tokens', 'once')), tab.Filename)];
            testCase.verifyEqual([tab.StageX_um_, tab.StageY_um_, tab.StageZ_um_], [idx * step * px, zeros(4, 1)], 'AbsTol', 1e-9);

            XR_matlab_stitching_wrapper(dataPath, csvFn, 'DS', false, 'DSR', false, 'IOScan', true, ...
                'xyPixelSize', px, 'resampleType', 'isotropic', 'axisOrder', 'x,y,z', 'channelPatterns', {'CamB_ch0'}, ...
                'xcorrShift', false, 'blendMethod', 'feather', 'blockSize', [64, 64, 1], 'parseCluster', false);

            dir_info = dir([dataPath, 'matlab_stitch/*.zarr']);
            testCase.assertNumElements(dir_info, 1);
            out = readzarr([dataPath, 'matlab_stitch/', dir_info(1).name]);
            testCase.verifyEqual(out(1 : 120, 1 : 120), gt);
        end

        function tilePositions3DNonStandardNames(testCase)
            % 2-photon style tiles (Tile_ch0_000x_000y_000z_0000t.tif) with
            % overlap in um and flipped x and z axes ('-x,y,-z'), as in
            % demo_phase_and_2photon_stitching.m
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            dz = 0.5;
            % 40 x 40 x 36 tiles (>= 30 slices for feather blending, see
            % StitchingTest/featherBlendingWithThinTiles) with 8 voxel
            % overlap in x/y and 4 slices in z
            gt = uint16(randi([1000, 60000], 40, 72, 68));
            overlapUm = [8 * px, 8 * px, 4 * dz];
            stepVox = [40 - 8, 36 - 4];
            for xi = 0 : 1
                for zi = 0 : 1
                    % with -x and -z the tile with the larger index is placed first
                    x0 = (1 - xi) * stepVox(1);
                    z0 = (1 - zi) * stepVox(2);
                    writetiff(gt(:, x0 + (1 : 40), z0 + (1 : 36)), sprintf('%sTile_ch0_%03dx_000y_%03dz_0000t.tif', dataPath, xi, zi));
                end
            end

            XR_generate_image_list_wrapper(dataPath, 'tile_position', 'xyPixelSize', px, 'dz', dz, ...
                'channelPatterns', {'Tile'}, 'tilePatterns', {'0000t', 'ch0', '000x', '000y', '000z'}, ...
                'overlapSize', overlapUm, 'overlapSizeType', 'um', 'objectiveScan', false, 'IOScan', true);

            csvFn = [dataPath, 'ImageList_from_tile_positions.csv'];
            tab = readtable(csvFn, 'Delimiter', 'comma');
            testCase.verifyTrue(ismember('mappedFilename', tab.Properties.VariableNames));
            testCase.verifyTrue(all(startsWith(tab.mappedFilename, 'Scan_Iter_0000_ch0_CAM1_stack0000_488nm_')));

            XR_matlab_stitching_wrapper(dataPath, csvFn, 'DS', false, 'DSR', false, 'IOScan', true, ...
                'xyPixelSize', px, 'dz', dz, 'resampleType', 'isotropic', 'axisOrder', '-x,y,-z', ...
                'channelPatterns', {'ch0'}, 'primaryCh', 'ch0', 'xcorrShift', false, 'blendMethod', 'feather', ...
                'save16bit', false, 'parseCluster', false);

            dir_info = dir([dataPath, 'matlab_stitch/*.zarr']);
            testCase.assertNumElements(dir_info, 1);
            out = readzarr([dataPath, 'matlab_stitch/', dir_info(1).name]);
            testCase.verifyEqual(size(out), [42, 74, 70]);
            testCase.verifyEqual(out(1 : 40, 1 : 72, 1 : 68), gt);
        end

        function tileNamesWithoutStandardChannel(testCase)
            % tile names that match neither LLSM pattern (no '_ch[0-9]_'):
            % the list must still be generated (regression for 881d319,
            % 'expression' was undefined), and the stitching channel, which
            % is only in mappedFilename, must select the process function
            % (regression for abff2bc, it was silently skipped)
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            dz = 0.5;
            gt = uint16(randi([1000, 50000], 40, 72, 32));
            for xi = 0 : 1
                writetiff(gt(:, xi * 32 + (1 : 40), :), sprintf('%sTile_c0_%03dx_000y_000z_0000t.tif', dataPath, xi));
            end

            XR_generate_image_list_wrapper(dataPath, 'tile_position', 'xyPixelSize', px, 'dz', dz, ...
                'channelPatterns', {'Tile'}, 'tilePatterns', {'0000t', 'c0', '000x', '000y', '000z'}, ...
                'overlapSize', [8, 8, 0], 'overlapSizeType', 'pixel', 'IOScan', true);

            csvFn = [dataPath, 'ImageList_from_tile_positions.csv'];
            tab = readtable(csvFn, 'Delimiter', 'comma');
            testCase.verifyEqual(sort(tab.StageX_um_), [0; 32 * px], 'AbsTol', 1e-9);

            usrFun = "@(x)x+1000";
            funFn = [dataPath, 'processFunction.mat'];
            save('-v7.3', funFn, 'usrFun');
            XR_matlab_stitching_wrapper(dataPath, csvFn, 'DS', false, 'DSR', false, 'IOScan', true, ...
                'xyPixelSize', px, 'dz', dz, 'resampleType', 'isotropic', 'axisOrder', 'x,y,z', ...
                'channelPatterns', {'ch0'}, 'primaryCh', 'ch0', 'xcorrShift', false, 'blendMethod', 'none', ...
                'processFunPath', {funFn}, 'parseCluster', false);

            dir_info = dir([dataPath, 'matlab_stitch/*.zarr']);
            testCase.assertNumElements(dir_info, 1);
            out = readzarr([dataPath, 'matlab_stitch/', dir_info(1).name]);
            testCase.verifyEqual(out(1 : 40, 1 : 72, 1 : 32), gt + 1000, ...
                'The process function was not applied to tiles whose channel is only in mappedFilename.');
        end
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function tilePositionStepsForNonSquareTiles(testCase)
            % KNOWN ISSUE: stitch_generate_imagelist_from_tile_positions.m
            % computes positions as (tileIndices - ts) .* sz_1([2, 1, 3])
            % with sz_1 already in x, y, z order, so the x step uses the tile
            % height and the y step the tile width. Square tiles hide it.
            % Fix: use sz_1 without reordering.
            dataPath = testCase.createTempFolder();
            px = testCase.XYPixelSize;
            % 48 x 64 (y x) tiles in a 2 x 2 grid
            for xi = 0 : 1
                for yi = 0 : 1
                    writetiff(ones(48, 64, 'uint16'), sprintf(['%sScan_Iter_0000_CamB_ch0_CAM1_stack0000_642nm_', ...
                        '0000000msec_0000000000msecAbs_%03dx_%03dy_000z_0000t.tif'], dataPath, xi, yi));
                end
            end

            XR_generate_image_list_wrapper(dataPath, 'tile_position', 'xyPixelSize', px, ...
                'channelPatterns', {'Scan'}, 'overlapSize', [8, 8, 1], 'overlapSizeType', 'pixel', 'IOScan', true);

            tab = readtable([dataPath, 'ImageList_from_tile_positions.csv'], 'Delimiter', 'comma');
            testCase.verifyEqual(max(tab.StageX_um_), (64 - 8) * px, 'AbsTol', 1e-9, 'Wrong tile step in x.');
            testCase.verifyEqual(max(tab.StageY_um_), (48 - 8) * px, 'AbsTol', 1e-9, 'Wrong tile step in y.');
        end
    end

    methods (Access = private)
        function [fns] = writeSliceLists(~, dataPath, mcs2)
            % one SliceList csv per tile (3 slices, stage moving in x); the
            % image list uses the first MCS2 position of each tile
            fsn = 'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_%010dmsecAbs_000x_%03dy_000z_0000t.tif';
            header = ['Filename,Channel,Is_CamB_,XFPGAPosition_um_,YFPGAPosition_um_,ZFPGAPosition_um_,', ...
                'XMCS2Position_um_,YMCS2Position_um_,ZMCS2Position_um_'];
            nT = size(mcs2, 1);
            fns = cell(nT, 1);
            for i = 1 : nT
                fns{i} = sprintf(fsn, 1000 * i, i - 1);
                fid = fopen(sprintf('%sSliceList_tile%d_part%04d.csv', dataPath, i, i - 1), 'w');
                fprintf(fid, '%s\n', header);
                for s = 1 : 3
                    fprintf(fid, '%s,0,0,%f,%f,%f,%f,%f,%f\n', fns{i}, mcs2(i, 1) + 0.2 * s, mcs2(i, 2), mcs2(i, 3), ...
                        mcs2(i, 1) + 0.2 * (s - 1), mcs2(i, 2), mcs2(i, 3));
                end
                fclose(fid);
            end
        end

        function [] = verifyEncoderImageList(testCase, dataPath, fns, mcs2)
            tab = readtable([dataPath, 'ImageList_from_encoder.csv'], 'Delimiter', 'comma');
            [~, order] = sort(tab.Filename);
            tab = tab(order, :);
            [sortedFns, fnOrder] = sort(fns);
            testCase.verifyEqual(tab.Filename, sortedFns);
            testCase.verifyEqual([tab.StageX_um_, tab.StageY_um_, tab.StageZ_um_], mcs2(fnOrder, :), 'AbsTol', 1e-6);
        end
    end
end
