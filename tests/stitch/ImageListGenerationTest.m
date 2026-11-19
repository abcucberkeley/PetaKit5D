classdef ImageListGenerationTest < pktest.TestCase
    % XR_generate_image_list_wrapper: image list csv files (stage
    % coordinates for stitching) from tile lists, tile positions, encoder
    % slice lists, and acquisition databases. The tile_position cases follow
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

        function tileListOfASingleTile(testCase)
            % the minimum tile index of a single tile was the minimum over
            % its axes, so the tile was placed away from the origin
            dataPath = testCase.createTempFolder();
            XR_generate_image_list_wrapper(dataPath, 'tile_list', 'tileFilenames', {'tile_a.tif'}, ...
                'tileIndices', [0, 0, 1, 3, 2], 'tileInterval', [50, 40, 10]);

            tab = readtable([dataPath, 'ImageList_from_tile_list.csv'], 'Delimiter', 'comma');
            testCase.verifyEqual(tab.Filename, {'tile_a.tif'});
            testCase.verifyEqual([tab.StageX_um_, tab.StageY_um_, tab.StageZ_um_], [0, 0, 0]);
        end

        function tilePositionsOfASingleTile(testCase)
            % the mode of the sizes of a single tile (a 1 x 3 row) was a
            % single number, and indexing it as [y, x, z] errored
            dataPath = testCase.createTempFolder();
            writetiff(ones(48, 64, 4, 'uint16'), [dataPath, 'Scan_Iter_0000_CamB_ch0_CAM1_stack0000_642nm_', ...
                '0000000msec_0000000000msecAbs_002x_003y_001z_0000t.tif']);

            XR_generate_image_list_wrapper(dataPath, 'tile_position', 'xyPixelSize', testCase.XYPixelSize, ...
                'channelPatterns', {'Scan'}, 'overlapSize', [8, 8, 1], 'overlapSizeType', 'pixel', 'IOScan', true);

            tab = readtable([dataPath, 'ImageList_from_tile_positions.csv'], 'Delimiter', 'comma');
            testCase.verifyEqual(height(tab), 1);
            testCase.verifyEqual([tab.StageX_um_, tab.StageY_um_, tab.StageZ_um_], [0, 0, 0]);
        end

        function tilePositionsOf2DZarrTiles(testCase)
            % the size of a 2D zarr tile (here made without the z axis
            % createzarr adds by default; the list reads only its header) is
            % [ny, nx], that of a 2D tiff [ny, nx, 1], and indexing it as
            % [y, x, z] errored; the list must be that of the same tiles as
            % tiff files
            px = testCase.XYPixelSize;
            formats = {'tif', 'zarr'};
            stages = cell(1, 2);
            for k = 1 : 2
                dataPath = testCase.createTempFolder();
                % 48 x 64 (y x) tiles in a 2 x 2 grid
                for xi = 0 : 1
                    for yi = 0 : 1
                        fn = sprintf(['%sScan_Iter_0000_CamB_ch0_CAM1_stack0000_642nm_0000000msec_', ...
                            '0000000000msecAbs_%03dx_%03dy_000z_0000t.%s'], dataPath, xi, yi, formats{k});
                        if k == 1
                            writetiff(ones(48, 64, 'uint16'), fn);
                        else
                            createzarr(fn, dataSize=[48, 64], blockSize=[48, 64], expand2dDim=false);
                        end
                    end
                end
                if k == 2
                    testCase.assertEqual(getImageSize(fn), [48, 64], 'The zarr tiles are not 2D.');
                end

                XR_generate_image_list_wrapper(dataPath, 'tile_position', 'xyPixelSize', px, ...
                    'channelPatterns', {'Scan'}, 'overlapSize', [8, 8, 1], 'overlapSizeType', 'pixel', 'IOScan', true, ...
                    'zarrFile', k == 2);

                tab = readtable([dataPath, 'ImageList_from_tile_positions.csv'], 'Delimiter', 'comma');
                testCase.assertEqual(height(tab), 4);
                tab = sortrows(tab, 'Filename');
                stages{k} = [tab.StageX_um_, tab.StageY_um_, tab.StageZ_um_];
            end
            testCase.verifyEqual(stages{2}, stages{1}, 'The 2D zarr tiles are placed unlike the same tiff tiles.');
            testCase.verifyEqual(max(stages{2}(:, 1 : 2)), [64 - 8, 48 - 8] * px, 'AbsTol', 1e-9);
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

        function sqliteViewWithoutMetadataHashReadsTheTables(testCase)
            % ImageListView_NoNull selects JSONsettings_MD5 from
            % EverythingView_NoNull, which joins on it but doesn't return
            % it, so the view fails with "no such column: JSONsettings_MD5"
            % (acquisition databases of 2026-09); the list is then read from
            % the image and slice tables with a warning, and must be the one
            % the view gives where the inner view returns the column
            broken = testCase.createTempFolder();
            fixed = testCase.createTempFolder();
            testCase.writeAcquisitionDb([broken, 'Scan_Iter_0000_slicelist.sqlite3'], 'missing');
            testCase.writeAcquisitionDb([fixed, 'Scan_Iter_0000_slicelist.sqlite3'], 'returned');

            for mapTilename = [true, false]
                testCase.verifyWarning(@() stitch_generate_imagelist_from_sqlite(broken, ...
                    'mapTilename', mapTilename), 'stitch_generate_imagelist_from_sqlite:missingMetadataHashInView');
                testCase.verifyWarningFree(@() stitch_generate_imagelist_from_sqlite(fixed, ...
                    'mapTilename', mapTilename));
                listBroken = fileread([broken, 'ImageList_from_sqlite.csv']);
                listFixed = fileread([fixed, 'ImageList_from_sqlite.csv']);
                testCase.verifyEqual(listBroken, listFixed, ...
                    'The list read from the tables is not the one the view gives.');
            end
            % the images without a metadata row or without slices are left
            % out, as the view's joins leave them out
            tab = readtable([broken, 'ImageList_from_sqlite.csv'], 'Delimiter', 'comma');
            testCase.verifyEqual(height(tab), 8);
            testCase.verifyFalse(any(contains(tab.Filename, {'nohash', 'nometa', 'noslice'})));
        end

        function sqliteViewWithOtherMissingColumnStillFails(testCase)
            % only the missing JSONsettings_MD5 is read around; another
            % missing column (here one whose name only starts with it) is
            % the view's error, without the warning
            dataPath = testCase.createTempFolder();
            testCase.writeAcquisitionDb([dataPath, 'Scan_Iter_0000_slicelist.sqlite3'], 'other');
            testCase.verifyWarningFree(@() testCase.verifyError(@() ...
                stitch_generate_imagelist_from_sqlite(dataPath), ?MException));
            testCase.verifyFalse(exist([dataPath, 'ImageList_from_sqlite.csv'], 'file') == 2);
        end

        function tilePositionStepsForNonSquareTiles(testCase)
            % the x and y steps used to be swapped for non-square tiles
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

        function [] = writeAcquisitionDb(~, fn, hash)
            % a made-up acquisition database: image, slice and metadata
            % tables, and the two views the image list reads. hash: the
            % inner view's JSONsettings_MD5, 'returned' (the view works),
            % 'missing' (joined on, not returned: the view fails as in the
            % databases of 2026-09) or 'other' (the outer view asks for
            % JSONsettings_MD5_old, which nothing returns)
            conn = sqlite(fn, 'create');
            cleanup = onCleanup(@() close(conn));
            execute(conn, ['CREATE TABLE MetadataList (MetadataList_id INTEGER PRIMARY KEY, ', ...
                'JSONsettings_MD5 INTEGER UNIQUE, JSONsettings TEXT)']);
            execute(conn, ['CREATE TABLE ImageList (ImageList_id INTEGER PRIMARY KEY, Filename TEXT UNIQUE, ', ...
                'Camera INTEGER, Channel INTEGER, Laser_nm INTEGER, File_Stack_Index INTEGER, Elapsed_ms REAL, ', ...
                'Absolute_ms INTEGER, Tile_index_X INTEGER, Tile_index_Y INTEGER, Tile_index_Z INTEGER, ', ...
                'Tile_index_T INTEGER, Frame_Type TEXT, JSONsettings_MD5 INTEGER)']);
            execute(conn, ['CREATE TABLE SliceList (SliceList_id INTEGER PRIMARY KEY, ImageList_id INTEGER, ', ...
                'Slice_Index INTEGER, X_Target_Position_um REAL, Y_Target_Position_um REAL, ', ...
                'Z_Target_Position_um REAL, Encoder_Timestamp_ms REAL)']);
            execute(conn, 'INSERT INTO MetadataList VALUES (1, 1001, ''{"a": 1}''), (2, 1002, ''{"a": 2}'')');
            % 8 images (2 cameras x 2 channels x 2 tiles, not in name order),
            % one without metadata, one whose hash has no metadata row and
            % one without slices
            imageRow = ['INSERT INTO ImageList VALUES (%d, ''acq_Cam%s_ch%d_%dnm_%s.tif'', ', ...
                '%d, %d, %d, 0, %g, %d, %d, 0, 0, 0, '''', %s)'];
            cams = 'AB';
            k = 0;
            for cam = [1, 0]
                for ch = [1, 0]
                    for x = [1, 0]
                        k = k + 1;
                        laser = 488 + 72 * ch;
                        execute(conn, sprintf(imageRow, k, cams(cam + 1), ch, laser, sprintf('%03dx', x), ...
                            cam, ch, laser, 0.5 * k, 5000 + 7 * k, x, num2str(1001 + mod(k, 2))));
                    end
                end
            end
            execute(conn, sprintf(imageRow, 9, 'A', 0, 488, 'nohash', 0, 0, 488, 0, 0, 5, 'NULL'));
            execute(conn, sprintf(imageRow, 10, 'A', 0, 488, 'nometa', 0, 0, 488, 0, 0, 6, '9999'));
            execute(conn, sprintf(imageRow, 11, 'A', 0, 488, 'noslice', 0, 0, 488, 0, 0, 7, '1001'));
            % three slices an image, the images' slices interleaved as the
            % cameras record them; one slice of image 3 without a target x
            for s = 0 : 2
                for k = [1 : 8, 9, 10]
                    xt = sprintf('%g', 100 * mod(k, 2) + 0.25 * s + 0.125);
                    if k == 3 && s == 1
                        xt = 'NULL';
                    end
                    execute(conn, sprintf('INSERT INTO SliceList VALUES (NULL, %d, %d, %s, %g, %g, %g)', ...
                        k, s, xt, -20.5 + 0.5 * k, 3 + 0.25 * s, 1000 * k + 10 * s));
                end
            end
            hashCol = '';
            outerHash = 'JSONsettings_MD5';
            if strcmp(hash, 'returned')
                hashCol = 'ImageList.JSONsettings_MD5 AS JSONsettings_MD5, ';
            elseif strcmp(hash, 'other')
                outerHash = 'JSONsettings_MD5_old';
            end
            execute(conn, ['CREATE VIEW EverythingView_NoNull AS SELECT ImageList.ImageList_id ImageList_id, ', ...
                hashCol, 'IFNULL(Filename, ''N/A'') AS Filename, IFNULL(Camera, ''N/A'') AS Camera, ', ...
                'IFNULL(Channel, ''N/A'') AS Channel, IFNULL(Laser_nm, ''N/A'') AS Laser_nm, ', ...
                'IFNULL(File_Stack_Index, ''N/A'') AS File_Stack_Index, IFNULL(Elapsed_ms, ''N/A'') AS Elapsed_ms, ', ...
                'IFNULL(Absolute_ms, ''N/A'') AS Absolute_ms, IFNULL(Tile_index_X, ''N/A'') AS Tile_index_X, ', ...
                'IFNULL(Tile_index_Y, ''N/A'') AS Tile_index_Y, IFNULL(Tile_index_Z, ''N/A'') AS Tile_index_Z, ', ...
                'IFNULL(Tile_index_T, ''N/A'') AS Tile_index_T, IFNULL(Frame_Type, ''N/A'') AS Frame_Type, ', ...
                'IFNULL(X_Target_Position_um, ''N/A'') AS X_Target_Position_um, ', ...
                'IFNULL(Y_Target_Position_um, ''N/A'') AS Y_Target_Position_um, ', ...
                'IFNULL(Z_Target_Position_um, ''N/A'') AS Z_Target_Position_um, ', ...
                'IFNULL(Encoder_Timestamp_ms, ''N/A'') AS Encoder_Timestamp_ms, ', ...
                'IFNULL(MetadataList.JSONsettings, ''N/A'') AS JSONsettings ', ...
                'FROM SliceList INNER JOIN ImageList ON ImageList.ImageList_id = SliceList.ImageList_id ', ...
                'INNER JOIN MetadataList ON MetadataList.JSONsettings_MD5 = ImageList.JSONsettings_MD5 ', ...
                'ORDER BY Filename, Encoder_Timestamp_ms']);
            execute(conn, ['CREATE VIEW ImageListView_NoNull AS SELECT ImageList_id, ', ...
                'IFNULL(Filename, ''N/A'') AS Filename, IFNULL(Camera, ''N/A'') AS Camera, ', ...
                'IFNULL(Channel, ''N/A'') AS Channel, IFNULL(Laser_nm, ''N/A'') AS Laser_nm, ', ...
                'IFNULL(File_Stack_Index, ''N/A'') AS File_Stack_Index, IFNULL(Elapsed_ms, ''N/A'') AS Elapsed_ms, ', ...
                'IFNULL(Absolute_ms, ''N/A'') AS Absolute_ms, IFNULL(Tile_index_X, ''N/A'') AS Tile_index_X, ', ...
                'IFNULL(Tile_index_Y, ''N/A'') AS Tile_index_Y, IFNULL(Tile_index_Z, ''N/A'') AS Tile_index_Z, ', ...
                'IFNULL(Tile_index_T, ''N/A'') AS Tile_index_T, IFNULL(Frame_Type, ''N/A'') AS Frame_Type, ', ...
                'IFNULL(', outerHash, ', ''N/A'') AS JSONsettings_MD5, JSONsettings, ', ...
                'avg(X_Target_Position_um) Mean_X_Target_Position_um, avg(Y_Target_Position_um) Mean_Y_Target_Position_um, ', ...
                'avg(Z_Target_Position_um) Mean_Z_Target_Position_um FROM EverythingView_NoNull ', ...
                'GROUP BY ImageList_id ORDER BY Filename, Encoder_Timestamp_ms']);
            clear cleanup;
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
