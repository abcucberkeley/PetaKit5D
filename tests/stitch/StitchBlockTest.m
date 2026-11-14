classdef StitchBlockTest < pktest.TestCase
    % processStitchBlock, which blends the tiles of each stitching batch
    % (block), on batches made here: feather blending of large files (distance
    % maps pooled in z, poolSize [1 1 2]), where the minor tiles (weights far
    % below the major tile's) are read only where the major tile has no data
    % or doesn't reach; the blend method's name; and blocks with a border
    % (BorderSize) through XR_stitching_frame_zarr_dev_v1.

    properties (Constant)
        XYPixelSize = 0.108
    end

    properties (TestParameter)
        borderBlend = {'feather', 'none'}
    end

    methods (Test, TestTags = {'Unit'})
        function largeFileHoleOfTheMajorTileGetsTheMinorTiles(testCase)
            % the box of the voxels without data in the major tile is in the
            % coordinates of the tiles' overlap region, but was used as the
            % batch's: the minor tile was read at another place (where the
            % major tile has data) and the hole stayed 0. Tiles cut from a
            % volume; the major tile has a hole that a minor tile fills.
            gt = uint16(randi([1000, 60000], [20, 24, 8]));
            boxes = [1, 1, 1, 20, 24, 8; 4, 4, 1, 11, 11, 8; 12, 14, 1, 19, 21, 8];
            tiles = StitchBlockTest.cut(gt, boxes);
            tiles{1}(6 : 9, 6 : 9, :) = 0;

            out = testCase.blendBatch(tiles, boxes, [9, 2, 2], [1, 1, 2], 'feather');

            testCase.verifyEqual(out, gt, 'The hole of the major tile was not filled from the minor tile.');
        end

        function largeFileMinorTileOutsideTheMajorTileIsRead(testCase)
            % the minor tiles the major tile doesn't cover were never read
            % (~any(major_cover(~minor_inds)) is always false): their parts
            % outside it stayed 0. Where they are inside, the major tile
            % wins; here the minor tile's values are 1000 above.
            gt = uint16(randi([1000, 50000], [19, 21, 8]));
            boxes = [1, 1, 1, 16, 20, 8; 4, 4, 1, 11, 11, 8; 12, 14, 1, 19, 21, 8];
            tiles = StitchBlockTest.cut(gt, boxes);
            tiles{3} = tiles{3} + 1000;
            expected = zeros(size(gt), 'uint16');
            expected(12 : 19, 14 : 21, :) = tiles{3};
            expected(1 : 16, 1 : 20, :) = gt(1 : 16, 1 : 20, :);

            out = testCase.blendBatch(tiles, boxes, [9, 2, 2], [1, 1, 2], 'feather');

            testCase.verifyEqual(out, expected, 'The minor tile outside the major tile was not blended in.');
        end

        function featherBlendMethodInAnyCase(testCase)
            % the feather steps compared the method with strcmpi and the
            % merge with strcmp, so 'Feather' computed the weights, then
            % failed in tile_blending_3d ('Unsupported blending method')
            gt = uint16(randi([1000, 60000], [16, 24, 8]));
            boxes = [1, 1, 1, 16, 14, 8; 1, 9, 1, 16, 24, 8];
            tiles = StitchBlockTest.cut(gt, boxes);
            tiles{2} = tiles{2} + 500;

            ref = testCase.blendBatch(tiles, boxes, [5, 3], [], 'feather');
            out = testCase.blendBatch(tiles, boxes, [5, 3], [], 'Feather');

            testCase.verifyEqual(out, ref);
        end
    end

    methods (Test, TestTags = {'Integration'})
        function blocksWithABorderGiveTheSameVolume(testCase, borderBlend)
            % with BorderSize, the blocks of a tile reached past its ends
            % (where its end or start lies in a neighbor's border) and
            % readzarr refused, and the last (smaller) batch was cropped to
            % the batch size; a border changes no voxel of the stitched
            % volume. Two by two tiles, 4 x 4 x 2 batches.
            dataPath = testCase.createTempFolder();
            gt = uint16(randi([1000, 60000], [100, 120, 32]));
            starts = [0, 0, 0; 36, 0, 0; 0, 48, 0; 36, 48, 0];
            tiles = struct('im', cell(4, 1), 'idx', [], 'xyz', []);
            for i = 1 : 4
                s = starts(i, :);
                tiles(i).im = gt(s(1) + 1 : s(1) + 64, s(2) + 1 : s(2) + 72, s(3) + 1 : s(3) + 32);
                tiles(i).idx = [s(2) > 0, s(1) > 0, 0];
                tiles(i).xyz = s([2, 1, 3]) * testCase.XYPixelSize;
            end
            [fsns, csvFn] = pktest.writeStitchTiles(dataPath, tiles);
            tab = readtable(csvFn, 'Delimiter', 'comma');
            tileFullpaths = cellfun(@(x) [dataPath, x, '.tif'], fsns, 'unif', 0);
            xyz = [tab.StageX_um_, tab.StageY_um_, tab.StageZ_um_];
            tileIdx = [cat(1, tiles.idx), ones(4, 1)];

            outs = cell(1, 2);
            borders = {[0, 0, 0], [5, 6, 3]};
            for i = 1 : 2
                resultFn = sprintf('%sstitch_%d/result.zarr', dataPath, i);
                mkdir(fileparts(resultFn));
                XR_stitching_frame_zarr_dev_v1(tileFullpaths, xyz, 'axisOrder', 'x,y,z', ...
                    'xyPixelSize', testCase.XYPixelSize, 'DSR', true, 'resampleType', 'isotropic', ...
                    'xcorrShift', false, 'tileIdx', tileIdx, 'BlendMethod', borderBlend, 'ResultPath', resultFn, ...
                    'blockSize', [32, 32, 16], 'batchSize', [32, 32, 16], 'BorderSize', borders{i}, ...
                    'parseCluster', false);
                outs{i} = readzarr(resultFn);
            end
            testCase.verifyEqual(outs{2}, outs{1}, 'The volume with block borders differs.');
            testCase.verifyEqual(outs{1}(1 : 100, 1 : 120, 1 : 32), gt);
        end
    end

    methods (Access = private)
        function [out] = blendBatch(testCase, tiles, boxes, mapValues, poolSize, blendMethod)
            % processStitchBlock of one batch, the whole volume, with the
            % tiles at boxes ([y0 x0 z0 y1 x1 z1] in the volume) and distance
            % maps of constant values (pooled by poolSize, or full size)
            dataPath = testCase.createTempFolder();
            nT = numel(tiles);
            nvSize = max(boxes(:, 4 : 6), [], 1);
            tileFns = arrayfun(@(j) sprintf('%stile_%d.zarr', dataPath, j), 1 : nT, 'unif', 0);
            mapFns = arrayfun(@(j) sprintf('%smap_%d.zarr', dataPath, j), 1 : nT, 'unif', 0);
            stchBlockInfo = struct('tileInd', num2cell(1 : nT), 'bCoords', [], 'bboxCoords', []);
            for j = 1 : nT
                writezarr(tiles{j}, tileFns{j}, blockSize=[16, 16, 16]);
                sz = size(tiles{j}, 1 : 3);
                if ~isempty(poolSize)
                    sz(3) = ceil(sz(3) / poolSize(3));
                end
                writezarr(mapValues(j) * ones(sz, 'single'), mapFns{j}, blockSize=[16, 16, 16]);
                stchBlockInfo(j).bCoords = boxes(j, :);
                stchBlockInfo(j).bboxCoords = [1, 1, 1, boxes(j, 4 : 6) - boxes(j, 1 : 3) + 1];
            end
            bkFn = [dataPath, 'block_info.mat'];
            save(bkFn, 'tileFns');
            stitchFn = [dataPath, 'stitched.zarr'];
            createzarr(stitchFn, dataSize=nvSize, blockSize=nvSize, dtype='uint16');
            mkdir([dataPath, 'flags']);

            processStitchBlock(1, bkFn, '', [dataPath, 'flags/blocks_1.mat'], stitchFn, stchBlockInfo, [], ...
                'imSize', nvSize, 'batchSize', nvSize, 'dtype', 'uint16', 'BlendMethod', blendMethod, ...
                'imdistFullpaths', mapFns, 'imdistFileIdx', (1 : nT)', 'poolSize', poolSize, 'weightDegree', 10);
            out = readzarr(stitchFn);
        end
    end

    methods (Static, Access = private)
        function [tiles] = cut(gt, boxes)
            % the tiles of gt at boxes ([y0 x0 z0 y1 x1 z1])
            tiles = cell(1, size(boxes, 1));
            for j = 1 : size(boxes, 1)
                b = boxes(j, :);
                tiles{j} = gt(b(1) : b(4), b(2) : b(5), b(3) : b(6));
            end
        end
    end
end
