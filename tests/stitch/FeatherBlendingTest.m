classdef FeatherBlendingTest < pktest.TestCase
    % feather_blending_3d, which blends the tiles of a batch of
    % processStitchBlock by their distance maps: the whole batch, or, for a
    % batch that isn't full, the common overlap region of the tiles, written
    % at bbox into the block that holds the tiles. The mex takes uint16 and
    % single tiles with single weights; the MATLAB code blends the others.

    properties (TestParameter)
        % data orders other than 'y,x,z' (z last, z first, a 3-cycle)
        dataOrder = {'x,y,z', 'z,y,x', 'z,x,y'}
        % classes of the tiles (and the block) and of the weights that the
        % MATLAB code blends
        classes = struct('uint8Tiles', {{'uint8', 'single'}}, 'doubleTiles', {{'double', 'single'}}, ...
            'uint16TilesDoubleWeights', {{'uint16', 'double'}}, 'singleTilesDoubleWeights', {{'single', 'double'}})
    end

    methods (Test, TestTags = {'Unit'})
        function overlapIsBlendedIntoTheBlock(testCase, classes)
            % the MATLAB code wrote its blend (single) into a copy of a block
            % of another class (indexing4d_mex rejects the mix), which it
            % dropped: the block came back with the tiles' values, so uint8
            % tiles in a batch that isn't full weren't blended
            [tileClass, weightClass] = classes{:};
            if strcmp(weightClass, 'double')
                % the mex rejects the weights as an input error
                testCase.expectMexFallback('feather_blending_3d_with_indexing_mex');
            end
            [f, d, block, bbox] = FeatherBlendingTest.randomBatch(tileClass, weightClass);
            expected = block;
            expected(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)) = ...
                cast(FeatherBlendingTest.blend(f, d), tileClass);

            [out, mex_compute] = feather_blending_3d(f, d, block, bbox);

            testCase.verifyClass(out, tileClass);
            testCase.verifyEqual(out, expected, 'The overlap is not blended into the block.');
            testCase.verifyFalse(mex_compute);
        end

        function overlapIsBlendedIntoTheBlockWithoutIndexingMex(testCase)
            % without indexing4d_mex (a stub that errors) indexing4d writes
            % into a copy of the block, which must be kept
            testCase.applyFixture(pktest.FunctionStubFixture('indexing4d_mex'));
            testCase.expectMexFallback('indexing4d_mex');
            [f, d, block, bbox] = FeatherBlendingTest.randomBatch('uint8', 'single');
            expected = block;
            expected(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)) = ...
                cast(FeatherBlendingTest.blend(f, d), 'uint8');

            out = feather_blending_3d(f, d, block, bbox);

            testCase.verifyEqual(out, expected, 'The overlap is not blended into the block.');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function perTileMapsFollowTheDataOrder(testCase, dataOrder)
            % compute_tile_bwdist, a tile's map without the overlap regions
            % (rows x, y, z): the regions were read as [y x z] of the tile
            % (turned to put z last), whose axes are x, y, z then (and for
            % 'x,y,z' data), so x and y were exchanged in the maps of data
            % orders other than 'y,x,z'. The maps of the same tiles in each
            % data order must be those of 'y,x,z', in that order.
            [tiles, regions] = FeatherBlendingTest.twoTiles();
            [~, ref] = FeatherBlendingTest.tileMap(testCase.createTempFolder(), tiles, regions, 'y,x,z', 1);
            [perm, im] = FeatherBlendingTest.tileMap(testCase.createTempFolder(), tiles, regions, dataOrder, 1);
            testCase.verifyEqual(im, permute(ref, perm), sprintf('The map of %s data.', dataOrder));
        end

        function largeFileMapsWithZNotLast(testCase)
            % compute_tile_bwdist_mip_slabs of data with z not last: the map
            % (turned back to the data's order) was written with the turned
            % size, so it failed unless the slab is a cube
            [tiles, regions] = FeatherBlendingTest.twoTiles();
            ref = FeatherBlendingTest.slabMap(testCase.createTempFolder(), tiles{1}, regions, 'y,x,z');
            for order = {'z,y,x', 'z,x,y'}
                [im, perm] = FeatherBlendingTest.slabMap(testCase.createTempFolder(), tiles{1}, regions, order{1});
                testCase.verifyEqual(im, permute(ref, perm), sprintf('The map of %s data.', order{1}));
            end
        end

        function singleMapWithoutLocationIds(testCase)
            % compute_tile_distance_transform with singleDistMap (its
            % default) and no locIds (its default) failed (index of an empty
            % cell): the tiles are of one location, as with locIds of ones
            [tiles, regions] = FeatherBlendingTest.twoTiles();
            [bkFn, tileFns] = FeatherBlendingTest.writeBlockInfo(testCase.createTempFolder(), tiles, regions, [1, 2, 3]);
            out = cell(1, 2);
            args = {{}, {'locIds', [1; 1]}};
            for i = 1 : 2
                stitchPath = testCase.createTempFolder();
                [fns, idx] = compute_tile_distance_transform(bkFn, stitchPath, tileFns, 'singleDistMap', true, ...
                    'blockSize', [16, 16, 16], 'parseCluster', false, args{i}{:});
                testCase.assertNumElements(fns, 1);
                testCase.verifyEqual(idx, 1);
                out{i} = readzarr(fns{1});
            end
            testCase.verifyEqual(out{1}, out{2});
        end

        function oneBoxForThePerTileMaps(testCase)
            % compute_tile_distance_transform with one distance box for
            % per-tile maps (as the frame passes a 1 x 6 distBboxes when it
            % turns singleDistMap off for tiles of other sizes): index error
            % for the second tile. The box is the box of each tile's map.
            [tiles, regions] = FeatherBlendingTest.twoTiles();
            [bkFn, tileFns] = FeatherBlendingTest.writeBlockInfo(testCase.createTempFolder(), tiles, regions, [1, 2, 3]);
            box = [3, 4, 1, 18, 27, 6];
            out = cell(2, 2);
            boxes = {box, [box; box]};
            for i = 1 : 2
                stitchPath = testCase.createTempFolder();
                fns = compute_tile_distance_transform(bkFn, stitchPath, tileFns, 'singleDistMap', false, ...
                    'distBboxes', boxes{i}, 'blockSize', [16, 16, 16], 'parseCluster', false);
                testCase.assertNumElements(fns, 2);
                out(i, :) = cellfun(@readzarr, fns, 'unif', 0);
            end
            testCase.verifyEqual(out(1, :), out(2, :));
        end
    end

    methods (Static, Access = private)
        function [tiles, regions] = twoTiles()
            % two [y x z] tiles of [20 30 6] with holes, overlapping in x
            % (columns 21-30 of tile 1 are columns 1-10 of tile 2), and the
            % overlap regions in each tile (rows x, y, z)
            tiles = cell(1, 2);
            for i = 1 : 2
                im = uint16(randi([1, 100], [20, 30, 6]));
                im(5 : 8, 3 + 20 * (i - 1) : 6 + 20 * (i - 1), :) = 0;
                im(12 : 13, 15 : 16, 2 : 4) = 0;
                tiles{i} = im;
            end
            regions = {[21, 30; 1, 20; 1, 6], [1, 10; 1, 20; 1, 6]};
        end

        function [dataOrderMat, perm] = dataOrderOf(dataOrder)
            % data_order_mat as the frame makes it, and the permutation of a
            % [y x z] array to the data's order
            [~, dataOrderMat] = sort(strrep(dataOrder, ',', ''));
            [~, perm] = sort(dataOrderMat([2, 1, 3]));
        end

        function [bkFn, tileFns] = writeBlockInfo(dataPath, tiles, regions, perm)
            % the tiles (zarr, permuted) and the block information file
            % with their overlap regions
            tileFns = cell(1, numel(tiles));
            for i = 1 : numel(tiles)
                tileFns{i} = sprintf('%stile_%d.zarr', dataPath, i);
                writezarr(permute(tiles{i}, perm), tileFns{i}, blockSize=[16, 16, 16]);
            end
            overlap_matrix = logical([0, 1; 0, 0]);
            ol_region_cell = {regions};
            overlap_map_mat = [1, 3];
            bkFn = [dataPath, 'block_info.mat'];
            save(bkFn, 'tileFns', 'overlap_matrix', 'ol_region_cell', 'overlap_map_mat');
        end

        function [perm, im] = tileMap(dataPath, tiles, regions, dataOrder, tileInd)
            % compute_tile_bwdist's map of a tile in the data order
            [dataOrderMat, perm] = FeatherBlendingTest.dataOrderOf(dataOrder);
            bkFn = FeatherBlendingTest.writeBlockInfo(dataPath, tiles, regions, perm);
            outFn = sprintf('%smap_%d.zarr', dataPath, tileInd);
            compute_tile_bwdist(bkFn, tileInd, outFn, 10, false, [16, 16, 16], [], 'zstd', [], dataOrderMat, false);
            im = readzarr(outFn);
        end

        function [im, perm] = slabMap(dataPath, tile, regions, dataOrder)
            % compute_tile_bwdist_mip_slabs's map of a tile whose z MIP slabs
            % (pooling [2 2 2]) are those of a [y x z] tile in the data order
            [dataOrderMat, perm] = FeatherBlendingTest.dataOrderOf(dataOrder);
            bkFn = FeatherBlendingTest.writeBlockInfo(dataPath, {tile, tile}, regions, perm);
            slab = tile(:, :, 1 : 2 : end);
            mkdir([dataPath, 'MIP_Slabs_pooling_2_2_2']);
            writezarr(permute(slab, perm), [dataPath, 'MIP_Slabs_pooling_2_2_2/tile_1_z.zarr'], blockSize=[16, 16, 16]);
            outFn = [dataPath, 'map.zarr'];
            compute_tile_bwdist_mip_slabs(bkFn, 1, outFn, 10, true, [16, 16, 16], [], 'zstd', [2, 2, 2], [], ...
                dataOrderMat, false);
            im = readzarr(outFn);
        end

        function [f, d, block, bbox] = randomBatch(tileClass, weightClass)
            % two tiles of the overlap region with values to 200 (tile 2
            % without data in places, tile 1 with data everywhere), their
            % weights, and a block whose region bbox ([y0, x0, z0, y1, x1,
            % z1]) is the overlap
            sz = [6, 7, 5];
            f = cast(randi([1, 200], [sz, 2]), tileClass);
            f(randperm(prod(sz), 40) + prod(sz)) = 0;
            d = cast(rand([sz, 2]) * 0.9 + 0.1, weightClass);
            block = cast(randi([1, 200], [10, 12, 8]), tileClass);
            bbox = [3, 4, 2, 2 + sz(1), 3 + sz(2), 1 + sz(3)];
        end

        function [im] = blend(f, d)
            % weighted mean of the tiles with data
            w = d .* (f ~= 0);
            im = sum(single(f) .* w, 4) ./ sum(w, 4);
        end
    end
end
