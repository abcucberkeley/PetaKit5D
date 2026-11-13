classdef StitchGeometryTest < pktest.TestCase
    % geometry helpers of the stitching: overlaps of tile cuboids
    % (check_cuboids_overlaps), the region of a minor tile that the major
    % tiles leave uncovered (check_major_tile_cover), and the axis order
    % check of XR_stitching_frame_zarr_dev_v1

    methods (Test, TestTags = {'Unit'})
        function cuboidsThatOnlyTouchDoNotOverlap(testCase)
            % cuboids that share only a face, an edge or a corner don't
            % overlap, also for a caller that asks for one output: the
            % stitching's recheck after the registration's shifts did, and
            % so kept such pairs, with half-overlap regions computed from a
            % zero overlap
            c = [0, 0, 0, 1, 1, 1];
            js = [1, 0, 0, 2, 1, 1; 0, 1, 0, 1, 2, 1; 0, 0, 1, 1, 1, 2; 1, 1, 0, 2, 2, 1; 1, 1, 1, 2, 2, 2; ...
                0.5, 0.5, 0.5, 2, 2, 2; 1.5, 0, 0, 2, 1, 1];
            expected = [false(5, 1); true; false];
            testCase.verifyEqual(check_cuboids_overlaps(c, js, false), expected, 'one output');
            [isOverlap, regions] = check_cuboids_overlaps(c, js, false);
            testCase.verifyEqual(isOverlap, expected, 'two outputs');
            testCase.verifyEqual(regions(6, :), [0.5, 0.5, 0.5, 1, 1, 1]);
            % in 2D a single z plane is no singularity
            c = [0, 0, 0, 1, 1, 0];
            js = [0.5, 0.5, 0, 2, 2, 0; 1, 0, 0, 2, 1, 0];
            testCase.verifyEqual(check_cuboids_overlaps(c, js, true), [true; false], '2D');
            testCase.verifyEqual(check_cuboids_overlaps(c, js, false), [false; false], '3D');
        end

        function majorTileCoverBoxHoldsEveryUncoveredPiece(testCase)
            % a major tile covers the middle of a minor tile: the box of the
            % minor tile's region that large-file feather blending reads
            % must hold the uncovered pieces on both sides (it ended at the
            % smallest end of the pieces, so the second was left out)
            [cover, box] = check_major_tile_cover([1, 1, 1, 10, 10, 10; 1, 4, 1, 10, 7, 10], [false; true]);
            testCase.verifyEqual(cover, [false; true]);
            testCase.verifyEqual(box(1, :), [1, 1, 1, 10, 10, 10]);
            % one uncovered piece, and none
            [cover, box] = check_major_tile_cover([1, 1, 1, 10, 10, 10; 1, 4, 1, 10, 10, 10], [false; true]);
            testCase.verifyEqual(cover, [false; true]);
            testCase.verifyEqual(box(1, :), [1, 1, 1, 10, 4, 10]);
            [cover, box] = check_major_tile_cover([2, 2, 2, 9, 9, 9; 1, 1, 1, 10, 10, 10], [false; true]);
            testCase.verifyEqual(cover, [true; true]);
            testCase.verifyEqual(box(1, :), zeros(1, 6));
        end

        function axisOrderMustBeAPermutationOfXYZ(testCase)
            % the check never raised (~regexpi(...) is false for no match
            % too), so the stitching sorted any order into some permutation;
            % 'z,x,y' was missing from it. It comes before the tiles are
            % read, so the tiles here don't exist: valid orders end in
            % another error later.
            dataPath = testCase.createTempFolder();
            for ao = {'x,x,z', 'a,b,c', 'x,y', 'x,y,z,y', 'x,y;z'}
                testCase.verifyError(@() testCase.frameWithAxisOrder(dataPath, ao{1}), ...
                    'XR_stitching_frame_zarr_dev_v1:axisOrder', ao{1});
            end
            for ao = {'x,y,z', 'y,x,z', 'z,y,x', 'x,z,y', 'z,x,y', 'y,z,x', '-x,y,z', 'x,-y,z', '-z,x,-y', 'yxz', ' y, x, z'}
                id = '';
                try
                    testCase.frameWithAxisOrder(dataPath, ao{1});
                catch ME
                    id = ME.identifier;
                end
                testCase.verifyNotEqual(id, 'XR_stitching_frame_zarr_dev_v1:axisOrder', ao{1});
            end
        end
    end

    methods (Access = private)
        function [] = frameWithAxisOrder(~, dataPath, axisOrder)
            fns = {[dataPath, 't1.tif'], [dataPath, 't2.tif']};
            evalc(['XR_stitching_frame_zarr_dev_v1(fns, [0, 0, 0; 1, 0, 0], ''axisOrder'', axisOrder, ', ...
                '''parseCluster'', false, ''ResultPath'', [dataPath, ''res/out.zarr''])']);
        end
    end
end
