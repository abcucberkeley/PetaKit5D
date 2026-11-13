classdef StitchShiftAssignmentTest < pktest.TestCase
    % stitch_global_grid_assignment (shiftMethod 'grid', the default): the
    % tile shifts from the pairwise registration shifts, by constrained
    % least squares (stitch_assignment_lsqlin). Above 400 tiles it uses
    % sparse matrices. Also stitch_global_assignment ('global'),
    % stitch_global_group_assignment ('group'), stitch_local_assignment
    % ('local') and stitch_shift_assignment on registration results written
    % to disk ('test', group files, the registration jobs' memory).

    methods (Test, TestTags = {'Unit'})
        function allPairsBelowTheThresholdKeepTheStagePositions(testCase)
            % lsqlin errored without rows in its objective
            args = testCase.gridProblem(4);
            args{2}(:, 3) = 0.1;
            nF = args{1};
            testCase.verifyEqual(stitch_global_grid_assignment(args{:}), zeros(nF, 3));
            testCase.verifyEqual(stitch_global_assignment(args{1 : 6}), zeros(nF, 3));
        end

        function onePairOfTiles(testCase)
            % one pair (one row in the objective and in each half of the
            % bounds): the second tile moves by the pair's shift (even, so
            % that the tiles at mean zero, at -shift / 2 and shift / 2, round
            % exactly)
            tileIdx = [0, 0, 0, 1; 1, 0, 0, 1];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            args = {2, [ti, tj, 0.9], [ti, tj, 8, -4, 2], overlap_matrix, [-50, -50, -10, 50, 50, 10], ...
                0.25, tileIdx, [1, 0.1, 10]};
            testCase.verifyEqual(stitch_global_grid_assignment(args{:}), [0, 0, 0; 8, -4, 2]);
            testCase.verifyEqual(stitch_global_assignment(args{1 : 6}), [0, 0, 0; 8, -4, 2]);
        end

        function freeTilesKeepTheirStagePositionsRelativeToTheirNeighbors(testCase)
            % tiles whose pairs are all below the threshold are free in the
            % objective: they were where lsqlin's iterations stopped. They
            % keep their stage positions relative to their neighbors: each
            % free tile (no two are neighbors) gets the mean of its
            % neighbors' shifts, also with asymmetric bounds (here -20 .. 8
            % voxels), the other tiles are at the optimum, all at mean zero.
            [args, trueOff, free] = testCase.freeTileProblem();
            ti = args{2}(:, 1);
            tj = args{2}(:, 2);
            x = trueOff;
            for k = free(:)'
                nb = [ti(tj == k); tj(ti == k)];
                x(k, :) = mean(x(nb, :), 1);
            end
            x = x - mean(x, 1);
            testCase.assertGreaterThan(min(abs(abs(x - floor(x)) - 0.5), [], 'all'), 1e-3);
            expected = round(x) - round(x(1, :));
            testCase.verifyEqual(stitch_global_grid_assignment(args{:}), expected);
            testCase.verifyEqual(stitch_global_assignment(args{1 : 6}), expected);
        end

        function freeTileHeldByItsBounds(testCase)
            % a free tile whose stage position relative to its neighbor is
            % outside its bounds goes to the nearest bound: three tiles in a
            % row, the last free, its pair bounded to 2 .. 8 voxels in x
            tileIdx = [(0 : 2)', zeros(3, 2), ones(3, 1)];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            mas = [-50, -50, -50, 50, 50, 50; 2, -50, -50, 8, 50, 50];
            args = {3, [ti, tj, [0.9; 0.1]], [ti, tj, zeros(2, 3)], overlap_matrix, mas, 0.25, tileIdx, [1, 1, 1]};
            % x = [-2, -2, 4] / 3: rounded [-1, -1, 1]
            testCase.verifyEqual(stitch_global_grid_assignment(args{:}), [0, 0, 0; 0, 0, 0; 2, 0, 0]);
        end

        function boundsOfAFreeTileHoldItsNeighbors(testCase)
            % a free tile that overlaps two tiles within 1 voxel each limits
            % how far apart they are: the least squares alone (10 voxels
            % apart) can't place it, so the bounded solve gives them 2 voxels
            % apart and the free tile in between
            tileIdx = [(0 : 2)', zeros(3, 2), ones(3, 1)];
            ti = [1; 1; 2];
            tj = [2; 3; 3];
            overlap_matrix = false(3);
            overlap_matrix(sub2ind([3, 3], ti, tj)) = true;
            rel = [10, 0, 0; 0, 0, 0; 0, 0, 0];
            mas = [-50, -50, -50, 4, 50, 50; -1, -50, -50, 1, 50, 50; -1, -50, -50, 1, 50, 50];
            args = {3, [ti, tj, [0.9; 0.1; 0.1]], [ti, tj, rel], overlap_matrix, mas, 0.25, tileIdx, [1, 1, 1]};
            d = stitch_global_grid_assignment(args{:});
            testCase.verifyEqual(d(:, 1), [0; 2; 1]);
        end

        function groupsWithoutPairsAboveTheThresholdKeepTheirStagePositions(testCase)
            % the group method (shiftMethod 'group') errored for a group
            % without a pair above the threshold to another group (it has no
            % nodes in step 2), or in lsqlin if no group had one. Six tiles in
            % a row in three groups of two: pairs (1, 2), (3, 4) and (5, 6)
            % within the groups, (2, 3) and (4, 5) between them.
            tileIdx = [(0 : 5)', zeros(6, 2), ones(6, 1)];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            sx = [2; 2; -2; 7; 4];
            mas = repmat([-50, -50, -50, 50, 50, 50], 5, 1);
            assign = @(xc) stitch_global_group_assignment(6, [ti, tj, xc], [ti, tj, sx, zeros(5, 2)], ...
                overlap_matrix, mas, 0.25, tileIdx, [1, 1, 1], [1; 1; 2; 2; 3; 3], [], '', false, 1, false, '');
            % groups 1 and 2 linked by (2, 3); group 3 free, with tile 5 at
            % tile 4's shift (its stage position relative to it)
            testCase.verifyEqual(assign([0.9; 0.9; 0.9; 0.1; 0.9]), [0, 2, 4, 2, 2, 6]' .* [1, 0, 0]);
            % no pair above the threshold between groups: each group at its
            % stage position relative to its neighbors (the first unmoved)
            testCase.verifyEqual(assign([0.9; 0.1; 0.9; 0.1; 0.9]), [0, 2, 2, 0, 0, 4]' .* [1, 0, 0]);
            % no pair above the threshold: the stage positions
            testCase.verifyEqual(assign(0.1 * ones(5, 1)), zeros(6, 3));
        end

        function groupLabelsOutOfTileOrderGiveTheSameShifts(testCase)
            % step 2 of the group method took only the pairs whose first tile
            % is in the group of the lower label, so with the labels swapped
            % the pair between the groups was dropped and the groups kept
            % their stage positions. Four tiles in a row, groups of two.
            tileIdx = [(0 : 3)', zeros(4, 2), ones(4, 1)];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            mas = repmat([-50, -50, -50, 50, 50, 50], 3, 1);
            assign = @(g) stitch_global_group_assignment(4, [ti, tj, 0.9 * ones(3, 1)], [ti, tj, [2; 2; -4], zeros(3, 2)], ...
                overlap_matrix, mas, 0.25, tileIdx, [1, 1, 1], g, [], '', false, 1, false, '');
            expected = [0; 2; 4; 0] .* [1, 0, 0];
            testCase.verifyEqual(assign([1; 1; 2; 2]), expected);
            testCase.verifyEqual(assign([2; 2; 1; 1]), expected);
            % and the pair's bounds turned around with it: tile 3 at most 1
            % voxel right of tile 2
            mas(2, 4) = 1;
            assign = @(g) stitch_global_group_assignment(4, [ti, tj, 0.9 * ones(3, 1)], [ti, tj, [2; 2; -4], zeros(3, 2)], ...
                overlap_matrix, mas, 0.25, tileIdx, [1, 1, 1], g, [], '', false, 1, false, '');
            expected = [0; 2; 3; -1] .* [1, 0, 0];
            testCase.verifyEqual(assign([1; 1; 2; 2]), expected, 'bounded');
            testCase.verifyEqual(assign([2; 2; 1; 1]), expected, 'bounded, labels swapped');
        end

        function groupsOfOneTileOnTheClusterPath(testCase)
            % the group method's cluster path (more than 100 tiles and more
            % than one group) errored for a group of one tile, which has no
            % task (an empty function string); here as the local path. 102
            % tiles in a row, the last in a group of its own.
            n = 102;
            tileIdx = [(0 : n - 1)', zeros(n, 2), ones(n, 1)];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            sx = randi([-3, 3], n - 1, 1);
            xcorrDir = testCase.createTempFolder();
            assign = @(parseCluster) stitch_global_group_assignment(n, [ti, tj, 0.9 * ones(n - 1, 1)], [ti, tj, sx, zeros(n - 1, 2)], ...
                overlap_matrix, repmat([-50, -50, -50, 50, 50, 50], n - 1, 1), 0.25, tileIdx, [1, 1, 1], [ones(n - 1, 1); 2], ...
                [], xcorrDir, parseCluster, 1, false, '');
            [~, d] = evalc('assign(true)');
            testCase.verifyEqual(d, [0; cumsum(sx)] .* [1, 0, 0]);
            testCase.verifyEqual(d, assign(false));
        end

        function localAssignmentPlacesEachTreeOfTheForest(testCase)
            % the local method placed only tile 1's tree (DFS from tile 1):
            % the tiles of the other trees stayed at their stage positions.
            % Each tree's first tile keeps its position, as tile 1 does.
            overlap_matrix = false(5);
            overlap_matrix(1, 2) = true;
            overlap_matrix(3, 4) = true;
            d = stitch_local_assignment(5, [1, 2, 0.9; 3, 4, 0.9], [1, 2, 5, 1, 0; 3, 4, -3, 2, 1], overlap_matrix, 0.25);
            testCase.verifyEqual(d, [0, 0, 0; 5, 1, 0; 0, 0, 0; -3, 2, 1; 0, 0, 0]);
            % a single tile (no pairs: a graph without nodes)
            testCase.verifyEqual(stitch_local_assignment(1, zeros(0, 3), zeros(0, 5), false, 0.25), zeros(1, 3));
        end

        function testMethodAssignsThePairShifts(testCase)
            % the 'test' method took the pairs' values at their index in
            % the list of all pairs (an index error for three tiles in a
            % row), the tile numbers as x and y shifts, and its bounds as
            % if every second tile followed the first; it gives the
            % registration's shifts, and the grid method's with its weights
            % (10, 1, 0.1 for y, x, z) where no bound holds them
            trueOff = [0, 0, 0; 6, -2, -1; 1, 1, -2];
            [args, p] = testCase.shiftAssignmentArgs('test', 3, 1, trueOff);
            [~, ~, dxyz] = evalc('stitch_shift_assignment(args{:})');
            testCase.verifyEqual(round(dxyz ./ p.vs, 6), trueOff, 'three tiles in a row');
            % tiles numbered against x, 6 voxels of overlap: the second tile
            % comes 8 voxels closer, within its bounds of -5 .. 12 voxels
            trueOff = [0, 0, 0; 8, 3, -2];
            [args, p] = testCase.shiftAssignmentArgs('test', 2, 1, trueOff, 'overlap', [6, 10], 'reverseX', true);
            [~, ~, dxyz] = evalc('stitch_shift_assignment(args{:})');
            testCase.verifyEqual(round(dxyz ./ p.vs, 6), trueOff, 'tiles numbered against x');
            % a tile whose pair is below the threshold keeps its stage
            % position relative to its neighbor (the solution made unique),
            % also where its bounds are -12 .. 2 voxels along x (3 voxels of
            % overlap), not symmetric about that position
            [args, p] = testCase.shiftAssignmentArgs('test', 3, 1, [0, 0, 0; 1, -2, -1; 5, 1, -2], 'overlap', [3, 10], ...
                'lowPairs', 2);
            [~, ~, dxyz] = evalc('stitch_shift_assignment(args{:})');
            testCase.verifyEqual(round(dxyz ./ p.vs, 6), [0, 0, 0; 1, -2, -1; 1, -2, -1], 'a tile with its pair below the threshold');
            % noisy shifts on a 3 x 2 grid
            args = testCase.shiftAssignmentArgs('test', 3, 2, [randi([-2, 2], 6, 2), randi([-1, 1], 6, 1)], 'noise', 0.4);
            [~, ~, dxyz] = evalc('stitch_shift_assignment(args{:})');
            args{13} = 'grid';
            args{15} = [10, 1, 0.1];
            [~, ~, expected] = evalc('stitch_shift_assignment(args{:})');
            testCase.verifyEqual(dxyz, expected, 'the grid method with its weights');
        end

        function registrationMemoryUsesTheTileSizeAlongXYZ(testCase)
            % the registration jobs' memory took the size of the overlap
            % plus the max offset along x, y, z, at most the tile's [ny nx
            % nz]: two tiles of 8192 x 1024 x 1000 (y, x, z) side by side
            % along x need 112 x 8192 x 1000 voxels (4 GB, single), 128 GB
            % with the factor of 32; they got 112 x 1024 x 1000 (1 GB, 32 GB)
            args = testCase.shiftAssignmentArgs('grid', 2, 1, zeros(2, 3), 'size', [8192, 1024, 1000], 'overlap', [100, 0]);
            logFn = testCase.recordFrameworkMemory();
            evalc('stitch_shift_assignment(args{:})');
            testCase.verifyEqual(str2double(readlines(logFn, 'EmptyLineRule', 'skip')), 128);
        end

        function groupFileWithoutDataFoldersGroupsEveryTile(testCase)
            % a group file of 4 columns (x, y, z, g) put its tiles in data
            % folder 0, so every tile was 1 away from its group's entry and
            % the method warned that tiles are not in the groups
            % (offsets whose differences are within the z bounds of 5 voxels)
            [args, p] = testCase.shiftAssignmentArgs('group', 4, 2, [randi([-3, 3], 8, 2), randi([-2, 2], 8, 1)], ...
                'groups', [1; 1; 2; 2; 1; 1; 2; 2]);
            testCase.verifyWarningFree(@() StitchShiftAssignmentTest.runQuietly(@() stitch_shift_assignment(args{:})));
            [~, ~, dxyz] = evalc('stitch_shift_assignment(args{:})');
            testCase.verifyEqual(round(dxyz ./ p.vs, 6), p.trueOff - p.trueOff(1, :));
            % without a group file the groups are the tiles' data folders
            % (one here)
            args{16} = '';
            [~, ~, dxyz] = evalc('stitch_shift_assignment(args{:})');
            testCase.verifyEqual(round(dxyz ./ p.vs, 6), p.trueOff - p.trueOff(1, :), 'without a group file');
        end

        function shiftsRoundTheExactOptimumAtAnyThreadCount(testCase)
            % the shifts are the exact optimum of the group (at mean zero)
            % rounded, at 1, 2, 4 and 24 threads: lsqlin's solution is
            % inexact and drifted along the group's constant by rounding
            % noise that depended on the thread count. Here the exact
            % optimum lies 1e-7 above or below half integers.
            [args, x] = testCase.halfIntegerProblem(10);
            expected = round(x) - round(x(1, :));
            nt = maxNumCompThreads;
            testCase.addTeardown(@maxNumCompThreads, nt);
            for n = [1, 2, 4, 24]
                maxNumCompThreads(n);
                testCase.verifyEqual(stitch_global_grid_assignment(args{:}), expected, sprintf('%d threads', n));
            end
        end

        function integerShiftsRoundAlike(testCase)
            % 2 x 2 tiles with integer pairwise shifts, one tile 2 voxels off
            % in x: at mean zero the exact shifts lie on half integers (0.5
            % and -1.5) and round alike, so that the tile stays 2 voxels off
            % relative to the others (they rounded away from zero, 1 and -2)
            tileIdx = [0, 0, 0, 1; 0, 1, 0, 1; 1, 0, 0, 1; 1, 1, 0, 1];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            pos = [0, 0, 0; 0, 0, 0; -2, 0, 0; 0, 0, 0];
            rel = pos(tj, :) - pos(ti, :);
            nP = numel(ti);
            args = {4, [ti, tj, 0.9 * ones(nP, 1)], [ti, tj, rel], overlap_matrix, ...
                repmat([-20, -20, -10, 20, 20, 10], nP, 1), 0.25, tileIdx, [1, 0.1, 10]};
            testCase.verifyEqual(stitch_global_grid_assignment(args{:}), pos);
            testCase.verifyEqual(stitch_global_assignment(args{1 : 6}), pos);
        end

        function activeBoundsGiveTheExactOptimum(testCase)
            % with a bound active, the shifts are the exact optimum: lsqlin's
            % solution refined on its active set (lsqlin alone stops about
            % 1e-8 away). Three tiles in a row, the first pair measured 10
            % voxels apart with a bound of u: the optimum along x is
            % [-2 u, u, u] / 3, here 1e-8 above a half integer.
            tileIdx = [(0 : 2)', zeros(3, 2), ones(3, 1)];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            u = 4.5 + 3e-8;
            mas = [-50, -50, -50, u, 50, 50; -50, -50, -50, 50, 50, 50];
            rel = [10, 0, 0; 0, 0, 0];
            args = {3, [ti, tj, [0.9; 0.9]], [ti, tj, rel], overlap_matrix, mas, 0.25, tileIdx, [1, 1, 1]};
            d = stitch_global_grid_assignment(args{:});
            x = [-2 * u; u; u] / 3;
            testCase.verifyEqual(d(:, 1), round(x) - round(x(1)));
        end

        function boundsHoldForTheirOwnPairs(testCase)
            % the bounds come in the order of max_xcorr_mat (the pairs by
            % first, then second tile); they bounded the pairs of
            % find(overlap_matrix)'s order instead, so a pair whose overlap
            % allows only 4 voxels got another pair's bounds
            gx = 4;
            gy = 3;
            [ix, iy] = meshgrid(0 : gx - 1, 0 : gy - 1);
            nF = gx * gy;
            tileIdx = [ix(:), iy(:), zeros(nF, 1), ones(nF, 1)];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            [fi, fj] = find(overlap_matrix);
            p = find(fi ~= ti | fj ~= tj, 1);
            testCase.assertNotEmpty(p);
            nP = numel(ti);
            rel = zeros(nP, 3);
            rel(p, 1) = 10;
            mas = repmat([-50, -50, -50, 50, 50, 50], nP, 1);
            mas(p, 4) = 4;
            args = {nF, [ti, tj, 0.9 * ones(nP, 1)], [ti, tj, rel], overlap_matrix, mas, 0.25, tileIdx, [1, 1, 1]};
            d = {stitch_global_grid_assignment(args{:}), stitch_global_assignment(args{1 : 6})};
            for k = 1 : 2
                pairShift = d{k}(tj, :) - d{k}(ti, :);
                testCase.verifyEqual(pairShift(p, 1), 4, sprintf('method %d', k));
                testCase.verifyTrue(all(pairShift >= mas(:, 1 : 3) & pairShift <= mas(:, 4 : 6), 'all'), sprintf('method %d', k));
            end
        end

        function separateGroupsOfALargeGridHaveMeanZeroEach(testCase)
            % above 400 tiles (sparse form) the solution of two groups of
            % overlapping tiles was moved to one mean, with each group's
            % constant arbitrary: each group is at mean zero, its shifts
            % (relative to its first tile) those of the group alone
            [args, groups] = testCase.twoGroupProblem(15, 14);
            d = stitch_global_grid_assignment(args{:});
            for g = 1 : 2
                [argsG, inds] = testCase.subProblem(args, groups == g);
                dG = stitch_global_grid_assignment(argsG{:});
                testCase.verifyEqual(d(inds, :) - d(inds(1), :), dG - dG(1, :), sprintf('group %d', g));
            end
        end

        function largeGridMatchesDenseSolver(testCase)
            % the sparse solution is moved to mean zero, as the dense solver
            % leaves it, so the rounded shifts agree except rarely by 1
            [args, trueOff] = testCase.gridProblem(21);
            d = stitch_global_grid_assignment(args{:});
            dense = testCase.denseAssignment(args{:});
            testCase.verifyLessThanOrEqual(max(abs(d - dense), [], 'all'), 1);
            testCase.verifyLessThanOrEqual(nnz(d ~= dense), 0.01 * numel(d), ...
                sprintf('%d of %d shifts differ from the dense solver.', nnz(d ~= dense), numel(d)));
            % and both find the simulated positions (noise 0.7 voxels)
            testCase.verifyLessThanOrEqual(max(abs(d - (trueOff - trueOff(1, :))), [], 'all'), 3);
        end
    end

    methods (Access = private)
        function [args, p] = shiftAssignmentArgs(testCase, method, gx, gy, trueOff, varargin)
            % stitch_shift_assignment's arguments for a gx x gy grid of tiles
            % ('y,x,z', [ny nx nz] = size, voxel sizes 0.108, 0.108 and 0.25
            % um) with overlap voxels between neighbors along x and y
            % (reverseX: the x positions decrease with the tile index) and
            % the registration results of the neighbor pairs on disk: the
            % shifts between the true offsets (voxels, x y z; plus noise
            % voxels of Gaussian noise), xcorr 0.9 (0.1 for the pairs lowPairs,
            % by first, then second tile); tile 1 as zarr metadata; groups: a
            % 4-column group file with these groups
            ip = inputParser;
            ip.addParameter('size', [64, 72, 20]);
            ip.addParameter('overlap', [12, 10]);
            ip.addParameter('reverseX', false);
            ip.addParameter('noise', 0);
            ip.addParameter('lowPairs', []);
            ip.addParameter('groups', []);
            ip.parse(varargin{:});
            pr = ip.Results;
            d = testCase.createTempFolder();
            xcorrDir = [d, 'xcorr/'];
            mkdir(xcorrDir);
            p.vs = [0.108, 0.108, 0.25];
            nF = gx * gy;
            [ix, iy] = ndgrid(0 : gx - 1, 0 : gy - 1);
            tileIdx = [ix(:), iy(:), zeros(nF, 1), ones(nF, 1)];
            ext = pr.size([2, 1, 3]);
            pos = tileIdx(:, 1 : 3);
            if pr.reverseX
                pos(:, 1) = gx - 1 - pos(:, 1);
            end
            xyz = pos .* (ext - 1 - [pr.overlap, 0]) .* p.vs;
            paths = arrayfun(@(k) sprintf('%stile_%d.zarr', d, k), (1 : nF)', 'unif', 0);
            createzarr(paths{1}, dataSize=pr.size, blockSize=[64, 64, 64], dtype='uint16');
            cuboids = [xyz, xyz + (ext - 1) .* p.vs];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            overlap_regions = zeros(nF * (nF - 1) / 2, 6);
            for k = 1 : numel(ti)
                [~, overlap_regions((ti(k) - 1) * nF - ti(k) * (ti(k) + 1) / 2 + tj(k), :)] = ...
                    check_cuboids_overlaps(cuboids(ti(k), :), cuboids(tj(k), :), false);
            end
            p.trueOff = trueOff;
            shifts = trueOff(tj, :) - trueOff(ti, :) + randn(numel(ti), 3) * pr.noise;
            order_flag = (0.5 - (xyz(ti, :) > xyz(tj, :))) * 2;
            xcorr = 0.9 * ones(numel(ti), 1);
            xcorr(pr.lowPairs) = 0.1;
            for f = unique(ti)'
                k = find(ti == f);
                relative_shift_mat = shifts(k, :) .* order_flag(k, :);
                max_xcorr_mat = xcorr(k);
                pair_indices = [ti(k), tj(k)];
                save(sprintf('%sxcorr_tile_%d.mat', xcorrDir, f), 'relative_shift_mat', 'max_xcorr_mat', 'pair_indices');
            end
            groupFile = '';
            if ~isempty(pr.groups)
                groupFile = [d, 'groups.csv'];
                writetable(table(tileIdx(:, 1), tileIdx(:, 2), tileIdx(:, 3), pr.groups(:), 'VariableNames', {'x', 'y', 'z', 'g'}), groupFile);
            end
            args = {paths, xcorrDir, repmat(pr.size, nF, 1), xyz, p.vs, [2, 1, 3], overlap_matrix, overlap_regions, ...
                [12, 12, 5], [1, 1, 1], 0.25, tileIdx, method, false, [1, 0.1, 10], groupFile, false, [], false, 1, false, '', 1};
        end

        function [logFn] = recordFrameworkMemory(testCase)
            % a stand-in for generic_computing_frameworks_wrapper that writes
            % each call's memAllocate to logFn and reports every task done
            folder = testCase.createTempFolder();
            logFn = [folder, 'memAllocate.txt'];
            fid = fopen([folder, 'generic_computing_frameworks_wrapper.m'], 'w');
            fprintf(fid, '%s\n', 'function [is_done_flag] = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, funcStrs, varargin)', ...
                'ip = inputParser;', 'ip.KeepUnmatched = true;', 'ip.addParameter(''memAllocate'', []);', 'ip.parse(varargin{:});', ...
                sprintf('writelines(string(ip.Results.memAllocate), ''%s'', ''WriteMode'', ''append'');', logFn), ...
                'is_done_flag = true(numel(outputFullpaths), 1);', 'end');
            fclose(fid);
            pktest.addFunctionFolder(folder, 'generic_computing_frameworks_wrapper');
            testCase.addTeardown(@() clear('generic_computing_frameworks_wrapper'));
            testCase.addTeardown(@() rmpath(folder));
        end
    end

    methods (Static, Access = private)
        function [varargout] = runQuietly(f)
            % the outputs of f(), without what it prints
            varargout = cell(1, nargout);
            [~, varargout{:}] = evalc('f()');
        end

        function [args, trueOff, free] = freeTileProblem()
            % a 5 x 5 grid of tiles with exact pairwise shifts, three tiles
            % (no two neighbors: an edge tile with one smaller and two larger
            % neighbors, one with two and one, an inner one) whose pairs are
            % below the xcorr threshold, bounds -20 .. 8 voxels
            g = 5;
            [ix, iy] = meshgrid(0 : g - 1, 0 : g - 1);
            nF = g * g;
            tileIdx = [ix(:), iy(:), zeros(nF, 1), ones(nF, 1)];
            trueOff = [randi([-3, 3], nF, 2), randi([-2, 2], nF, 1)];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            nP = numel(ti);
            free = [3; 15; 17];
            xc = 0.9 * ones(nP, 1);
            xc(ismember(ti, free) | ismember(tj, free)) = 0.1;
            rel = trueOff(tj, :) - trueOff(ti, :);
            mas = repmat([-20, -20, -20, 8, 8, 8], nP, 1);
            args = {nF, [ti, tj, xc], [ti, tj, rel], overlap_matrix, mas, 0.25, tileIdx, [1, 0.1, 10]};
        end

        function [args, x] = halfIntegerProblem(g)
            % a g x g grid (g even) with exact pairwise shifts of positions x
            % (mean zero) 1e-7 above or below half integers
            [ix, iy] = meshgrid(0 : g - 1, 0 : g - 1);
            nF = g * g;
            tileIdx = [ix(:), iy(:), zeros(nF, 1), ones(nF, 1)];
            v = randi([-10, 10], nF / 2, 3) + 0.5 + (2 * (rand(nF / 2, 3) > 0.5) - 1) * 1e-7;
            x = [v; -v];
            x = x(randperm(nF), :);
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            nP = numel(ti);
            rel = x(tj, :) - x(ti, :);
            mas = repmat([-100, -100, -100, 100, 100, 100], nP, 1);
            args = {nF, [ti, tj, 0.9 * ones(nP, 1)], [ti, tj, rel], overlap_matrix, mas, 0.25, tileIdx, [1, 0.1, 10]};
        end

        function [args, groups] = twoGroupProblem(gx, gy)
            % two gx x gy grids of tiles that don't overlap each other, with
            % noisy pairwise shifts, some below the xcorr threshold
            [ix, iy] = meshgrid(0 : gx - 1, 0 : gy - 1);
            n1 = gx * gy;
            tileIdx = [ix(:), iy(:), zeros(n1, 1), ones(n1, 1)];
            tileIdx = [tileIdx; tileIdx + [gx + 2, 0, 0, 0]];
            groups = [ones(n1, 1); 2 * ones(n1, 1)];
            nF = 2 * n1;
            trueOff = [randi([-15, 15], nF, 2), randi([-5, 5], nF, 1)];
            [ti, tj, overlap_matrix] = StitchShiftAssignmentTest.neighborPairs(tileIdx);
            nP = numel(ti);
            rel = trueOff(tj, :) - trueOff(ti, :) + randn(nP, 3) * 0.7;
            xc = 0.3 + 0.7 * rand(nP, 1);
            xc(rand(nP, 1) < 0.05) = 0.1;
            mas = repmat([-100, -100, -30, 97, 97, 27], nP, 1);
            args = {nF, [ti, tj, xc], [ti, tj, rel], overlap_matrix, mas, 0.25, tileIdx, [1, 0.1, 10]};
        end

        function [args, inds] = subProblem(args, mask)
            % the problem of the tiles in mask alone (renumbered)
            [nF, mx, ash, om, mas, thr, tileIdx, aw] = args{:};
            inds = find(mask);
            newIdx = zeros(nF, 1);
            newIdx(inds) = 1 : numel(inds);
            keep = mask(mx(:, 1)) & mask(mx(:, 2));
            mx = mx(keep, :);
            mx(:, 1 : 2) = newIdx(mx(:, 1 : 2));
            ash = ash(keep, :);
            ash(:, 1 : 2) = newIdx(ash(:, 1 : 2));
            args = {numel(inds), mx, ash, om(inds, inds), mas(keep, :), thr, tileIdx(inds, :), aw};
        end

        function [ti, tj, overlap_matrix] = neighborPairs(tileIdx)
            % the pairs of neighboring tiles (i < j, sorted by i, then j)
            nF = size(tileIdx, 1);
            overlap_matrix = false(nF);
            for i = 1 : nF
                for j = i + 1 : nF
                    overlap_matrix(i, j) = sum(abs(tileIdx(i, 1 : 3) - tileIdx(j, 1 : 3))) == 1;
                end
            end
            [ti, tj] = find(overlap_matrix);
            [~, sinds] = sortrows([ti, tj]);
            ti = ti(sinds);
            tj = tj(sinds);
        end

        function [args, trueOff] = gridProblem(g)
            % a g x g grid of tiles with random offsets and noisy pairwise
            % shifts between neighbors, some below the xcorr threshold
            rs = RandStream('twister', 'Seed', 1);
            [ix, iy] = meshgrid(0 : g - 1, 0 : g - 1);
            nF = g * g;
            tileIdx = [ix(:), iy(:), zeros(nF, 1), ones(nF, 1)];
            xyz = [tileIdx(:, 1) * 900, tileIdx(:, 2) * 900, zeros(nF, 1)];
            trueOff = [randi(rs, [-15, 15], nF, 2), randi(rs, [-5, 5], nF, 1)];
            overlap_matrix = false(nF);
            for i = 1 : nF
                for j = i + 1 : nF
                    overlap_matrix(i, j) = sum(abs(tileIdx(i, 1 : 3) - tileIdx(j, 1 : 3))) == 1;
                end
            end
            [ti, tj] = find(overlap_matrix);
            [~, sinds] = sortrows([ti, tj]);
            ti = ti(sinds);
            tj = tj(sinds);
            nP = numel(ti);
            rel = trueOff(tj, :) - trueOff(ti, :) + randn(rs, nP, 3) * 0.7;
            xc = 0.3 + 0.7 * rand(rs, nP, 1);
            xc(rand(rs, nP, 1) < 0.05) = 0.1;
            max_xcorr_mat = [ti, tj, xc];
            absolute_shift_mat = [ti, tj, rel];
            max_shift_l = -ones(nP, 1) .* [100, 100, 30];
            max_shift_u = repmat([97, 97, 27], nP, 1);
            order_mat = xyz(tj, :) >= xyz(ti, :);
            max_allow_shift = [max_shift_l .* order_mat + max_shift_u .* (order_mat - 1), ...
                max_shift_u .* order_mat + max_shift_l .* (order_mat - 1)];
            args = {nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, max_allow_shift, 0.25, tileIdx, [1, 0.1, 10]};
        end

        function [d_shift] = denseAssignment(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, max_allow_shift, xcorr_thresh, tileIdx, axisWeight)
            % the dense solve stitch_global_grid_assignment does up to 400 tiles
            [n_i, n_j] = find(overlap_matrix);
            neq = numel(n_i);
            R = zeros(neq, nF);
            R(sub2ind(size(R), 1 : neq, n_i')) = -1;
            R(sub2ind(size(R), 1 : neq, n_j')) = 1;
            filt_inds = max_xcorr_mat(:, 3) < xcorr_thresh;
            mx = max_xcorr_mat(~filt_inds, :);
            nP = size(mx, 1);
            % weights of neighbors along tile index x, y, z
            axisWeight_xyz = axisWeight([2, 1, 3]);
            w = zeros(nP, 1);
            for i = 1 : nP
                dIdx = abs(tileIdx(mx(i, 1), 1 : 3) - tileIdx(mx(i, 2), 1 : 3));
                if sum(dIdx) == 1
                    w(i) = axisWeight_xyz(dIdx == 1);
                end
            end
            R_w = zeros(nP, nF);
            R_w(sub2ind(size(R_w), 1 : nP, mx(:, 1)')) = -1;
            R_w(sub2ind(size(R_w), 1 : nP, mx(:, 2)')) = 1;
            W = diag(w);
            R_w = W .^ 0.5 * R_w;
            d_w = W .^ 0.5 * absolute_shift_mat(~filt_inds, 3 : 5);
            d_shift = zeros(nF, 3);
            for i = 1 : 3
                b = [max_allow_shift(:, 3 + i); -max_allow_shift(:, i)];
                evalc('d_shift(:, i) = lsqlin(R_w, d_w(:, i), [R; -R], b);');
            end
            d_shift = round(d_shift);
            d_shift = d_shift - d_shift(1, :);
        end
    end
end
