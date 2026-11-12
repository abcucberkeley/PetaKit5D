classdef StitchShiftAssignmentTest < pktest.TestCase
    % stitch_global_grid_assignment (shiftMethod 'grid', the default): the
    % tile shifts from the pairwise registration shifts, by constrained
    % least squares (stitch_assignment_lsqlin). Above 400 tiles it uses
    % sparse matrices. Also stitch_global_assignment ('global') and
    % stitch_global_group_assignment ('group').

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

    methods (Static, Access = private)
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
