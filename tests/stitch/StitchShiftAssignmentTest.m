classdef StitchShiftAssignmentTest < pktest.TestCase
    % stitch_global_grid_assignment (shiftMethod 'grid', the default): the
    % tile shifts from the pairwise registration shifts, by constrained
    % least squares. Above 400 tiles it uses sparse matrices.

    methods (Test, TestTags = {'Unit'})
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
