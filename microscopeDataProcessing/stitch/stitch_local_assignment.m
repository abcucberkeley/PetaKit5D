function [d_shift] = stitch_local_assignment(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, xcorr_thresh)
% use max xcorr as weights to construct a graph, and use MST to trim the
% graph, and use DFS to assign the absolute shift based on a tile's
% predecessor. 


d_shift = zeros(nF, 3);
% also remove the pair with very small max corr, i.e., <0.5 
max_xcorr_mat(max_xcorr_mat(:, 3) < xcorr_thresh, 3) = 0.001;        
G = graph(max_xcorr_mat(:, 1), max_xcorr_mat(:, 2), -max_xcorr_mat(:, 3));
T = minspantree(G, 'type', 'forest');
aj = full(adjacency(T));
if any(size(aj) ~= nF)
    aj = padarray(aj, nF - size(aj), 0, 'post');
end
% overlap_matrix = overlap_matrix .* aj;
[inds_i, inds_j] = find(overlap_matrix .* aj);
absolute_shift_mat_orig = absolute_shift_mat;
absolute_shift_mat = absolute_shift_mat * 0;        
if numel(inds_i) > 0
    tinds = (inds_i - 1) * nF - inds_i .* (inds_i + 1) / 2 + inds_j;
    absolute_shift_mat(tinds, :) = absolute_shift_mat_orig(tinds, :);
end

% xruan (11/05/2021) change to use DFS search to just shift to the tile's precessor
% search from node 1
v = dfsearch(T, 1);
visit_flag = false(numel(v), 1);
for i = 1 : numel(v)
    n_i = v(i);
    if i == 1
        visit_flag(n_i) = true;
        continue;
    end

    n_nbs = neighbors(T, n_i);
    n_pr = n_nbs(visit_flag(n_nbs));
    if n_i > n_pr
        s = n_pr;
        t = n_i;
        st_sign = 1;
    else
        s = n_i;
        t = n_pr;
        st_sign = -1;
    end
    ind = (s - 1) * nF - s * (s + 1) / 2 + t;

    d_shift(n_i, :) = d_shift(n_pr, :) + st_sign .* absolute_shift_mat(ind, :);
    visit_flag(n_i) = true;
end

end
