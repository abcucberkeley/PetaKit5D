function [d_shift] = stitch_local_assignment(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, xcorr_thresh)
% use max xcorr as weights to construct a graph, and use MST to trim the
% graph, and use DFS to assign the absolute shift based on a tile's
% predecessor. 


d_shift = zeros(nF, 3);
% also remove the pair with very small max corr, i.e., <0.5 
max_xcorr_mat(max_xcorr_mat(:, 3) < xcorr_thresh, 3) = 0.001;        
G = graph(max_xcorr_mat(:, 1), max_xcorr_mat(:, 2), -max_xcorr_mat(:, 3));
T = minspantree(G, 'type', 'forest');
% absolute_shift_mat has a row [i, j, x, y, z] for each tile pair with an
% xcorr, and only the pairs in the tree are used below

% xruan (11/05/2021) change to use DFS search to just shift to the tile's precessor
% search from node 1
v = dfsearch(T, 1);
% indexed by tile, and tiles without overlaps are not in v
visit_flag = false(nF, 1);
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
    ind = find(absolute_shift_mat(:, 1) == s & absolute_shift_mat(:, 2) == t, 1);

    d_shift(n_i, :) = d_shift(n_pr, :) + st_sign .* absolute_shift_mat(ind, 3 : 5);
    visit_flag(n_i) = true;
end

end
