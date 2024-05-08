function [d_shift] = stitch_global_group_assignment_test_1(nF, max_xcorr_mat, absolute_shift_mat, ...
    overlap_matrix, max_allow_shift, xcorr_thresh, tileIdx, axisWeight, grpIdx, ...
    cuboid_overlap_ij_mat, xcorrDir, parseCluster, nodeFactor, mccMode, configFile)
% Two-step optimization: 
% step 1: grid method for tiles within each group
% step 2: weighted constrained linear least square for shifts across groups
% with the results from the first step as constrains. 
% 
% the weight is the function of max shift (currently just y=x). 
% max allowed shift is based on the maxShift parameter and the number of
% overlap between tiles.
%
%
% xruan (04/13/2023): add support for cluster computing for step 1;
% simpilify and reduce variables for the optimization problem in step 2 (only
% include tiles with edges across tiles). 


% step 1: grid method for tiles within the same group

uniq_grp = unique(grpIdx);
ng = numel(uniq_grp);
parseCluster = parseCluster && nF > 100 && ng > 1;  

if parseCluster
    uuid = get_uuid();    
    xcorrOptDir = sprintf('%s/xcorr_group_optimization/', xcorrDir);
    if ~exist(xcorrOptDir, 'dir')
        mkdir(xcorrOptDir);
    end
end

d_shift_pre = zeros(nF, 3);
all_gp_inds = false(size(max_xcorr_mat, 1), 1);
gind_cell = cell(ng, 1);
inputFullpaths = cell(ng, 1);
outputFullpaths = cell(ng, 1);
funcStrs = cell(ng, 1);
for g = 1 : ng
    gind = find(grpIdx == uniq_grp(g));
    gind_cell{g} = gind;
    nF_g = numel(gind);
    % new file indices within the group
    new_finds = cumsum(grpIdx == uniq_grp(g));
    
    [~, gp_inds] = knnsearch(gind, [max_xcorr_mat(:, 1); max_xcorr_mat(:, 2)]);
    gp_inds = all(reshape(gp_inds , [], 2) == 0, 2);
    all_gp_inds = all_gp_inds | gp_inds;

    max_xcorr_mat_g = max_xcorr_mat(gp_inds, :);
    new_pinds = new_finds(max_xcorr_mat_g(:, 1 : 2)); 
    max_xcorr_mat_g(:, 1 : 2) = new_pinds;
    absolute_shift_mat_g = absolute_shift_mat(gp_inds, :);
    absolute_shift_mat_g(:, 1 : 2) = new_pinds;
    overlap_matrix_g = overlap_matrix(gind, gind);
    max_allow_shift_g = max_allow_shift(gp_inds, :);
    tileIdx_g = tileIdx(gind, :);
    
    if parseCluster
        inputFullpath = sprintf('%s/xcorr_group_%d_opt_information.mat', xcorrOptDir, g);
        inputTmpname = sprintf('%s/xcorr_group_%d_opt_information_%s.mat', xcorrOptDir, g, uuid);
        save('-v7.3', inputTmpname, 'nF_g', 'max_xcorr_mat_g', 'absolute_shift_mat_g', ...
            'overlap_matrix_g', 'max_allow_shift_g', 'tileIdx_g');
        movefile(inputTmpname, inputFullpath);
        
        inputFullpaths{g} = inputFullpath;
        outputFullpaths{g} = sprintf('%s/xcorr_group_%d_opt_result.mat', xcorrOptDir, g);
        funcStrs{g} = sprintf(['stitch_global_grid_assignment_wrapper(''%s'',', ...
            '''%s'',%d,[%s],''%s'')'], inputFullpath, outputFullpaths{g}, xcorr_thresh, ...
            strrep(num2str(axisWeight, '%.20d,'), ' ', ''), uuid);
    else
        [d_shift_g] = stitch_global_grid_assignment(nF_g, max_xcorr_mat_g, absolute_shift_mat_g, ...
            overlap_matrix_g, max_allow_shift_g, xcorr_thresh, tileIdx_g, axisWeight);
        d_shift_pre(gind, :) = d_shift_g;
    end        
end

% submit jobs and collect results
if parseCluster
    ng_max = sum(grpIdx == mode(grpIdx));
    memAllocate = max(20, ng_max^2 * 0.0005) * nodeFactor;
    maxTrialNum_xcorr = 2;
    for i = 1 : 3
        is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'maxTrialNum', maxTrialNum_xcorr, 'parseCluster', parseCluster, ...
            'memAllocate', memAllocate * 2^(i - 1), 'mccMode', mccMode, 'configFile', configFile);
        if all(is_done_flag)
            break;
        end
    end
    
    if ~all(is_done_flag)
        error('group registartion cannot be fully done!')
    end
    for g = 1 : ng
        a = load(outputFullpaths{g}, 'd_shift_g');
        d_shift_g = a.d_shift_g;
        d_shift_pre(grpIdx == uniq_grp(g), :) = d_shift_g;        
    end
end

% if there is only one group, reduce to grid method
if ng == 1
    d_shift = round(d_shift_pre);
    d_shift = d_shift - d_shift(1, :);
    return;
end

% step 2: shifts across groups: only keep overlaps across groups
% define a max problem size, and only keep the strongest pairs across
% groups propotional to the number of pairs across groups

gn_i = max_xcorr_mat(all_gp_inds, 1);
gn_j = max_xcorr_mat(all_gp_inds, 2);
overlap_matrix(sub2ind([nF, nF], gn_i, gn_j)) = 0;
max_xcorr_mat(all_gp_inds, :) = [];
absolute_shift_mat(all_gp_inds, :) = [];
max_allow_shift(all_gp_inds, :) = [];

max_xcorr_mat_filt = max_xcorr_mat;
filt_inds = max_xcorr_mat_filt(:, 3) < xcorr_thresh;
max_xcorr_mat_filt(filt_inds, :) = [];
absolute_shift_mat_filt = absolute_shift_mat(~filt_inds, :);
max_allow_shift_filt = max_allow_shift(~filt_inds, :);
cuboid_overlap_mat_filt = cuboid_overlap_ij_mat(~filt_inds, :);

% infer the distances across groups
% the distances is the weighted averaged distances based on the absolute
% shifts across groups
group_overlap_matrix = false(ng, ng);
group_pair_sizes = zeros(ng, ng);
group_overlap_volumes = zeros(ng, ng);
grp_dist_mat = zeros(ng, ng, 3);
grp_dist_pre_mat = zeros(ng, ng, 3);
for gi = 1 : ng - 1
    grpIdx_i = gind_cell{gi};    
    for gj = i + 1 : ng
        grpIdx_j = gind_cell{gj};
        xcorr_pair_ij = ismember(max_xcorr_mat_filt(:, 1), grpIdx_i) & ismember(max_xcorr_mat_filt(:, 2), grpIdx_j);
        group_pair_sizes(gi, gj) = sum(xcorr_pair_ij);
        group_overlap_matrix(gi, gj) = any(overlap_matrix(grpIdx_i, grpIdx_j), 'all');
        if group_pair_sizes(gi, gj) == 0
            continue;
        end

        cuboid_overlap_mat_filt_ij = cuboid_overlap_mat_filt(xcorr_pair_ij, :);

        group_overlap_volumes(gi, gj) = sum(prod(cuboid_overlap_mat_filt_ij(:, 4 : 6) - cuboid_overlap_mat_filt_ij(:, 1 : 3), 2));

        g_w_ij = max_xcorr_mat_filt(xcorr_pair_ij, 3);
        grp_dist_mat(gi, gj, :) = sum(absolute_shift_mat_filt(xcorr_pair_ij, 3 : 5) .* g_w_ij, 1) / sum(g_w_ij);
        
        shift_pre_ij = d_shift_pre(max_xcorr_mat_filt(xcorr_pair_ij, 2), :) - d_shift_pre(max_xcorr_mat_filt(xcorr_pair_ij, 1), :);
        grp_dist_pre_mat(gi, gj, :) = sum(shift_pre_ij .* g_w_ij, 1) / sum(g_w_ij);
    end
end

% mapping between tile indices to the nodes in optimization problem

nP = sum(group_pair_sizes > 0, 'all');
% w = ones(nP, 1) * 0.00;
w = ones(nP, 1) * 0.01;

R_w = zeros(nP, ng);

g_pair_inds = find(group_pair_sizes);
[np_i, np_j] = ind2sub([ng, ng], g_pair_inds);
inds_i = sub2ind(size(R_w), 1 : nP, np_i');
inds_j = sub2ind(size(R_w), 1 : nP, np_j');
R_w(inds_i) = -1;
R_w(inds_j) = 1;

W = diag(w);
R_w = W.^0.5 * R_w;
% R_w = R;

% d_w = absolute_shift_mat_filt(:, 3 : 5);
d_w = reshape(grp_dist_mat, [], 3);
d_w = d_w(g_pair_inds, :);
d_w = W.^0.5 * d_w;

% inequalities
neq = sum(group_overlap_matrix(:));

R = zeros(neq, ng);

[n_i, n_j] = find(group_overlap_matrix);

inds_i = sub2ind(size(R), 1 : neq, n_i');
inds_j = sub2ind(size(R), 1 : neq, n_j');

R(inds_i) = -1; 
R(inds_j) = 1;

group_max_allow_shift = zeros(neq, 6);
for gi = 1 : ng - 1
    grpIdx_i = gind_cell{gi};    
    for gj = i + 1 : ng
        grpIdx_j = gind_cell{gj};
        xcorr_pair_ij = ismember(max_xcorr_mat_filt(:, 1), grpIdx_i) & ismember(max_xcorr_mat_filt(:, 2), grpIdx_j);
        if group_pair_sizes(gi, gj) == 0
            continue;
        end

        group_max_allow_shift(n_i == gi & n_j == gj, :) = [min(max_allow_shift_filt(xcorr_pair_ij, 1 : 3), [], 1), ...
            max(max_allow_shift_filt(xcorr_pair_ij, 4 : 6), [], 1)];
    end
end

% Aeq and beq: empty
Aeq = [];
beq = [];

d_shift_opt = zeros(ng, 3);
for i = 1 : 3
    C = R_w;
    d = d_w(:, i);
    A = [R; -R];

    l = group_max_allow_shift(:, i);
    u = group_max_allow_shift(:, 3 + i);
    b = [u; -l];

    [x, resnorm, residual, exitflag, output, lambda] = lsqlin(C, d, A, b, Aeq, beq);
    d_shift_opt(:, i) = x;
end

step 3: reconstruct the d_shift for all tiles
use MST to assign the absolute shifts for tiles from group to group
G = graph(np_i, np_j, -group_pair_sizes(g_pair_inds));
G_1 = graph(np_i, np_j, -group_overlap_volumes(g_pair_inds));
T = minspantree(G, 'type', 'forest');

%  use DFS search to just shift to the tile's precessor search from group 1
d_shift = zeros(nF, 3);
v = dfsearch(T, 1);
visit_flag = false(numel(v), 1);
for i = 1 : numel(v)
    n_i = v(i);
    if i == 1
        grpIdx_i = gind_cell{n_i};
        d_shift(grpIdx_i, :) = d_shift_pre(grpIdx_i, :);
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
    
    grpIdx_i = gind_cell{n_i};    
    % grpIdx_pr = gind_cell{n_pr};    
    
    grp_dist_opt = d_shift_opt(n_i, :) - d_shift_opt(n_pr, :);
    d_shift(grpIdx_i, :) = d_shift_pre(grpIdx_i, :) - st_sign .* squeeze(grp_dist_pre_mat(s, t, :))' + grp_dist_opt;
    visit_flag(n_i) = true;
end
for i = 1 : ng
    d_shift(grpIdx_i, :) = d_shift_pre(grpIdx_i, :) - 
end


% round to integers and normalize for the first tile.
d_shift = round(d_shift);
d_shift = d_shift - d_shift(1, :);

end


