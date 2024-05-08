function [d_shift] = stitch_global_group_assignment(nF, max_xcorr_mat, absolute_shift_mat, ...
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
    if nF_g == 1
        continue;
    end

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
    memAllocate = max(20, ng_max^2 * 0.0001) * nodeFactor;
    maxTrialNum_xcorr = 2;
    minTaskJobNum = 1;
    if ng >= 4
        minTaskJobNum = max(2, round(ng / 2));
    end
    if ng_max > 1000
        minTaskJobNum = ng;
    end
    for i = 1 : 3
        is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'maxTrialNum', maxTrialNum_xcorr, 'parseCluster', parseCluster, ...
            'memAllocate', memAllocate * 2^(i - 1), 'minTaskJobNum', minTaskJobNum, ...
            'mccMode', mccMode, 'configFile', configFile);
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

% infer the distances across groups
% treat the overall overlap region between two groups as a super node
% the distances is the weighted averaged distances based on the absolute
% shifts across groups
group_overlap_matrix = false(ng, ng);
group_pair_sizes = zeros(ng, ng);
grp_dist_mat = zeros(ng, ng, 3);
grp_pair_centroids = zeros(ng, ng, 6);
for gi = 1 : ng - 1
    grpIdx_i = gind_cell{gi};    
    for gj = gi + 1 : ng
        grpIdx_j = gind_cell{gj};
        xcorr_pair_ij = ismember(max_xcorr_mat_filt(:, 1), grpIdx_i) & ismember(max_xcorr_mat_filt(:, 2), grpIdx_j);
        group_pair_sizes(gi, gj) = sum(xcorr_pair_ij);
        group_overlap_matrix(gi, gj) = any(overlap_matrix(grpIdx_i, grpIdx_j), 'all');
        if group_pair_sizes(gi, gj) == 0
            continue;
        end

        g_w_ij = max_xcorr_mat_filt(xcorr_pair_ij, 3);
        grp_dist_mat(gi, gj, :) = sum(absolute_shift_mat_filt(xcorr_pair_ij, 3 : 5) .* g_w_ij, 1) / sum(g_w_ij);
        
        grp_pair_centroids(gi, gj, :) = [mean(d_shift_pre(unique(max_xcorr_mat_filt(xcorr_pair_ij, 1)), :), 1), mean(d_shift_pre(unique(max_xcorr_mat_filt(xcorr_pair_ij, 2)), :), 1)];
    end
end

% construct the optimization problem

nP = sum(group_pair_sizes > 0, 'all');
n_node = nP * 2;
% w = ones(nP, 1) * 0.00;
w = ones(nP, 1) * 0.01;

R_w = zeros(nP, n_node);

g_pair_inds = find(group_pair_sizes);
[np_i, np_j] = ind2sub([ng, ng], g_pair_inds);
g_node_pairs = [(1 : nP)', (nP + 1 : 2 * nP)', np_i, np_j];
inds_i = sub2ind(size(R_w), 1 : nP, (1 : nP));
inds_j = sub2ind(size(R_w), 1 : nP, (nP + 1 : 2 * nP));
R_w(inds_i) = -1;
R_w(inds_j) = 1;

W = diag(w);
R_w = W.^0.5 * R_w;
% R_w = R;

% d_w = absolute_shift_mat_filt(:, 3 : 5);
d_w = reshape(grp_dist_mat, [], 3);
d_w = d_w(g_pair_inds, :);
d_w = W.^0.5 * d_w;

% inequalities A
neq = nP;

R = zeros(neq, n_node);

[n_i, n_j] = find(group_pair_sizes);

inds_i = sub2ind(size(R), 1 : neq, 1 : nP);
inds_j = sub2ind(size(R), 1 : neq, nP + 1 : 2 * nP);

R(inds_i) = -1; 
R(inds_j) = 1;

group_max_allow_shift = zeros(neq, 6);
for gi = 1 : ng - 1
    grpIdx_i = gind_cell{gi};    
    for gj = gi + 1 : ng
        grpIdx_j = gind_cell{gj};
        xcorr_pair_ij = ismember(max_xcorr_mat_filt(:, 1), grpIdx_i) & ismember(max_xcorr_mat_filt(:, 2), grpIdx_j);
        if group_pair_sizes(gi, gj) == 0
            continue;
        end

        group_max_allow_shift(n_i == gi & n_j == gj, :) = [min(max_allow_shift_filt(xcorr_pair_ij, 1 : 3), [], 1), ...
            max(max_allow_shift_filt(xcorr_pair_ij, 4 : 6), [], 1)];
    end
end

% Aeq and beq: distances for the centroids of nodes within a group
Aeq_cell = cell(ng, 1);
Beq_cell = cell(ng, 1);
node_centroids = zeros(ng, 3);
for g = 1 : ng
    node_g = unique([g_node_pairs(g_node_pairs(:, 3) == g, 1); g_node_pairs(g_node_pairs(:, 4) == g, 2)]);
    n_node_g = numel(node_g);
    if n_node_g == 0
        continue;
    end

    nE_g = n_node_g * (n_node_g - 1) / 2;
    Aeq_g = [];
    Beq_g = [];
    if nE_g > 0
        Aeq_g = zeros(nE_g, n_node);
        [gn_i, gn_j] = find(triu(ones(n_node_g), 1));
        
        einds_i = sub2ind(size(Aeq_g), 1 : nE_g, node_g(gn_i)');
        einds_j = sub2ind(size(Aeq_g), 1 : nE_g, node_g(gn_j)');
        Aeq_g(einds_i) = -1; 
        Aeq_g(einds_j) = 1;
    end
    
    node_centroids_g = zeros(n_node_g, 3);
    for k = 1 : n_node_g
        node_g_k = node_g(k);
        g_node_pair_k = g_node_pairs(any(g_node_pairs(:, 1 : 2) == node_g_k, 2), :);
        node_centroids_g_k = grp_pair_centroids(g_node_pair_k(3), g_node_pair_k(4), :);
        if node_g_k <= nP
            node_centroids_g(k, :) = node_centroids_g_k(1 : 3);
        else
            node_centroids_g(k, :) = node_centroids_g_k(4 : 6);
        end
    end
    if nE_g > 0
        Beq_g = node_centroids_g(gn_j, :) - node_centroids_g(gn_i, :);
    end

    Aeq_cell{g} = Aeq_g;
    Beq_cell{g} = Beq_g;
    node_centroids(node_g, :) = node_centroids_g;
end

Aeq = cat(1, Aeq_cell{:});
Beq = cat(1, Beq_cell{:});

% solve the optimization problem
d_shift_opt = zeros(n_node, 3);
for i = 1 : 3
    C = R_w;
    d = d_w(:, i);
    A = [R; -R];

    l = group_max_allow_shift(:, i);
    u = group_max_allow_shift(:, 3 + i);
    b = [u; -l];
    
    if isempty(Beq)
        beq = [];        
    else
        beq = Beq(:, i);
    end

    [x, resnorm, residual, exitflag, output, lambda] = lsqlin(C, d, A, b, Aeq, beq);
    d_shift_opt(:, i) = x;
end

% step 3: reconstruct the d_shift for all tiles
d_shift = zeros(nF, 3);
for g = 1 : ng
    grpIdx_i = gind_cell{g};
    node_g = unique([g_node_pairs(g_node_pairs(:, 3) == g, 1); g_node_pairs(g_node_pairs(:, 4) == g, 2)]);

    d_shift(grpIdx_i, :) = d_shift_pre(grpIdx_i, :) - node_centroids(node_g(1), :) + d_shift_opt(node_g(1), :);
end

% round to integers and normalize for the first tile.
d_shift = round(d_shift);
d_shift = d_shift - d_shift(1, :);

end

