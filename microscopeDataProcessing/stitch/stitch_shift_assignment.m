function [xyz_shift, d_shift] = stitch_shift_assignment(zarrFullpaths, xcorrDir, imSizes, xyz, ...
    px, xyz_factors, overlap_matrix, overlap_regions, MaxOffset, xcorrDownsample, xcorrThresh, tileIdx, assign_method, ...
    stitch2D, axisWeight, groupFile, parseCluster, nodeFactor)
% main function for stitch shift assignment 
% The main code is taken from XR_stitching_frame_zarr_dev_v1.m (to simplify
% the function).
%
% First compute xcorr between overlapping tiles, and then perform
% assignment of shifts either by local or global assignment method
% local assignment: MST and DFS 
% global assignment: weighted constrained linear least square
% 
% author: Xiongtao Ruan (11/05/2021)
% xruan (02/06/2022): for 'test' assignment method, not compute xcorr for tiles overlap from corners.
% xruan (05/09/2022): add support for 2d stitch
% xruan (05/15/2022): add support for group-based stitch
% xruan (05/21/2022): change absolute_shift_mat to finds pari, shifts NX5 
% xruan (06/20/2022): add input variable axisWeight for user defined weights for optimization
% xruan (08/03/2022): fix bug for max offset bounds based on the coordinate orders
% xruan (12/13/2022): change xcorr thresh as user defined parameter with default 0.25


fprintf('Compute cross-correlation based registration between overlap tiles...\n');

if isempty(axisWeight)
    axisWeight = [1, 0.1, 10]; % order y, x, z
end

xcorr_thresh = xcorrThresh;

if ~exist(xcorrDir, 'dir')
    mkdir_recursive(xcorrDir, true);
end

xf = xyz_factors(1);
yf = xyz_factors(2);
zf = xyz_factors(3);

nF = numel(zarrFullpaths);
absolute_shift_mat = zeros(nF * (nF - 1) / 2, 5); % order: x, y, z
max_xcorr_mat = zeros(nF * (nF - 1) / 2, 3);

% refactor the code with slurm generic framework
overlap_matrix_orig = overlap_matrix;
[ti, tj] = ind2sub(size(overlap_matrix), find(overlap_matrix));
switch assign_method
    case {'grid', 'test'}
        if ~all(tileIdx(:, 4) == tileIdx(1, 4))
            warning('The tiles are not from the same data folder, please make sure they are continuous grids when using grid method.');
        end
        
        % grid_overlap_inds = sum(abs(tileIdx(ti, :) - tileIdx(tj, :)), 2) == 1;
        corner_inds = sum(abs(tileIdx(ti, 1 : 3) - tileIdx(tj, 1 : 3)), 2) ~= 1;

        % remove non-grid overlaps (corner overlaps)
        ti_r = ti(corner_inds);
        tj_r = tj(corner_inds);
        overlap_matrix(sub2ind(size(overlap_matrix), ti_r, tj_r)) = 0; 
    case 'group'
        if isempty(groupFile)
            if all(tileIdx(:, 4) == tileIdx(1, 4))
                warning('The input group file is not provided, and the tiles are from the same data folder, treat all tiles as one group.');
            else
                sprintf('The input group file is not provided, treat tiles in each subfolder as a group.\n');
            end
            gIdx = tileIdx;
            gid = tileIdx(:, 4);
        else
            t = readtable(groupFile);
            if size(t, 2) == 4
               gIdx = [t.x, t.y, t.z, zeros(size(t, 1), 1)];
            else
               gIdx = [t.x, t.y, t.z, t.did];
            end
            gid = t.g;
        end
        
        [kidx, d] = knnsearch(gIdx, tileIdx, 'k', 1);
        
        if any(d ~= 0)
            warning('Some tiles are not included in the groups');
        end

        grpIdx = gid(kidx);
        uniqGrps = unique(grpIdx);
        
        for g = 1 : numel(uniqGrps)
            ginds = grpIdx == uniqGrps(g);
            ij_inds = ginds(ti) & ginds(tj); 
            ti_g = ti(ij_inds);
            tj_g = tj(ij_inds);

            corner_inds = sum(abs(tileIdx(ti_g, :) - tileIdx(tj_g, :)), 2) ~= 1;
            ti_g = ti_g(corner_inds);
            tj_g = tj_g(corner_inds);
            
            overlap_matrix(sub2ind(size(overlap_matrix), ti_g, tj_g)) = 0; 
        end
end
[ti, tj] = ind2sub(size(overlap_matrix), find(overlap_matrix));

inputFullpaths = zarrFullpaths(ti);
outputFullpaths = arrayfun(@(x) sprintf('%s/xcorr_tile_%d_tile_%d.mat', xcorrDir, ti(x), tj(x)), 1 : numel(ti), 'unif', 0);

cuboid_mat = [xyz, xyz + (imSizes(:, [2, 1, 3]) - 1) .* [xf, yf, zf] * px];

pinds = (ti - 1) * nF - ti .* (ti + 1) / 2 + tj;
cuboid_overlap_ij_mat = overlap_regions(pinds, :);

fprintf('Compute pairwise cross correlation between overlap tiles...\n');
if stitch2D
    funcStrs = arrayfun(@(x) sprintf(['multires_cross_correlation_registration_2d(''%s'',''%s'',''%s'',', ...
                                '[%s],[%s],[%s],%0.20d,[%s],''downSample'',[%s],''MaxOffset'',%s);toc;'], ...
                                zarrFullpaths{ti(x)}, zarrFullpaths{tj(x)}, outputFullpaths{x}, strrep(mat2str(cuboid_mat(ti(x), :)), ' ', ','), ...
                                strrep(mat2str(cuboid_mat(tj(x), :)), ' ', ','), strrep(mat2str(cuboid_overlap_ij_mat(x, :)), ' ', ','), ...
                                px, sprintf('%.20d;%.20d;%.20d', xf, yf, zf), strrep(num2str(xcorrDownsample, '%.20d,'), ' ', ''), ...
                                strrep(mat2str(MaxOffset), ' ', ',')), 1 : numel(ti), 'unif', 0);
else    
    funcStrs = arrayfun(@(x) sprintf(['multires_cross_correlation_registration_imblock_test(''%s'',''%s'',''%s'',', ...
                                '[%s],[%s],[%s],%0.20d,[%s],''downSample'',[%s],''MaxOffset'',%s);toc;'], ...
                                zarrFullpaths{ti(x)}, zarrFullpaths{tj(x)}, outputFullpaths{x}, strrep(mat2str(cuboid_mat(ti(x), :)), ' ', ','), ...
                                strrep(mat2str(cuboid_mat(tj(x), :)), ' ', ','), strrep(mat2str(cuboid_overlap_ij_mat(x, :)), ' ', ','), ...
                                px, sprintf('%.20d;%.20d;%.20d', xf, yf, zf), strrep(num2str(xcorrDownsample, '%.20d,'), ' ', ''), ...
                                strrep(mat2str(MaxOffset), ' ', ',')), 1 : numel(ti), 'unif', 0);
end

rawImageSizes = prod((cuboid_overlap_ij_mat(:, 4 : 6) - cuboid_overlap_ij_mat(:, 1 : 3))' ./ (px * [xf; yf; zf])) * 8 / 1024^3;
cpusPerTask_xcorr = min(24, max(1, prctile(ceil(rawImageSizes * 10 / 20), 90)) * nodeFactor);

maxTrialNum_xcorr = 2;
is_done_flag = slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, funcStrs, ...
    'cpusPerTask', cpusPerTask_xcorr, 'maxTrialNum', maxTrialNum_xcorr, 'parseCluster', parseCluster);    

maxTrialNum_xcorr = 2;    
if ~all(is_done_flag)
    is_done_flag = slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, funcStrs, ...
        'cpusPerTask', cpusPerTask_xcorr * 2, 'maxTrialNum', maxTrialNum_xcorr, 'parseCluster', parseCluster);
end

maxTrialNum_xcorr = 2;
if ~all(is_done_flag)
    is_done_flag = slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, funcStrs, ...
        'cpusPerTask', cpusPerTask_xcorr * 4, 'maxTrialNum', maxTrialNum_xcorr, 'parseCluster', parseCluster);
end

% collect results
if all(is_done_flag)
    d_w = zeros(numel(outputFullpaths), 6);
    for f = 1 : numel(outputFullpaths)
        xcorrFullpath = outputFullpaths{f};
        ind = pinds(f);
        i = ti(f);
        j = tj(f);

        if exist(xcorrFullpath, 'file')
            order_flag = (0.5 - (xyz(i, :) > xyz(j, :))) * 2;
            % order_flag = [1, 1, 1];
            a = load(xcorrFullpath);
            absolute_shift_mat(ind, :) = [i, j, a.relative_shift .* order_flag];
            max_xcorr_mat(ind, :) = [i, j, a.max_xcorr];
            d_w(f, :) = [i, j, a.max_xcorr, absolute_shift_mat(ind, 3 : 5)];
        end
    end
    absolute_shift_mat(absolute_shift_mat(:, 1) == 0 | absolute_shift_mat(:, 2) == 0, :) = [];
    max_xcorr_mat(max_xcorr_mat(:, 1) == 0 | max_xcorr_mat(:, 2) == 0, :) = [];
else
    error('Some xcorr files are missing!')
end

fprintf('Optimize registration for all tiles with %s method ...\n', assign_method);
fprintf('xcorr threshold: %f , %d / %d pairs are below the threashold.\n', ...
    xcorr_thresh, sum(max_xcorr_mat(:, 3) < xcorr_thresh), size(max_xcorr_mat, 1));

% perform assignment
switch assign_method
    case 'local'
        [d_shift] = stitch_local_assignment(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, xcorr_thresh);
    case 'global'
        neq = size(max_xcorr_mat, 1);
        max_shift_l = -ones(neq, 1) .* MaxOffset;
        max_shift_u = (cuboid_overlap_ij_mat(:, 4 : 6) - cuboid_overlap_ij_mat(:, 1 : 3)) ./ (px .* [xf, yf, zf]);
        max_shift_u = min(max_shift_u - 1, MaxOffset);
        
        % max_allow_shift = [max_shift_l, max_shift_u];
        order_mat = xyz(max_xcorr_mat(:, 2), :) >= xyz(max_xcorr_mat(:, 1), :);
        max_allow_shift = [max_shift_l .* order_mat + max_shift_u .* (order_mat - 1), max_shift_u .* order_mat + max_shift_l .* (order_mat - 1)];        
        
        [d_shift] = stitch_global_assignment(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, max_allow_shift, xcorr_thresh);
    case 'grid'
        neq = size(max_xcorr_mat, 1);
        max_shift_l = -ones(neq, 1) .* MaxOffset;
        max_shift_u = (cuboid_overlap_ij_mat(:, 4 : 6) - cuboid_overlap_ij_mat(:, 1 : 3)) ./ (px .* [xf, yf, zf]);
        % max_shift_u = min(max_shift_u - 1, MaxOffset);
        max_shift_u = min(max_shift_u - 3, MaxOffset);

        if stitch2D
            max_shift_u(:, 3) = 0;
        end
        
        % max_allow_shift = [max_shift_l, max_shift_u];
        order_mat = xyz(max_xcorr_mat(:, 2), :) >= xyz(max_xcorr_mat(:, 1), :);
        max_allow_shift = [max_shift_l .* order_mat + max_shift_u .* (order_mat - 1), max_shift_u .* order_mat + max_shift_l .* (order_mat - 1)];

        [d_shift] = stitch_global_grid_assignment(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, max_allow_shift, xcorr_thresh, tileIdx, axisWeight);
    case 'group'
        neq = size(max_xcorr_mat, 1);
        max_shift_l = -ones(neq, 1) .* MaxOffset;
        max_shift_u = (cuboid_overlap_ij_mat(:, 4 : 6) - cuboid_overlap_ij_mat(:, 1 : 3)) ./ (px .* [xf, yf, zf]);
        max_shift_u = min(max_shift_u - 1, MaxOffset);

        if stitch2D
            max_shift_u(:, 3) = 0;
        end
        
        % max_allow_shift = [max_shift_l, max_shift_u];
        order_mat = xyz(max_xcorr_mat(:, 2), :) >= xyz(max_xcorr_mat(:, 1), :);
        max_allow_shift = [max_shift_l .* order_mat + max_shift_u .* (order_mat - 1), max_shift_u .* order_mat + max_shift_l .* (order_mat - 1)];        
        % increase bound for inter-group tiles
        inter_ginds = grpIdx(max_xcorr_mat(:, 1)) ~= grpIdx(max_xcorr_mat(:, 2));
        max_allow_shift(inter_ginds, :) = [-ones(sum(inter_ginds), 1) .* MaxOffset, ones(sum(inter_ginds), 1) .* MaxOffset] * 2; 
        
        [d_shift] = stitch_global_group_assignment(nF, max_xcorr_mat, absolute_shift_mat, ...
            overlap_matrix, max_allow_shift, xcorr_thresh, tileIdx, axisWeight, grpIdx);        
    case 'test'
        neq = size(max_xcorr_mat, 1);
        max_shift_l = -ones(neq, 1) .* MaxOffset;
        max_shift_u = (cuboid_overlap_ij_mat(:, 4 : 6) - cuboid_overlap_ij_mat(:, 1 : 3)) ./ (px .* [xf, yf, zf]);
        max_shift_u = min(max_shift_u - 1, MaxOffset);
        
        max_allow_shift = [max_shift_l, max_shift_u];
        
        [d_shift] = stitch_global_grid_assignment_test(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, max_allow_shift, xcorr_thresh, tileIdx);
end

% for debug d_shift
if ~false
    absolute_shift_mat_1 = zeros(nF * (nF - 1) / 2, 3);
    for f = 1 : numel(ti)
        ind = pinds(f);
        i = ti(f);
        j = tj(f);
        absolute_shift_mat_1(ind, :) = d_shift(j, :) - d_shift(i, :);
    end
end

xyz_shift = xyz + d_shift .* [xf, yf, zf] .* px;

end


% unify input and output for different assignment methods

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

% calculate absolute shifts of tiles
% xyz_shift = xyz;
% for i = 1 : nF - 1
%     for j = i + 1 : nF
%         if ~overlap_matrix(i, j)
%             continue;
%         end
%         % ind = 0.5 * i * (2 * j - i - 1);
%         ind = (i - 1) * nF - i * (i + 1) / 2 + j;
%         if all(absolute_shift_mat(ind, :) == 0) 
%             continue;
%         end
%         % order_flag = (0.5 - (xyz(i, :) > xyz(j, :))) * 2;
%         % xyz_shift(j, :) = (xyz_shift(i, :) - xyz(i, :)) + (xyz_shift(j, :) + order_flag .* relative_shift_mat(ind, :) .* [xf, yf, zf] * px);
%         xyz_shift(j, :) = (xyz_shift(i, :) - xyz(i, :)) + (xyz_shift(j, :) + absolute_shift_mat(ind, :) .* [xf, yf, zf] * px);
%     end
% end
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


function [d_shift] = stitch_global_assignment(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, max_allow_shift, xcorr_thresh)
% solve weighted constrained linear least square problem for the assignment
% the weight is the function of max shift (currently just y=x). 
% max allowed shift is based on the maxShift parameter and the number of
% overlap between tiles.
% 02/03/2022: add threshold for objective to remove the samples with small xcorr

neq = sum(overlap_matrix(:));

R = zeros(neq, nF);

[n_i, n_j] = find(overlap_matrix);

inds_i = sub2ind(size(R), 1 : neq, n_i');
inds_j = sub2ind(size(R), 1 : neq, n_j');

R(inds_i) = -1; 
R(inds_j) = 1;

max_xcorr_mat_filt = max_xcorr_mat;
filt_inds = max_xcorr_mat_filt(:, 3) < xcorr_thresh;
max_xcorr_mat_filt(filt_inds, :) = [];
w = max_xcorr_mat_filt(:, 3);

nP = size(max_xcorr_mat_filt, 1);
R_w = zeros(nP, nF);

np_i = max_xcorr_mat_filt(:, 1);
np_j = max_xcorr_mat_filt(:, 2);
inds_i = sub2ind(size(R_w), 1 : nP, np_i');
inds_j = sub2ind(size(R_w), 1 : nP, np_j');
R_w(inds_i) = -1; 
R_w(inds_j) = 1;

W = diag(w);
R_w = W.^0.5 * R_w;
% R_w = R;

% c_inds = (np_i - 1) * nF - np_i .* (np_i + 1) / 2 + np_j;
% d_w = absolute_shift_mat(c_inds, :);
d_w = absolute_shift_mat(~filt_inds, 3 : 5);
d_w = W.^0.5 * d_w;

d_shift = zeros(nF, 3);
for i = 1 : 3
    C = R_w;
    d = d_w(:, i);
    A = [R; -R];

    l = max_allow_shift(:, i);
    u = max_allow_shift(:, 3 + i);
    b = [u; -l];

    [x,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b);
    d_shift(:, i) = x;
end

% round to integers and normalize for the first tile.
d_shift = round(d_shift);
d_shift = d_shift - d_shift(1, :);

end


function [d_shift] = stitch_global_grid_assignment(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, max_allow_shift, xcorr_thresh, tileIdx, axisWeight)
% solve weighted constrained linear least square problem for the assignment
% the weight is the function of max shift (currently just y=x). 
% max allowed shift is based on the maxShift parameter and the number of
% overlap between tiles.
% 02/03/2022: add threshold for objective to remove the samples with small xcorr

neq = sum(overlap_matrix(:));

R = zeros(neq, nF);

[n_i, n_j] = find(overlap_matrix);

inds_i = sub2ind(size(R), 1 : neq, n_i');
inds_j = sub2ind(size(R), 1 : neq, n_j');

R(inds_i) = -1; 
R(inds_j) = 1;

max_xcorr_mat_filt = max_xcorr_mat;
filt_inds = max_xcorr_mat_filt(:, 3) < xcorr_thresh;
max_xcorr_mat_filt(filt_inds, :) = [];
% w = max_xcorr_mat_filt(:, 3);
nP = size(max_xcorr_mat_filt, 1);
w = ones(nP, 1) * 0.00;
for i = 1 : nP
    s = max_xcorr_mat_filt(i, 1);
    t = max_xcorr_mat_filt(i, 2);
    
    % check if two tiles are neighboring tiles
    if sum(abs(tileIdx(s, 1 : 3) - tileIdx(t, 1 : 3))) == 1
        aind = find(abs(tileIdx(s, 1 : 3) - tileIdx(t, 1 : 3)));
        switch aind
            case 1
                % x
                w(i) = axisWeight(2);
            case 2
                % y
                w(i) = axisWeight(1);
            case 3
                % z
                w(i) = axisWeight(3);
        end
    end    
end

R_w = zeros(nP, nF);

np_i = max_xcorr_mat_filt(:, 1);
np_j = max_xcorr_mat_filt(:, 2);
inds_i = sub2ind(size(R_w), 1 : nP, np_i');
inds_j = sub2ind(size(R_w), 1 : nP, np_j');
R_w(inds_i) = -1; 
R_w(inds_j) = 1;

W = diag(w);
R_w = W.^0.5 * R_w;
% R_w = R;

% c_inds = (np_i - 1) * nF - np_i .* (np_i + 1) / 2 + np_j;
d_w = absolute_shift_mat(~filt_inds, 3 : 5);
d_w = W.^0.5 * d_w;

d_shift = zeros(nF, 3);
for i = 1 : 3
    C = R_w;
    d = d_w(:, i);
    A = [R; -R];

    l = max_allow_shift(:, i);
    u = max_allow_shift(:, 3 + i);
    b = [u; -l];

    [x,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b);
    d_shift(:, i) = x;
end

% round to integers and normalize for the first tile.
d_shift = round(d_shift);
d_shift = d_shift - d_shift(1, :);

end


function [d_shift] = stitch_global_group_assignment(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, max_allow_shift, xcorr_thresh, tileIdx, axisWeight, grpIdx)
% Two-step optimization: 
% step 1: grid method for tiles within each group
% step 2: weighted constrained linear least square for shifts across groups
% with the results from the first step as constrains. 
% 
% the weight is the function of max shift (currently just y=x). 
% max allowed shift is based on the maxShift parameter and the number of
% overlap between tiles.


% step 1: grid method for tiles within the same group

uniq_grp = unique(grpIdx);
ng = numel(uniq_grp);

all_gp_inds = false(size(max_xcorr_mat, 1), 1);
d_shift_pre = zeros(nF, 3);
for g = 1 : ng
    gind = find(grpIdx == uniq_grp(g));
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

    [d_shift_g] = stitch_global_grid_assignment(nF_g, max_xcorr_mat_g, absolute_shift_mat_g, ...
        overlap_matrix_g, max_allow_shift_g, xcorr_thresh, tileIdx_g, axisWeight);
    d_shift_pre(gind, :) = d_shift_g;
end

% if there is only one group, reduce to grid method
if ng == 1
    d_shift = round(d_shift_pre);
    d_shift = d_shift - d_shift(1, :);
    return;
end

% step 2: shifts across groups: only keep overlaps across groups

gn_i = max_xcorr_mat(all_gp_inds, 1);
gn_j = max_xcorr_mat(all_gp_inds, 2);
overlap_matrix(sub2ind([nF, nF], gn_i, gn_j)) = 0;
max_xcorr_mat(all_gp_inds, :) = [];
absolute_shift_mat(all_gp_inds, :) = [];
max_allow_shift(all_gp_inds, :) = [];

neq = sum(overlap_matrix(:));

R = zeros(neq, nF);

[n_i, n_j] = find(overlap_matrix);

inds_i = sub2ind(size(R), 1 : neq, n_i');
inds_j = sub2ind(size(R), 1 : neq, n_j');

R(inds_i) = -1; 
R(inds_j) = 1;

max_xcorr_mat_filt = max_xcorr_mat;
filt_inds = max_xcorr_mat_filt(:, 3) < xcorr_thresh;
max_xcorr_mat_filt(filt_inds, :) = [];
% w = max_xcorr_mat_filt(:, 3);
nP = size(max_xcorr_mat_filt, 1);
% w = ones(nP, 1) * 0.00;
w = ones(nP, 1) * 0.01;

R_w = zeros(nP, nF);

np_i = max_xcorr_mat_filt(:, 1);
np_j = max_xcorr_mat_filt(:, 2);
inds_i = sub2ind(size(R_w), 1 : nP, np_i');
inds_j = sub2ind(size(R_w), 1 : nP, np_j');
R_w(inds_i) = -1; 
R_w(inds_j) = 1;

W = diag(w);
R_w = W.^0.5 * R_w;
% R_w = R;

d_w = absolute_shift_mat(~filt_inds, 3 : 5);
d_w = W.^0.5 * d_w;

% Aeq and beq: distances from tiles within the same groups.
nE = numel(gn_i); 
Aeq = zeros(nE, nF);
einds_i = sub2ind(size(Aeq), 1 : nE, gn_i');
einds_j = sub2ind(size(Aeq), 1 : nE, gn_j');
Aeq(einds_i) = -1; 
Aeq(einds_j) = 1;
Beq = d_shift_pre(gn_j, :) - d_shift_pre(gn_i, :);

d_shift = zeros(nF, 3);
for i = 1 : 3
    C = R_w;
    d = d_w(:, i);
    A = [R; -R];

    l = max_allow_shift(:, i);
    u = max_allow_shift(:, 3 + i);
    b = [u; -l];

    beq = Beq(:, i);

    [x,resnorm,residual,exitflag,output,lambda] = lsqlin(C, d, A, b, Aeq, beq);
    d_shift(:, i) = x;
end

% round to integers and normalize for the first tile.
d_shift = round(d_shift);
d_shift = d_shift - d_shift(1, :);

end


function [d_shift] = stitch_global_grid_assignment_test(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, max_allow_shift, xcorr_thresh, tileIdx)
% solve weighted constrained linear least square problem for the assignment
% the weight is the function of max shift (currently just y=x). 
% max allowed shift is based on the maxShift parameter and the number of
% overlap between tiles.
% 02/03/2022: add threshold for objective to remove the samples with small xcorr

neq = sum(overlap_matrix(:));

R = zeros(neq, nF);

[n_i, n_j] = find(overlap_matrix);

inds_i = sub2ind(size(R), 1 : neq, n_i');
inds_j = sub2ind(size(R), 1 : neq, n_j');

R(inds_i) = -1; 
R(inds_j) = 1;

max_xcorr_mat_filt = max_xcorr_mat;
max_xcorr_mat_filt(max_xcorr_mat_filt(:, 3) < xcorr_thresh, :) = [];
% w = max_xcorr_mat_filt(:, 3);
nP = size(max_xcorr_mat_filt, 1);
w = ones(nP, 1) * 0.00;
for i = 1 : nP
    s = max_xcorr_mat_filt(i, 1);
    t = max_xcorr_mat_filt(i, 2);
    
    % check if two tiles are neighboring tiles
    if sum(abs(tileIdx(s, :) - tileIdx(t, :))) == 1
        aind = find(abs(tileIdx(s, :) - tileIdx(t, :)));
        switch aind
            case 1
                w(i) = 1;
            case 2
                w(i) = 10;
            case 3
                w(i) = 0.1;
        end
    end    
end

R_w = zeros(nP, nF);

np_i = max_xcorr_mat_filt(:, 1);
np_j = max_xcorr_mat_filt(:, 2);
inds_i = sub2ind(size(R_w), 1 : nP, np_i');
inds_j = sub2ind(size(R_w), 1 : nP, np_j');
R_w(inds_i) = -1; 
R_w(inds_j) = 1;

W = diag(w);
R_w = W.^0.5 * R_w;
% R_w = R;

c_inds = (np_i - 1) * nF - np_i .* (np_i + 1) / 2 + np_j;
d_w = absolute_shift_mat(c_inds, :);
d_w = W.^0.5 * d_w;

d_shift = zeros(nF, 3);
for i = 1 : 3
    C = R_w;
    d = d_w(:, i);
    A = [R; -R];

    l = max_allow_shift(:, i);
    u = max_allow_shift(:, 3 + i);
    b = [u; -l];

    [x,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b);
    d_shift(:, i) = x;
end

% round to integers and normalize for the first tile.
d_shift = round(d_shift);
d_shift = d_shift - d_shift(1, :);

end

