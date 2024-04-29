function [xyz_shift, dxyz_shift] = stitch_shift_assignment(zarrFullpaths, xcorrDir, imSizes, xyz, ...
    px, xyz_factors, overlap_matrix, overlap_regions, MaxOffset, xcorrDownsample, xcorrThresh, tileIdx, assign_method, ...
    stitch2D, axisWeight, groupFile, largeZarr, poolSize, parseCluster, nodeFactor, mccMode, ConfigFile)
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
% xruan (04/07/2023): put xcorr results to subfolders if greater than 10k
% xruan (04/15/2023): combine xcorrs by the leading tile indices. 
% xruan (05/01/2023): add support for large zarr xcorr. 
% xruan (09/06/2023): change output d_shift to dxyz_shift to include pixel size info 


fprintf('Compute cross-correlation based registration between overlap tiles...\n');
t0 = tic();

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

% compute MIP slabs for large zarr files
mipDirStr = '';
% poolSize = [];
if largeZarr
    axis = [1, 1, 1];
    BatchSize = [1024, 1024, 1024];
    % poolSize = [5, 5, 5];
    if numel(poolSize) == 3
        mipDirStr = sprintf('MIP_Slabs_pooling_%d_%d_%d', poolSize(1), poolSize(2), poolSize(3));
    else
        mipDirStr = sprintf('MIP_Slabs_pooling_%d_%d_%d_%d_%d_%d', poolSize(1), ...
            poolSize(2), poolSize(3), poolSize(4), poolSize(5), poolSize(6));
    end
    mipSlab = true;
    parseParfor = false;
    
    [dataPaths, fsns] = fileparts(zarrFullpaths);
    uniq_dataPaths = unique(dataPaths);
    for d = 1 : numel(uniq_dataPaths)
        mipPath = [uniq_dataPaths{d}, '/', mipDirStr];
        if ~exist(mipPath, 'dir')
            mkdir(mipPath);
        end
    end

    funcStrs = cell(nF, 1);
    outputFullpaths = cell(nF, 1);
    for f = 1 : nF
        funcStrs{f} = sprintf(['XR_MIP_zarr(''%s'',''mipDirStr'',''%s'',''axis'',%s,', ...
            '''BatchSize'',%s,''poolSize'',%s,''mipSlab'',%s,''parseCluster'',%s,', ...
            '''parseParfor'',%s,''mccMode'',%s,''ConfigFile'',''%s'')'], ...
            zarrFullpaths{f}, mipDirStr, strrep(mat2str(axis), ' ', ','), strrep(mat2str(BatchSize), ' ', ','), ...
            strrep(mat2str(poolSize), ' ', ','), string(mipSlab), string(parseCluster), ...
            string(parseParfor), string(mccMode), ConfigFile);
        outputFullpaths{f} = sprintf('%s/%s/%s_z.zarr', dataPaths{f}, mipDirStr, fsns{f});
    end

    inputFullpaths = zarrFullpaths;
    memAllocate = prod(BatchSize) * 4 / 2^30 * 3;
    minTaskJobNum = round(nF / 2);
    for i = 1 : 3
        is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'cpusPerTask', 1, 'maxTrialNum', 2, 'parseCluster', parseCluster, ...
            'memAllocate', memAllocate * 2^(i - 1), 'minTaskJobNum', minTaskJobNum, ...
            'mccMode', mccMode, 'ConfigFile', ConfigFile);
        if all(is_done_flag)
            break;
        end
    end
end

% pairwise cross correlation computing
inputFullpaths = zarrFullpaths;

batchSize = 10000;
if nF > batchSize
    numBatch = ceil(nF ./ batchSize);    
    for i = 1 : numBatch
        s = (i - 1) * batchSize + 1;
        t = min(i * batchSize, nF);        
        xcorrSubdir = sprintf('%s/xcorr_tiles_%06d_%06d/', xcorrDir, s, t);
        if ~exist(xcorrSubdir, 'dir')
            mkdir(xcorrSubdir);
        end
        outputFullpaths(s : t) = arrayfun(@(x) sprintf('%s/xcorr_tile_%d.mat', xcorrSubdir, x), s : t, 'unif', 0);
    end
else
    outputFullpaths = arrayfun(@(x) sprintf('%s/xcorr_tile_%d.mat', xcorrDir, x), 1 : nF, 'unif', 0);
end

[ti, tj] = ind2sub(size(overlap_matrix), find(overlap_matrix));
cuboid_mat = [xyz, xyz + (imSizes(:, [2, 1, 3]) - 1) .* [xf, yf, zf] * px];
% for large zarr xcorr, don't group xcorr calculation by tiles
if largeZarr
    inputFullpaths = zarrFullpaths(ti);
    outputFullpaths = arrayfun(@(x) sprintf('%s/xcorr_tile_%d_%d.mat', xcorrDir, ti(x), tj(x)), 1 : numel(ti), 'unif', 0);

    funcStrs = cell(numel(ti), 1);
    for f = 1 : numel(ti)
        ti_f = ti(f);
        tj_f = tj(f);
        pair_indices = [ti_f, tj_f];
        
        pinds_f = (ti_f - 1) * nF - ti_f .* (ti_f + 1) / 2 + tj_f;
        cuboid_overlap_ij_mat_f = overlap_regions(pinds_f, :);
        zarrFullpaths_j_str = sprintf('{''%s''}', strjoin(zarrFullpaths(tj_f), ''','''));
    
        funcStrs{f} = sprintf(['cross_correlation_registration_wrapper(''%s'',%s,''%s'',', ...
            '[%s],[%s],[%s],[%s],%0.20d,[%s],''Stitch2D'',%s,''downSample'',[%s],', ...
            '''MaxOffset'',%s,''largeZarr'',%s,''mipDirStr'',''%s'',''poolSize'',%s,', ...
            '''parseCluster'',%s,''mccMode'',%s,''ConfigFile'',''%s'')'], ...
            zarrFullpaths{ti_f}, zarrFullpaths_j_str, outputFullpaths{f}, strrep(mat2str(pair_indices), ' ', ','), ...
            strrep(mat2str(cuboid_mat(ti_f, :)), ' ', ','), strrep(mat2str(cuboid_mat(tj_f, :)), ' ', ','), ...
            strrep(mat2str(cuboid_overlap_ij_mat_f), ' ', ','), px, sprintf('%.20d,%.20d,%.20d', xf, yf, zf), ...
            string(stitch2D), strrep(num2str(xcorrDownsample, '%.20d,'), ' ', ''), strrep(mat2str(MaxOffset), ' ', ','), ...
            string(largeZarr), mipDirStr, strrep(mat2str(poolSize), ' ', ','), string(parseCluster), ...
            string(mccMode), ConfigFile);
    end
    minTaskJobNum = numel(ti);
else
    funcStrs = cell(nF, 1);
    for f = 1 : nF
        inds_f = ti == f;
        if ~any(inds_f)
            continue;
        end
        ti_f = ti(inds_f);
        tj_f = tj(inds_f);
        pair_indices = [ti_f, tj_f];
        
        pinds_f = (ti_f - 1) * nF - ti_f .* (ti_f + 1) / 2 + tj_f;
        cuboid_overlap_ij_mat_f = overlap_regions(pinds_f, :);
        zarrFullpaths_j_str = sprintf('{''%s''}', strjoin(zarrFullpaths(tj_f), ''','''));
    
        funcStrs{f} = sprintf(['cross_correlation_registration_wrapper(''%s'',%s,''%s'',', ...
            '[%s],[%s],[%s],[%s],%0.20d,[%s],''Stitch2D'',%s,''downSample'',[%s],', ...
            '''MaxOffset'',%s,''largeZarr'',%s,''mipDirStr'',''%s'',''poolSize'',%s,', ...
            '''parseCluster'',%s,''mccMode'',%s,''ConfigFile'',''%s'')'], ...
            zarrFullpaths{f}, zarrFullpaths_j_str, outputFullpaths{f}, strrep(mat2str(pair_indices), ' ', ','), ...
            strrep(mat2str(cuboid_mat(f, :)), ' ', ','), strrep(mat2str(cuboid_mat(tj_f, :)), ' ', ','), ...
            strrep(mat2str(cuboid_overlap_ij_mat_f), ' ', ','), px, sprintf('%.20d,%.20d,%.20d', xf, yf, zf), ...
            string(stitch2D), strrep(num2str(xcorrDownsample, '%.20d,'), ' ', ''), strrep(mat2str(MaxOffset), ' ', ','), ...
            string(largeZarr), mipDirStr, strrep(mat2str(poolSize), ' ', ','), string(parseCluster), ...
            string(mccMode), ConfigFile);
    end
    minTaskJobNum = 1;
end

include_inds = ~cellfun(@isempty, funcStrs);
inputFullpaths = inputFullpaths(include_inds);
outputFullpaths = outputFullpaths(include_inds);
funcStrs = funcStrs(include_inds);

fprintf('Compute pairwise cross correlation between overlap tiles...\n');

sz = getImageSize(zarrFullpaths{1});
pinds = (ti - 1) * nF - ti .* (ti + 1) / 2 + tj;
cuboid_overlap_ij_mat = overlap_regions(pinds, :);
rawImageSizes = prod(min((cuboid_overlap_ij_mat(:, 4 : 6) - cuboid_overlap_ij_mat(:, 1 : 3))' ./ (px * [xf; yf; zf]) + MaxOffset(:), sz(:))) * 4 / 1024^3;
memAllocate = prctile(ceil(rawImageSizes) * (3.5 + 50 / prod(xcorrDownsample)), 99) * nodeFactor;
cpusPerTask_xcorr = 2;
maxTrialNum_xcorr = 2;

for i = 1 : 3
    is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', cpusPerTask_xcorr * 2^(i - 1), 'maxTrialNum', maxTrialNum_xcorr, ...
        'parseCluster', parseCluster, 'memAllocate', memAllocate * 2^(i - 1), ...
        'minTaskJobNum', minTaskJobNum, 'mccMode', mccMode, 'ConfigFile', ConfigFile);
    if all(is_done_flag)
        break;
    end
end

% collect results
if all(is_done_flag)
    d_w = zeros(nF * (nF - 1) / 2, 6);
    for f = 1 : numel(outputFullpaths)
        xcorrFullpath = outputFullpaths{f};
        a = load(xcorrFullpath);
        relative_shift_mat_f = a.relative_shift_mat;
        max_xcorr_mat_f = a.max_xcorr_mat;
        pair_indices_f = a.pair_indices;

        ti_f = pair_indices_f(:, 1);
        tj_f = pair_indices_f(:, 2);
        pinds_f = (ti_f - 1) * nF - ti_f .* (ti_f + 1) / 2 + tj_f;

        order_flag = (0.5 - (xyz(ti_f, :) > xyz(tj_f, :))) * 2;
        % order_flag = [1, 1, 1];
        absolute_shift_mat(pinds_f, :) = [ti_f, tj_f, relative_shift_mat_f .* order_flag];
        max_xcorr_mat(pinds_f, :) = [ti_f, tj_f, max_xcorr_mat_f];
        d_w(pinds_f, :) = [max_xcorr_mat(pinds_f, :), absolute_shift_mat(pinds_f, 3 : 5)];
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
            overlap_matrix, max_allow_shift, xcorr_thresh, tileIdx, axisWeight, ...
            grpIdx, cuboid_overlap_ij_mat, xcorrDir, parseCluster, nodeFactor, mccMode, ConfigFile);        
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

dxyz_shift = d_shift .* [xf, yf, zf] .* px;
xyz_shift = xyz + dxyz_shift;

fprintf('xcorr registration is done!\n')
toc(t0);
fprintf('\n');

end


% unify input and output for different assignment methods

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

