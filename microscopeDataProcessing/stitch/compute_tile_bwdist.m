function [] = compute_tile_bwdist(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, singleDistMap, blockSize, shardSize, compressor, distBbox, dataOrderMat, Overwrite)
% compute distance transform for a tile after removing overlap regions.
% 
% Author: Xiongtao Ruan (10/29/2020)
% xruan (11/02/2020): only compute distance in overlap regions, and also
% only consider distance in xy plane.
% xruan (11/19/2022): compute feather power within the distance map to save
% time for stitching processing


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInfoFullname', @ischar);
ip.addRequired('tileInd', @isscalar);
ip.addRequired('bwdistFullpath', @ischar);
ip.addRequired('weightDegree', @isscalar);
ip.addRequired('singleDistMap', @islogical);
ip.addRequired('blockSize', @isvector);
ip.addRequired('shardSize', @(x) isempty(x) || isvector(x));
ip.addRequired('compressor', @ischar);
ip.addRequired('distBbox', @(x) isempty(x) || isvector(x));
ip.addRequired('dataOrderMat', @(x) isvector(x));
ip.addOptional('Overwrite', false, @islogical);

ip.parse(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, singleDistMap, blockSize, shardSize, compressor, distBbox, dataOrderMat, Overwrite);

pr = ip.Results;
Overwrite = pr.Overwrite;

uuid = get_uuid();

if exist(bwdistFullpath, 'dir') 
    return;
end

persistent tileFns overlap_matrix ol_region_cell overlap_map_mat bkFn datenum
dir_info = dir(blockInfoFullname);
if ~(strcmpi(bkFn, blockInfoFullname) && datenum == dir_info.datenum) || isempty(tileFns)
    bkFn = blockInfoFullname;
    a = load(blockInfoFullname, 'tileFns', 'overlap_matrix', 'ol_region_cell', 'overlap_map_mat');
    tileFns = a.tileFns;
    overlap_matrix = a.overlap_matrix;
    ol_region_cell = a.ol_region_cell;
    overlap_map_mat = a.overlap_map_mat;
    dir_info = dir(bkFn);
    datenum = dir_info.datenum;
    clear a;
end

nF = numel(tileFns);
i = tileInd;

im_i = readzarr(tileFns{i});
if dataOrderMat(3) ~= 3
    im_i = permute(im_i, dataOrderMat);
end

im_i_orig = im_i ~= 0;
im_i([1, end], :, :) = 0;
im_i(:, [1, end], :) = 0;
if ndims(im_i) == 3
    im_i(:, :, [1, end]) = 0;
end
im_i = im_i == 0;

% remove isolated empty pixels 
% stat = regionprops3(im_i, 'volume');
% if numel(stat.Volume) > 1 && any(stat.Volume < 5)
%     im_i = bwareaopen(im_i, 5);
% end

sz = size(im_i, 1 : 3);
if singleDistMap
    im_i_c = im_i(:, :, round((sz(3) + 1) / 2));
    im_dist_c = fastPower((bwdist(im_i_c) + 1) / 10, weightDegree);
    im_dist = repmat(im_dist_c, 1, 1, sz(3));
    same_z_inds = squeeze(all(im_i == im_i_c, [1, 2]));
    for z = 1 : sz(3)
        if same_z_inds(z)
            continue;
        end
        im_i_z = im_i(:, :, z);
        im_dist(:, :, z) = fastPower((bwdist(im_i_z) + 1) / 10, weightDegree);            
    end
else
    im_dist = ones(sz, 'single');
    for j = 1 : nF
        if ~overlap_matrix(i, j) && ~overlap_matrix(j, i)
            continue;
        end
        if i < j
            ind = (j - 1) * nF + i;
            mregion = ol_region_cell{overlap_map_mat(:, 2) == ind}{1};
        else
            ind = (i - 1) * nF + j;
            mregion = ol_region_cell{overlap_map_mat(:, 2) == ind}{2};
        end

        midx = mregion;
        % dsr_ol = dsr(mregion_f(2, 1) : mregion_f(2, 2), mregion_f(1, 1) : mregion_f(1, 2), mregion_f(3, 1) : mregion_f(3, 2));
        % dsr(midx(2, 1) : midx(2, 2), midx(1, 1) : midx(1, 2), midx(3, 1) : midx(3, 2)) = 0;

        mbbox = [midx(2, 1), midx(1, 1), midx(3, 1), midx(2, 2), midx(1, 2), midx(3, 2)];
        if any(mbbox(4 : 6) < mbbox(1 : 3))
            continue;
        end

        mbbox_pad = [max(mbbox(1 : 3) - 1, 1), min(mbbox(4 : 6) + 1, sz)];
        try 
            im_ij = crop3d_mex(im_i, mbbox_pad);
        catch ME
            disp(ME)
            im_ij = im_i(mbbox_pad(1) : mbbox_pad(4), mbbox_pad(2) : mbbox_pad(5), mbbox_pad(3) : mbbox_pad(6));
        end
        % skip it when im_ij contains no object
        if ~any(im_ij(:))
            continue;
        end
        % distance in all three dimensions
        % im_dist_ij = single(bwdist(im_ij) + 1);
        % distance in xy plane
        im_dist_ij = zeros(size(im_ij), 'single');
        for z = 1 : size(im_ij, 3)
            % im_dist_ij(:, :, z) = bwdist(im_ij(:, :, z)) + 1;
            im_dist_ij(:, :, z) = fastPower((bwdist(im_ij(:, :, z)) + 1) / 10, weightDegree);
        end

        s = mbbox(1 : 3) - mbbox_pad(1 : 3) + 1;
        t = mbbox(4 : 6) - mbbox_pad(1 : 3) + 1;
        try 
            indexing3d_mex(im_dist, mbbox, crop3d_mex(im_dist_ij, [s', t']));
        catch ME
            disp(ME)
            im_dist(mbbox(1) : mbbox(4), mbbox(2) : mbbox(5), mbbox(3) : mbbox(6)) = im_dist_ij(s(1) : t(1), s(2) : t(2), s(3) : t(3));
        end
    end
end

% apply a window in z direction
if ismatrix(im_i)
    win_z = 1;
else
    win_z = tukeywin(sz(3) * 1.1, 0.5);
    win_z = win_z(round(sz(3) * 0.05) : round(sz(3) * 0.05) + sz(3) - 1);
end
im_dist = im_dist .* permute(fastPower(win_z, weightDegree), [2, 3, 1]);

im_dist = im_dist .* im_i_orig;
clear im_i_orig im_i;

if dataOrderMat(3) ~= 3
    [~, data_order_reverse_mat] = sort(dataOrderMat);
    im_dist = permute(im_dist, data_order_reverse_mat);
    sz = size(im_dist, 1:3);
end

if ~isempty(distBbox)
    bufferSize = 100;
    dfactor = 0.99;
    winType = 'hann';
    dist_y = distance_weight_single_axis(sz(1), distBbox([1, 4]), bufferSize, dfactor, winType);
    dist_x = distance_weight_single_axis(sz(2), distBbox([2, 5]), bufferSize, dfactor, winType);
    dist_z = distance_weight_single_axis(sz(3), distBbox([3, 6]), bufferSize, dfactor, winType);

    im_dist_wt = dist_y .* permute(dist_x, [2, 1]) .* permute(dist_z, [2, 3, 1]);
    im_dist_wt = fastPower(im_dist_wt, weightDegree);
    if dfactor > 0
        im_dist_wt = max(im_dist_wt, 1e-45);
    end
    im_dist = im_dist .* im_dist_wt;
    clear im_dist_wt;
end

% write to zarr
zarrFilename = bwdistFullpath;
tmpFilename = [zarrFilename '_' uuid];
% write(bim_dist, tmpFilename, "BlockSize", bim_i.BlockSize, "Adapter", ZarrAdapter);
% createZarrFile(tmpFilename, 'chunks', blockSize, 'dtype', 'f4', 'order', 'F', ...
%     'shape', size(im_dist), 'cname', 'zstd', 'level', 2);
if ~exist(tmpFilename, 'dir')
    dimSeparator = '.';
    if prod(ceil(sz ./ blockSize)) > 10000
        dimSeparator = '/';
    end
    createzarr(tmpFilename, dataSize=sz, blockSize=blockSize, dtype='single', dimSeparator=dimSeparator);
end
writezarr(im_dist, tmpFilename);

% mv tmp result folder to output folder
if exist(zarrFilename, 'dir')
    if ~Overwrite
        fprintf('Distance Zarr result %s already exists, skip it!\n', zarrFilename);
        rmdir(tmpFilename, 's');
        return;
    else
        rmdir(zarrFilename, 's');
    end
end
movefile(tmpFilename, zarrFilename);

end

