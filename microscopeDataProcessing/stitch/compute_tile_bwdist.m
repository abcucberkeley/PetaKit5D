function [] = compute_tile_bwdist(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, singleDistMap, blockSize, compressor, Overwrite)
% compute distance transform for a tile after removing overlap regions.
% 
% Author: Xiongtao Ruan (10/29/2020)
% xruan (11/02/2020): only compute distance in overlap regions, and also
% only consider distance in xy plane.
% xruan (11/19/2022): compute feather power within the distance map to save
% time for stitching processing

if nargin < 8
    Overwrite = false;
end

uuid = get_uuid();

if exist(bwdistFullpath, 'dir') 
    return;
end

persistent zarrHeaders overlap_matrix ol_region_cell bkFn datenum
dir_info = dir(blockInfoFullname);
if ~(strcmpi(bkFn, blockInfoFullname) && datenum == dir_info.datenum) || isempty(zarrHeaders)
    bkFn = blockInfoFullname;
    a = load(blockInfoFullname,'zarrHeaders', 'overlap_matrix', 'ol_region_cell');
    zarrHeaders = a.zarrHeaders;
    overlap_matrix = a.overlap_matrix;
    ol_region_cell = a.ol_region_cell;
    dir_info = dir(bkFn);
    datenum = dir_info.datenum;
    clear a;
end

nF = numel(zarrHeaders);
i = tileInd;
bim_i = zarrHeaders{i};
% [~, fsname] = fileparts(bim_i.Source);

% im_i = gather(bim_i);
im_i = bim_i.Adapter.getIORegion([1, 1, 1], bim_i.Size);
im_i_orig = im_i ~= 0;
im_i([1, end], :, :) = 0;
im_i(:, [1, end], :) = 0;
if ndims(im_i) == 3
    im_i(:, :, [1, end]) = 0;
end
im_i = im_i == 0;

counter = 1;
sz = size(im_i, [1 : 3]);
im_dist = ones(sz, 'single');
if singleDistMap
    im_i_c = im_i(:, :, round((sz(3) + 1) / 2));
    im_dist_c = ((bwdist(im_i_c) + 1) / 10) .^ weightDegree;
    for z = 1 : sz(3)
        im_i_z = im_i(:, :, z);
        if sum(im_i_z ~= im_i_c, 'all') == 0
            im_dist(:, :, z) = im_dist_c;
        else
            im_dist(:, :, z) = ((bwdist(im_i(:, :, z)) + 1) / 10) .^ weightDegree;
        end
    end
else
    for j = 1 : nF
        if ~overlap_matrix(i, j) && ~overlap_matrix(j, i)
            continue;
        end
        if i < j
            mregion = ol_region_cell{i, j}{1};
        else
            mregion = ol_region_cell{j, i}{2};
        end

        midx = mregion;
        % dsr_ol = dsr(mregion_f(2, 1) : mregion_f(2, 2), mregion_f(1, 1) : mregion_f(1, 2), mregion_f(3, 1) : mregion_f(3, 2));
        % dsr(midx(2, 1) : midx(2, 2), midx(1, 1) : midx(1, 2), midx(3, 1) : midx(3, 2)) = 0;

        mbbox = [midx(2, 1), midx(1, 1), midx(3, 1), midx(2, 2), midx(1, 2), midx(3, 2)];
        if any(mbbox(4 : 6) < mbbox(1 : 3))
            continue;
        end

        mbbox_pad = [max(mbbox(1 : 3) - 1, 1), min(mbbox(4 : 6) + 1, sz)];
        im_ij = im_i(mbbox_pad(1) : mbbox_pad(4), mbbox_pad(2) : mbbox_pad(5), mbbox_pad(3) : mbbox_pad(6));
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
            im_dist_ij(:, :, z) = ((bwdist(im_ij(:, :, z)) + 1) / 10) .^ weightDegree;
        end

        s = mbbox(1 : 3) - mbbox_pad(1 : 3) + 1;
        t = mbbox(4 : 6) - mbbox_pad(1 : 3) + 1;
        im_dist(mbbox(1) : mbbox(4), mbbox(2) : mbbox(5), mbbox(3) : mbbox(6)) = im_dist_ij(s(1) : t(1), s(2) : t(2), s(3) : t(3));
        % im_i(mbbox(1) : mbbox(4), mbbox(2) : mbbox(5), mbbox(3) : mbbox(6)) = false;

        counter = counter + 1;
    end
end

% apply a window in z direction
if ismatrix(im_i)
    win_z = 1;
else
    win_z = tukeywin(sz(3) * 1.1, 0.5);
    win_z = win_z(round(sz(3) * 0.05) : round(sz(3) * 0.05) + sz(3) - 1);
end
im_dist = im_dist .* (permute(win_z, [2, 3, 1]) .^ weightDegree);

im_dist = im_dist .* im_i_orig;
clear im_i_orig im_i;

% write to zarr
zarrFilename = bwdistFullpath;
tmpFilename = [zarrFilename '_' uuid];
% write(bim_dist, tmpFilename, "BlockSize", bim_i.BlockSize, "Adapter", ZarrAdapter);
% createZarrFile(tmpFilename, 'chunks', blockSize, 'dtype', 'f4', 'order', 'F', ...
%     'shape', size(im_dist), 'cname', 'zstd', 'level', 2);
try
    createZarrFile(tmpFilename, 'chunks', blockSize, 'dtype', 'f4', 'order', 'F', ...
        'shape', size(im_dist), 'cname', compressor, 'level', 1);    
    % bim = blockedImage(tmpFilename, sz, blockSize, zeros(1, 'single'), "Adapter", CZarrAdapter, 'mode', 'w');
catch ME
    disp(ME);
    bim = blockedImage(tmpFilename, sz, blockSize, zeros(1, 'single'), "Adapter", ZarrAdapter, 'mode', 'w');
    bim.Adapter.close();    
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
