function [] = compute_tile_bwdist_mip_slabs(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, singleDistMap, blockSize, shardSize, compressor, poolSize, distBbox, Overwrite)
% compute distance transform for a tile for large zarr file with MIP slab
% (only z for now). 
%
% copied from compute_tile_bwdist.m
% 
% Author: Xiongtao Ruan (05/02/2023)

if nargin < 10
    Overwrite = false;
end

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

tileFn = tileFns{i};
[dataPath, fsn] = fileparts(tileFn);
switch numel(poolSize) 
    case 3
        MIPZFn = sprintf('%s/MIP_Slabs_pooling_%d_%d_%d/%s_z.zarr', dataPath, poolSize(1), poolSize(2), poolSize(3), fsn);
        poolSize = [poolSize, 1, 1, 1];
    case {6, 9}
        MIPZFn = sprintf('%s/MIP_Slabs_pooling_%d_%d_%d_%d_%d_%d/%s_z.zarr', dataPath, ...
            poolSize(1), poolSize(2), poolSize(3), poolSize(4), poolSize(5), poolSize(6), fsn);        
end

im_i = readzarr(MIPZFn);
im_i_orig = im_i ~= 0;
im_i([1, end], :, :) = 0;
im_i(:, [1, end], :) = 0;
if ndims(im_i) == 3
    im_i(:, :, [1, end]) = 0;
end
im_i = im_i == 0;

sz = size(im_i, [1 : 3]);
im_dist = ones(sz, 'single');

% distance in xy
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

if ~isempty(distBbox)
    distBbox = distBbox ./ [poolSize(4 : 5), poolSize(3), poolSize(4 : 5), poolSize(3)];
    distBbox = max(1, round(distBbox));
    
    bufferSize = 100;
    % bufferSize = 50;
    bufferSizes = max(1, round(bufferSize ./ [poolSize(4 : 5), poolSize(3)]));
    dfactor = 0.99;
    winType = 'hann';
    dist_y = distance_weight_single_axis(sz(1), distBbox([1, 4]), bufferSizes(1), dfactor, winType);
    dist_x = distance_weight_single_axis(sz(2), distBbox([2, 5]), bufferSizes(2), dfactor, winType);
    dist_z = distance_weight_single_axis(sz(3), distBbox([3, 6]), bufferSizes(3), dfactor, winType);

    im_dist_wt = (dist_y .^ weightDegree) .* permute(dist_x .^ weightDegree, [2, 1]) .* permute(dist_z .^ weightDegree, [2, 3, 1]);
    if dfactor > 0
        im_dist_wt = max(im_dist_wt, 1e-40);
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
try
    zarrSubSize = [20, 20, 20];
    createzarr(tmpFilename, dataSize=size(im_dist), blockSize=blockSize, shardSize=shardSize, ...
        dtype='single', compressor=compressor, zarrSubSize=zarrSubSize);
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
