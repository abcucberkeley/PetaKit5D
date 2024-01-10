function [outputFullpaths, imdistFileIdx] = compute_tile_distance_transform(blockInfoFullname, stitchPath, tileFullpaths, varargin)
% wrapper for compute distance transform for a tile after removing overlap regions.
% 
% 
% xruan (12/10/2020): check whether the overlap regions are the same if the
% distance files exist. 
% xruan (10/25/2021): add support for a single distance map for all tiles in 
% feather blending (save time for computing).
% xruan (05/26/2023): add support for a single distance map for each location 
% xruan (08/20/2023): add support for bounding boxes for distance transform. 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInfoFullname', @ischar);
ip.addRequired('stitchPath', @ischar);
ip.addRequired('tileFullpaths', @iscell);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('blendWeightDegree', 10, @isnumeric);
ip.addParameter('singleDistMap', true, @islogical);
ip.addParameter('locIds', [], @isnumeric);
ip.addParameter('distBboxes', [], @isnumeric);
ip.addParameter('blockSize', [500, 500, 500], @isnumeric);
ip.addParameter('shardSize', [], @isnumeric);
ip.addParameter('compressor', 'lz4', @ischar);
ip.addParameter('largeZarr', false, @islogical);
ip.addParameter('poolSize', [], @isnumeric);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);
ip.addParameter('uuid', '', @ischar);

ip.parse(blockInfoFullname, stitchPath, tileFullpaths, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
blendWeightDegree = pr.blendWeightDegree;
singleDistMap = pr.singleDistMap;
locIds = pr.locIds;
distBboxes = pr.distBboxes;
blockSize = pr.blockSize;
shardSize = pr.shardSize;
compressor = pr.compressor;
largeZarr =  pr.largeZarr;
poolSize =  pr.poolSize;
parseCluster = pr.parseCluster;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;
uuid = pr.uuid;
if isempty(uuid)
    uuid = get_uuid();
end

% define input, output and function handle for the cluster computing framework
nF = numel(tileFullpaths);
[~, fsn] = fileparts(tileFullpaths{1});
if ispc
    distPath = sprintf('%s/imdist/', stitchPath);
else
    distPath = sprintf('%s/imdist/%s/', stitchPath, fsn);
end
if ~exist(distPath, 'dir')
    mkdir_recursive(distPath);
end

% load overlap info
a = load(blockInfoFullname, 'overlap_matrix', 'ol_region_cell', 'overlap_map_mat');
overlap_matrix = a.overlap_matrix;
ol_region_cell = a.ol_region_cell;
overlap_map_mat = a.overlap_map_mat;
clear a;

% check if the parameter file exist, if so, check if it is the same setting
% as current running. 
paramFullname = [distPath, 'parameters.mat'];
if exist(paramFullname, 'file')
    b = load(paramFullname, 'overlap_matrix', 'ol_region_cell', 'overlap_map_mat');
    isSameSetting = singleDistMap || (all(size(overlap_matrix) == size(b.overlap_matrix)) && ...
        all(overlap_matrix == b.overlap_matrix, 'all') && ...
        all(cellfun(@isequal, ol_region_cell, b.ol_region_cell), 'all'));
    if ~isSameSetting
        disp('The existing distance transforms are not in the same setting as current running, delete them!')
        rmdir(distPath, 's');
    end
end
mkdir(distPath);

paramTmpName = sprintf('%s/parameters_%s.mat', distPath, uuid);
save('-v7.3', paramTmpName, 'ip', 'overlap_matrix', 'ol_region_cell', 'overlap_map_mat');
movefile(paramTmpName, paramFullname);

imdistFileIdx = (1 : nF)';
distTileIdx = (1 : nF)';
if singleDistMap
    [uniq_locIds, uniq_inds] = unique(locIds);    
    tileFullpaths = tileFullpaths(uniq_inds);
    imdistFileIdx = zeros(nF, 1);
    imSizes = zeros(numel(uniq_locIds), 3);
    for i = 1 : numel(uniq_locIds)
        imdistFileIdx(locIds == uniq_locIds(i)) = i;
        imSizes(i, :) = getImageSize(tileFullpaths{i});
    end
    distTileIdx = uniq_inds;

    % merge different locations with same image size
    if numel(uniq_locIds) > 1
        mapped_inds = 1 : numel(uniq_locIds);
        for i = 2 : numel(uniq_locIds)
            if any(all(imSizes(i, :) == imSizes(1 : i - 1, :), 2))
                mapped_inds(i) = find(all(imSizes(i, :) == imSizes(1 : i - 1, :), 2), 1, 'first');
            end
        end
        uniq_mapped_inds = unique(mapped_inds);
        tileFullpaths = tileFullpaths(uniq_mapped_inds);
        for i = 1 : numel(uniq_locIds)
            imdistFileIdx(locIds == uniq_locIds(i)) = find(mapped_inds(i) == uniq_mapped_inds);
        end
        distTileIdx = uniq_inds(uniq_mapped_inds);
    end
    if ~isempty(distBboxes)
        if size(distBboxes, 1) == 1
            if numel(distTileIdx) > 1
                distBboxes = repmat(distBboxes, numel(distTileIdx), 1);
            end
        elseif size(distBboxes, 1) ~= numel(distTileIdx)
            error('Number of distance bounding box does not match the number of unique tile image sizes!')
        end
    end    
end
nF = numel(tileFullpaths);
inputFullpaths = tileFullpaths;
[~, fsnames] = fileparts(tileFullpaths);
if nF == 1
    fsnames = {fsnames};
end

if isempty(distBboxes)
    distBbox_strs = repmat({'[]'}, nF, 1);
else
    distBbox_strs = arrayfun(@(x) strrep(mat2str(distBboxes(x, :)), ' ', ','), 1 : nF, 'unif', 0);
end
    
if largeZarr
    switch numel(poolSize)
        case 3
            poolSize_str = sprintf('pooling_%d_%d_%d', poolSize(1), poolSize(2), poolSize(3));
        case {6, 9}
            poolSize_str = sprintf('pooling_%d_%d_%d_%d_%d_%d', poolSize(1), ...
                poolSize(2), poolSize(3), poolSize(4), poolSize(5), poolSize(6));
    end
    outputFullpaths = cellfun(@(x) sprintf('%s/%s_z_%s_wr_%d.zarr', distPath, x, ...
        poolSize_str, blendWeightDegree), fsnames, 'unif', 0);
    funcStrs = arrayfun(@(x) sprintf('compute_tile_bwdist_mip_slabs(''%s'',%d,''%s'',%d,%s,%s,%s,''%s'',%s,%s,%s)', ...
        blockInfoFullname, distTileIdx(x), outputFullpaths{x}, blendWeightDegree, ...
        string(singleDistMap), strrep(mat2str(blockSize), ' ', ','), strrep(mat2str(shardSize), ' ', ','), ...
        compressor, strrep(mat2str(poolSize), ' ', ','), distBbox_strs{x}, string(Overwrite)), 1 : nF, 'unif', 0);
else
    outputFullpaths = cellfun(@(x) sprintf('%s/%s_wr_%d.zarr', distPath, x, blendWeightDegree), fsnames, 'unif', 0);
    funcStrs = arrayfun(@(x) sprintf('compute_tile_bwdist(''%s'',%d,''%s'',%d,%s,%s,%s,''%s'',%s,%s)', ...
        blockInfoFullname, distTileIdx(x), outputFullpaths{x}, blendWeightDegree, ...
        string(singleDistMap), strrep(mat2str(blockSize), ' ', ','), strrep(mat2str(shardSize), ' ', ','), ...
        compressor, distBbox_strs{x}, string(Overwrite)), 1 : nF, 'unif', 0);
end

imSize = getImageSize(tileFullpaths{1});

% memory needed
memAllocate = prod(imSize) * 4 / 1024^3 * 10;

minTaskJobNum = 1;
if memAllocate > 10
    minTaskJobNum = numel(funcStrs);
end

generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, funcStrs, ...
    'memAllocate', memAllocate, 'minTaskJobNum', minTaskJobNum, 'parseCluster', parseCluster, ...
    mccMode=mccMode, ConfigFile=ConfigFile);

end

