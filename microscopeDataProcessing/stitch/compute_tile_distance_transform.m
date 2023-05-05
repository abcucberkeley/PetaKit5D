function [outputFullpaths] = compute_tile_distance_transform(blockInfoFullname, stitchPath, tileFullpaths, varargin)
% wrapper for compute distance transform for a tile after removing overlap regions.
% 
% 
% xruan (12/10/2020): check whether the overlap regions are the same if the
% distance files exist. 
% xruan (10/25/2021): add support for a single distance map for all tiles in 
% feather blending (save time for computing).


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInfoFullname', @ischar);
ip.addRequired('stitchPath', @ischar);
ip.addRequired('tileFullpaths', @iscell);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('blendWeightDegree', 10, @isnumeric);
ip.addParameter('singleDistMap', true, @islogical);
ip.addParameter('zarrHeaders', {}, @iscell);
ip.addParameter('blockSize', [500, 500, 500], @isnumeric);
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
zarrHeaders = pr.zarrHeaders;
blockSize = pr.blockSize;
compressor = pr.compressor;
largeZarr =  pr.largeZarr;
poolSize =  pr.poolSize;
parseCluster = pr.parseCluster;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

% define input, output and function handle for the cluster computing framework
nF = numel(tileFullpaths);
distPath = [stitchPath, 'imdist/'];

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
else
    rmdir(distPath, 's');
end
mkdir(distPath);

save('-v7.3', paramFullname, 'ip', 'overlap_matrix', 'ol_region_cell', 'overlap_map_mat');

if singleDistMap
    tileFullpaths = tileFullpaths(1);
end
inputFullpaths = tileFullpaths;
[~, fsnames] = fileparts(tileFullpaths);

if largeZarr
    switch numel(poolSize)
        case 3
            poolSize_str = sprintf('pooling_%d_%d_%d', poolSize(1), poolSize(2), poolSize(3));
        case 6
            poolSize_str = sprintf('pooling_%d_%d_%d_%d_%d_%d', poolSize(1), ...
                poolSize(2), poolSize(3), poolSize(4), poolSize(5), poolSize(6));
    end
    if singleDistMap
        outputFullpaths = {sprintf('%s/%s_z_%s_wr_%d.zarr', distPath, fsnames, ...
            poolSize_str, blendWeightDegree)};
        funcStrs = arrayfun(@(x) sprintf('compute_tile_bwdist_mip_slabs(''%s'',%d,''%s'',%d,%s,%s,''%s'',%s,%s)', ...
            blockInfoFullname, x, outputFullpaths{x}, blendWeightDegree, string(singleDistMap), ...
            strrep(mat2str(blockSize), ' ', ','), compressor, strrep(mat2str(poolSize), ' ', ','), ...
            string(Overwrite)), 1 : 1, 'unif', 0);    
    else
        outputFullpaths = cellfun(@(x) sprintf('%s/%s_z_%s_wr_%d.zarr', distPath, x, ...
            poolSize_str, blendWeightDegree), fsnames, 'unif', 0);
        funcStrs = arrayfun(@(x) sprintf('compute_tile_bwdist_mip_slabs(''%s'',%d,''%s'',%d,%s,%s,''%s'',%s,%s)', ...
            blockInfoFullname, x, outputFullpaths{x}, blendWeightDegree, string(singleDistMap), ...
            strrep(mat2str(blockSize), ' ', ','), compressor, strrep(mat2str(poolSize), ' ', ','), ...
            string(Overwrite)), 1 : nF, 'unif', 0);
    end
else
    if singleDistMap
        outputFullpaths = {sprintf('%s/%s_wr_%d.zarr', distPath, fsnames, blendWeightDegree)};
        funcStrs = arrayfun(@(x) sprintf('compute_tile_bwdist(''%s'',%d,''%s'',%d,%s,%s,''%s'',%s)', ...
            blockInfoFullname, x, outputFullpaths{x}, blendWeightDegree, string(singleDistMap), ...
            strrep(mat2str(blockSize), ' ', ','), compressor, string(Overwrite)), 1 : 1, 'unif', 0);    
    else
        outputFullpaths = cellfun(@(x) sprintf('%s/%s_wr_%d.zarr', distPath, x, blendWeightDegree), fsnames, 'unif', 0);
        funcStrs = arrayfun(@(x) sprintf('compute_tile_bwdist(''%s'',%d,''%s'',%d,%s,%s,''%s'',%s)', ...
            blockInfoFullname, x, outputFullpaths{x}, blendWeightDegree, string(singleDistMap), ...
            strrep(mat2str(blockSize), ' ', ','), compressor,string(Overwrite)), 1 : nF, 'unif', 0);
    end
end

if isempty(zarrHeaders)
    a = load(blockInfoFullname, 'tileFns');
    tileFns = a.tileFns;
end

imSize = getImageSize(tileFns{1});

% memory needed
memAllocate = prod(imSize) * 4 / 1024^3 * 10;

generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, funcStrs, ...
    'memAllocate', memAllocate, 'parseCluster', parseCluster, mccMode=mccMode, ...
    ConfigFile=ConfigFile);


end

