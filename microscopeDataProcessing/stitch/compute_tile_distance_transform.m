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
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('uuid', '', @ischar);

ip.parse(blockInfoFullname, stitchPath, tileFullpaths, varargin{:});

pr = ip.Results;
blendWeightDegree = pr.blendWeightDegree;
singleDistMap = pr.singleDistMap;
zarrHeaders = pr.zarrHeaders;
blockSize = pr.blockSize;
compressor = pr.compressor;
parseCluster = pr.parseCluster;


% define input, output and function handle for the cluster computing framework
nF = numel(tileFullpaths);
distPath = [stitchPath, 'imdist/'];

% load overlap info
a = load(blockInfoFullname, 'overlap_matrix', 'ol_region_cell');
overlap_matrix = a.overlap_matrix;
ol_region_cell = a.ol_region_cell;
clear a;

% check if the parameter file exist, if so, check if it is the same setting
% as current running. 
paramFullname = [distPath, 'parameters.mat'];
if exist(paramFullname, 'file')
    b = load(paramFullname, 'overlap_matrix', 'ol_region_cell');
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

save('-v7.3', paramFullname, 'ip', 'overlap_matrix', 'ol_region_cell');

if singleDistMap
    tileFullpaths = tileFullpaths(1);
end
inputFullpaths = tileFullpaths;
[~, fsnames] = fileparts(tileFullpaths);
if singleDistMap
    outputFullpaths = {sprintf('%s/%s_wr_%d.zarr', distPath, fsnames, blendWeightDegree)};
    funcStrs = arrayfun(@(x) sprintf('compute_tile_bwdist(''%s'',%d,''%s'',%d,%s,%s,''%s'')', ...
        blockInfoFullname, x, outputFullpaths{x}, blendWeightDegree, string(singleDistMap), ...
        strrep(mat2str(blockSize), ' ', ','), compressor), 1 : 1, 'unif', 0);    
else
    outputFullpaths = cellfun(@(x) sprintf('%s/%s_wr_%d.zarr', distPath, x, blendWeightDegree), fsnames, 'unif', 0);
    funcStrs = arrayfun(@(x) sprintf('compute_tile_bwdist(''%s'',%d,''%s'',%d,%s,%s,''%s'')', ...
        blockInfoFullname, x, outputFullpaths{x}, blendWeightDegree, string(singleDistMap), ...
        strrep(mat2str(blockSize), ' ', ','), compressor), 1 : nF, 'unif', 0);
end

if isempty(zarrHeaders)
    a = load(blockInfoFullname, 'zarrHeaders');
    zarrHeaders = a.zarrHeaders;
end

imSize = zarrHeaders{1}.Size;

% abc cluster
cpusPerTask = min(24, ceil(prod(imSize) * 4 / 1024^3 * 10 / 20));

slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, funcStrs, ...
    'cpusPerTask', cpusPerTask, 'parseCluster', parseCluster);


end

