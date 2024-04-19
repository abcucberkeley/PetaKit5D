function [block_info_fullname, PerBlockInfoFullpaths] = stitch_process_block_info(int_xyz_shift, imSizes, nvSize, blockSize, overlap_matrix, ol_region_cell, half_ol_region_cell, overlap_map_mat, BorderSize, tileFns, stichInfoPath, nv_fsname, isPrimaryCh, varargin)
% move the block information processing code to this file to enable cluster
% processing for large scale data with many tiles
%  
%
% Author: Xiongtao Ruan (04/10/2023)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('int_xyz_shift', @isnumeric);
ip.addRequired('imSizes', @isnumeric);
ip.addRequired('nvSize', @isnumeric);
ip.addRequired('blockSize', @isnumeric);
ip.addRequired('overlap_matrix', @islogical);
ip.addRequired('ol_region_cell', @iscell);
ip.addRequired('half_ol_region_cell', @iscell);
ip.addRequired('overlap_map_mat', @isnumeric);
ip.addRequired('BorderSize', @isnumeric);
ip.addRequired('tileFns', @iscell);
ip.addRequired('stichInfoPath', @ischar);
ip.addRequired('nv_fsname', @ischar);
ip.addRequired('isPrimaryCh', @islogical);
ip.addParameter('stitchInfoFullpath', '', @ischar);
ip.addParameter('stitch2D', false, @islogical);
ip.addParameter('taskSize', 10, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(int_xyz_shift, imSizes, nvSize, blockSize, overlap_matrix, ol_region_cell, half_ol_region_cell, overlap_map_mat, BorderSize, tileFns, stichInfoPath, nv_fsname, isPrimaryCh, varargin{:});

pr = ip.Results;
stitchInfoFullpath = pr.stitchInfoFullpath;
stitch2D = pr.stitch2D;
taskSize = pr.taskSize;
uuid = pr.uuid;
parseCluster = pr.parseCluster;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

% t0 = tic();
fprintf('Process stitch block information...\n')

bSubSz = ceil(nvSize ./ blockSize);
numBlocks = prod(bSubSz);
nF = size(imSizes, 1);
nxs = nvSize(2);
nys = nvSize(1);
nzs = nvSize(3);

bboxes = zeros(nF, 6);
bboxStart_mat = zeros(nF, 3);
bboxEnd_mat = zeros(nF, 3);

for i = 1 : nF
    st_idx = int_xyz_shift(i, :);
    % if any(st_idx < 1) || any(st_idx > [nxs, nys, nzs])
    % bound idx to the positions of the stitched image
    sx = imSizes(i, 2);
    sy = imSizes(i, 1);
    sz = imSizes(i, 3);
    xridx = max(1, 1 - st_idx(1)) : min(sx, nxs - st_idx(1));
    yridx = max(1, 1 - st_idx(2)) : min(sy, nys - st_idx(2));
    zridx = max(1, 1 - st_idx(3)) : min(sz, nzs - st_idx(3));

    xidx = st_idx(1) + xridx;
    yidx = st_idx(2) + yridx;
    zidx = st_idx(3) + zridx;
    if stitch2D
        zridx = 1;
        zidx = 1;
    end    
    if isempty(xidx) || isempty(yidx) || isempty(zidx)
        continue;
    end
        
    bboxes(i, :) = [yidx(1), xidx(1), zidx(1), yidx(end), xidx(end), zidx(end)];
    bboxStart_mat(i, :) = [yridx(1), xridx(1), zridx(1)];
    bboxEnd_mat(i, :) = [yridx(end), xridx(end), zridx(end)];
end

block_info_tmp_fullname = sprintf('%s/%s_block_info_task_size_%d_%s.mat', stichInfoPath, nv_fsname, taskSize, uuid);
block_info_fullname = sprintf('%s/%s_block_info_task_size_%d.mat', stichInfoPath, nv_fsname, taskSize);

save('-v7.3', block_info_tmp_fullname, 'tileFns', 'int_xyz_shift', 'imSizes', 'overlap_matrix', 'half_ol_region_cell', ...
    'ol_region_cell', 'overlap_map_mat', 'bboxes', 'bboxStart_mat', 'bboxEnd_mat', 'nvSize', 'blockSize', 'bSubSz', 'BorderSize');
movefile(block_info_tmp_fullname, block_info_fullname);

% organize the block info w.r.t the blocks and directly save them to json
% files with cluster computing

if isPrimaryCh 
    PerBlockInfoPath = sprintf('%s/%s/', stichInfoPath, nv_fsname);
    mkdir(PerBlockInfoPath);
    PerBlockInfoFlagPath = sprintf('%s/block_flags/', PerBlockInfoPath);
    mkdir(PerBlockInfoFlagPath);    
else
    [pstr, fsn] = fileparts(stitchInfoFullpath);
    PerBlockInfoPath = [pstr, '/', fsn];
end

numTasks = ceil(numBlocks / taskSize);
PerBlockInfoFullpaths = cell(numTasks, 1);
for t = 1 : numTasks
    bs = (t - 1) * taskSize + 1; 
    bt = min(t * taskSize, numBlocks);
    
    PerBlockInfoFullpaths{t} = sprintf('%s/stitch_block_info_blocks_%d_%d.json', PerBlockInfoPath, bs, bt);
end

if ~isPrimaryCh
    return;
end

blk_taskSize = max(1e5, min(1e6, taskSize * 500));
blk_taskSize = ceil(blk_taskSize / taskSize) * taskSize;
numTasks = ceil(numBlocks / blk_taskSize);
FlagFullpaths = cell(numTasks, 1);
funcStrs = cell(numTasks, 1);
for t = 1 : numTasks
    bs = (t - 1) * blk_taskSize + 1; 
    bt = min(t * blk_taskSize, numBlocks);

    % save block info searately for each task for faster access
    FlagFullpath = sprintf('%s/stitch_block_info_task_size_%d_blocks_%d_%d.mat', PerBlockInfoFlagPath, taskSize, bs, bt);
    FlagFullpaths{t} = FlagFullpath;
    funcStrs{t} = sprintf(['stitch_organize_block_info([%s],%d,''%s'',''%s'',', ...
        '''%s'',''BorderSize'',[%s])'], strrep(mat2str([bs, bt]), ' ', ','), taskSize, ...
        block_info_fullname, PerBlockInfoPath, FlagFullpath, strrep(mat2str(BorderSize), ' ', ','));
end

if ~isPrimaryCh
    return;
end

% cluster setting
cpusPerTask = 2;
masterCompute = true;
memAllocate = blk_taskSize * 0.0001;
maxTrialNum = 2;

if numTasks <= 2
    parseCluster = false;
end
minTaskJobNum = 1;
if numTasks >= 4
    minTaskJobNum = 2;
end

inputFullpaths = repmat({block_info_fullname}, numTasks, 1);
outputFullpaths = FlagFullpaths;

for i = 1 : 3
    is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', cpusPerTask, 'memAllocate', memAllocate, 'minTaskJobNum', minTaskJobNum, ...
        'maxTrialNum', maxTrialNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster, ...
        'mccMode', mccMode, 'ConfigFile', ConfigFile);
    if all(is_done_flag)
        break;
    end
end

% fprintf('Done!\n')
% toc(t0);

end

