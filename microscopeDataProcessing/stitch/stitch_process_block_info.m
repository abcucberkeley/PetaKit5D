function [block_info_fullname, PerBlockInfoFullpaths, block_info_bytes] = stitch_process_block_info(int_xyz_shift, imSizes, nvSize, blockSize, overlap_matrix, ol_region_cell, half_ol_region_cell, overlap_map_mat, BorderSize, tileFns, stichInfoPath, nv_fsname, isPrimaryCh, varargin)
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
ip.addParameter('maxFileNumPerFolder', 20000, @(x)isscalar(x));  
ip.addParameter('uuid', '', @ischar);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(int_xyz_shift, imSizes, nvSize, blockSize, overlap_matrix, ol_region_cell, half_ol_region_cell, overlap_map_mat, BorderSize, tileFns, stichInfoPath, nv_fsname, isPrimaryCh, varargin{:});

pr = ip.Results;
stitchInfoFullpath = pr.stitchInfoFullpath;
stitch2D = pr.stitch2D;
taskSize = pr.taskSize;
maxFileNumPerFolder =  pr.maxFileNumPerFolder;
uuid = pr.uuid;
parseCluster = pr.parseCluster;
mccMode = pr.mccMode;
configFile = pr.configFile;

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
    PerBlockInfoPath = sprintf('%s/%s_task_size_%d/', stichInfoPath, nv_fsname, taskSize);
    mkdir(PerBlockInfoPath);
    PerBlockInfoFlagPath = sprintf('%s/block_flags/', PerBlockInfoPath);
    mkdir(PerBlockInfoFlagPath);    
else
    [pstr, fsn] = fileparts(stitchInfoFullpath);
    PerBlockInfoPath = sprintf('%s/%s_task_size_%d/', pstr, fsn, taskSize);
end

numTasks = ceil(numBlocks / taskSize);
% if the number of task greater than 20000, save them in separate folders
nFile = maxFileNumPerFolder;
PerBlockInfoPath_cell = {PerBlockInfoPath};
if numTasks > nFile
    nFolder = ceil(numTasks / nFile);
    PerBlockInfoPath_cell = cell(nFolder, 1);
    for t = 1 : nFolder
        bs = (t - 1) * nFile + 1;
        bt = min(t * nFile, numTasks);
        PerBlockInfoPath_t = sprintf('%s/task_size_%d_tasks_%d_%d/', PerBlockInfoPath, taskSize, bs, bt);
        mkdir(PerBlockInfoPath_t);
        PerBlockInfoPath_cell{t} = PerBlockInfoPath_t;
    end
end
PerBlockInfoFullpaths = cell(numTasks, 1);
for t = 1 : numTasks
    bs = (t - 1) * taskSize + 1;
    bt = min(t * taskSize, numBlocks);
    
    if numTasks > nFile
        PerBlockInfoFullpaths{t} = sprintf('%s/stitch_block_info_blocks_%d_%d.json', PerBlockInfoPath_cell{ceil(t / nFile)}, bs, bt);        
    else
        PerBlockInfoFullpaths{t} = sprintf('%s/stitch_block_info_blocks_%d_%d.json', PerBlockInfoPath, bs, bt);        
    end
end

block_info_bytes = [];
if ~isPrimaryCh
    % collect the block info data size if using cluster
    if parseCluster
        block_info_bytes = get_block_info_file_bytes(PerBlockInfoPath_cell);
    end    
    return;
end

blk_taskSize = max(1e5, min(1e6, taskSize * 500));
blk_taskSize = ceil(blk_taskSize / taskSize) * taskSize;
if numTasks > nFile
    blk_taskSize = min(blk_taskSize, nFile * taskSize);
    if blk_taskSize < nFile * taskSize
        for i = 1 : ceil(log2(nFile * taskSize / blk_taskSize))
            cur_blk_taskSize = nFile * taskSize / 2^i;
            if blk_taskSize > cur_blk_taskSize || floor(cur_blk_taskSize) ~= cur_blk_taskSize
                break;
            end
        end
        blk_taskSize = nFile * taskSize / 2^(i-1);
    end
end
blk_numTasks = ceil(numBlocks / blk_taskSize);
FlagFullpaths = cell(blk_numTasks, 1);
funcStrs = cell(blk_numTasks, 1);
for t = 1 : blk_numTasks
    bs = (t - 1) * blk_taskSize + 1; 
    bt = min(t * blk_taskSize, numBlocks);

    % save block info searately for each task for faster access
    FlagFullpath = sprintf('%s/stitch_block_info_task_size_%d_blocks_%d_%d.mat', PerBlockInfoFlagPath, taskSize, bs, bt);
    FlagFullpaths{t} = FlagFullpath;
    PerBlockInfoPath_t = PerBlockInfoPath;
    if numTasks > nFile
        PerBlockInfoPath_t = PerBlockInfoPath_cell{ceil(bs / (nFile * taskSize))};
    end
    funcStrs{t} = sprintf(['stitch_organize_block_info([%s],%d,''%s'',''%s'',', ...
        '''%s'',''BorderSize'',[%s])'], strrep(mat2str([bs, bt]), ' ', ','), taskSize, ...
        block_info_fullname, PerBlockInfoPath_t, FlagFullpath, strrep(mat2str(BorderSize), ' ', ','));
end

% cluster setting
cpusPerTask = 2;
masterCompute = true;
memAllocate = blk_taskSize * 0.0001;
maxTrialNum = 2;

minTaskJobNum = 1;
if blk_numTasks >= 4
    minTaskJobNum = 2;
end

inputFullpaths = repmat({block_info_fullname}, blk_numTasks, 1);
outputFullpaths = FlagFullpaths;

for i = 1 : 3
    is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', cpusPerTask, 'memAllocate', memAllocate, 'minTaskJobNum', minTaskJobNum, ...
        'maxTrialNum', maxTrialNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster, ...
        'mccMode', mccMode, 'configFile', configFile);
    if all(is_done_flag)
        break;
    end
end

if ~all(is_done_flag)
    error('Some block (%d / %d) info file cannot finished!', sum(~is_done_flag), numel(is_done_flag));
end

% collect the block info data size if using cluster
if parseCluster
    block_info_bytes = get_block_info_file_bytes(PerBlockInfoPath_cell);
end

end


function [block_info_bytes] = get_block_info_file_bytes(PerBlockInfoPath_cell)

nF = numel(PerBlockInfoPath_cell);
s_mat_cell = cell(nF, 1);
block_info_bytes_cell = cell(nF, 1);
for t = 1 : numel(PerBlockInfoPath_cell)
    PerBlockInfoPath = PerBlockInfoPath_cell{t};
    dir_info = dir([PerBlockInfoPath, '*.json']);
    fsn = {dir_info.name};
    s_mat = regexp(fsn, 'stitch_block_info_blocks_(\d+)_\d+.json', 'tokens');
    s_mat = cellfun(@(x) str2double(x{1}{1}), s_mat);
    block_info_bytes = [dir_info.bytes];

    s_mat_cell{t} = s_mat(:);
    block_info_bytes_cell{t} = block_info_bytes(:);
end

s_mat = cat(1, s_mat_cell{:});
block_info_bytes = cat(1, block_info_bytes_cell{:});
[~, inds] = sort(s_mat);
block_info_bytes = block_info_bytes(inds);

end

