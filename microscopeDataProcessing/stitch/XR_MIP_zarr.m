function XR_MIP_zarr(zarrFullpath, varargin)
% save MIP for large scale zarr file. The idea is to first generate MIPs
% for each batch for all three axis, and then generate final MIPs. 
% 
% 
% Author: Xiongtao Ruan (02/17/2022)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addParameter('mipDirStr', 'MIPs', @ischar); % y, x, z
ip.addParameter('axis', [1, 1, 1], @isnumeric); % y, x, z
ip.addParameter('BatchSize', [2048, 2048, 2048] , @isnumeric); % in y, x, z
ip.addParameter('poolSize', [] , @isnumeric); % pooling size for mips
ip.addParameter('zarrSubSize', [20, 20, 20] , @isnumeric); % in y, x, z
ip.addParameter('mipSlab', false, @islogical); % compute MIP slabs (without the final MIPs)
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('jobLogDir', '../job_logs/', @ischar);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(zarrFullpath, varargin{:});

pr = ip.Results;
mipDirStr = pr.mipDirStr;
axis = pr.axis;
BatchSize = pr.BatchSize;
poolSize = pr.poolSize;
zarrSubSize = pr.zarrSubSize;
mipSlab = pr.mipSlab;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;

zarrFullpath = strip(zarrFullpath, 'right', filesep);
[dataPath, fsname, ext] = fileparts(zarrFullpath);

MIPPath = sprintf('%s/%s/', dataPath, mipDirStr);
if ~exist(MIPPath, 'dir')
    mkdir(MIPPath);
end

axis_strs = {'y', 'x', 'z'};
MIPFullpaths = cellfun(@(x) sprintf('%s/%s_MIP_%s.tif', MIPPath, fsname, x), axis_strs, 'unif', 0);
MIPZarrpaths = cellfun(@(x) sprintf('%s/%s_%s.zarr', MIPPath, fsname, x), axis_strs, 'unif', 0);

done_flag = false(3, 1);
for i = 1 : 3
    done_flag(i) = (axis(i) == 0) || ((~mipSlab && exist(MIPFullpaths{i}, 'file')) || (mipSlab && exist(MIPZarrpaths{i}, 'file')));
end
if all(done_flag)
    disp('The output results exist, skip it!');
    return;
end

fprintf('Start Large-file MIP for %s...\n', fsname);

tic
zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', MIPPath, fsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end 

try
    bim = blockedImage(zarrFullpath, 'Adapter', CZarrAdapter);
catch ME
    disp(ME);
    bim = blockedImage(zarrFullpath, 'Adapter', ZarrAdapter);    
end
imSize = bim.Size;
dtype = bim.ClassUnderlying;

% pool size for other axes
poolSize_1 = [1, 1, 1];
if isempty(poolSize)
    inBlockSize = bim.BlockSize;
else
    if numel(poolSize) == 6
        poolSize_1 = poolSize(4 : 6);
        poolSize = poolSize(1 : 3);        
    end
    inBlockSize = lcm(poolSize, poolSize_1);
end

% MIPs for each block
BatchSize = min(imSize, max(BatchSize, ceil(BatchSize ./ inBlockSize) .* inBlockSize));
bSubSz = ceil(imSize ./ BatchSize);
numBatch = prod(bSubSz);

% in case of out size too large for very large data that causes oom for the main 
% job, increase batch size if necessary.
byteNum = dataTypeToByteNumber(dtype);
if ~mipSlab
    outVolSizes = zeros(3, 1);
    for i = 1 : 3
        outVolSizes(i) = prod([imSize(setdiff(1 : 3, i)), bSubSz(i)]) * byteNum / 1024^3;
    end
    
    % if the max intermediate MIP files is greater than 100 GB, increase the BatchSize
    % by the blockSize in the axis with largest outVolSizes
    while any(outVolSizes > 100) && prod(BatchSize) * byteNum / 1024^3 < 100
        [~, ind] = max(outVolSizes);
        BatchSize = min(imSize, BatchSize + inBlockSize .* ((1 : 3) == ind));
        bSubSz = ceil(imSize ./ BatchSize);
        numBatch = prod(bSubSz);
    
        for i = 1 : 3
            outVolSizes(i) = prod([imSize(setdiff(1 : 3, i)), bSubSz(i)]) * byteNum / 1024^3;
        end
    end
end

[Y, X, Z] = ndgrid(1 : bSubSz(1), 1 : bSubSz(2), 1 : bSubSz(3));
bSubs = [Y(:), X(:), Z(:)];
clear Y X Z

batchBBoxes = zeros(numBatch, 6);
batchBBoxes(:, 1 : 3) = (bSubs - 1) .* BatchSize + 1; 
batchBBoxes(:, 4 : 6) = min(batchBBoxes(:, 1 : 3) + BatchSize - 1, imSize);

MIPZarrTmppaths = cellfun(@(x) sprintf('%s/%s_%s_%s.zarr', MIPPath, fsname, x, uuid), axis_strs, 'unif', 0);

zarr_done_flag = false(3, 1);
for i = 1 : 3
    zarr_done_flag(i) = (axis(i) == 0) | exist(MIPZarrpaths{i}, 'dir');
end

% if all mip zarr files exist, directly generate final MIPs
if all(zarr_done_flag) && ~mipSlab
    % collect results and generate MIPs 
    for i = 1 : 3
        if axis(i) == 0
            continue;
        end
        saveMIP_zarr(MIPZarrpaths{i}, MIPFullpaths{i}, dtype, (1 : 3) == i);
    end
    return;
end

if isempty(poolSize)
    poolSize = BatchSize;
end

% initialize zarr files
for i = 1 : 3
    if exist(MIPZarrTmppaths{i}, 'dir')
        continue;
    end
    axis_flag = false(3, 1);
    axis_flag(i) = true;
    
    outSize = imSize;
    outSize(axis_flag) = ceil(outSize(axis_flag) / poolSize(axis_flag));
    outSize(~axis_flag) = ceil(outSize(~axis_flag) ./ poolSize_1(~axis_flag));
    BlockSize_i = BatchSize;
    BlockSize_i(axis_flag) = ceil(BatchSize(axis_flag) / poolSize(axis_flag));
    BlockSize_i(~axis_flag) = ceil(BatchSize(~axis_flag) ./ poolSize_1(~axis_flag));
    
    createzarr(MIPZarrTmppaths{i}, dataSize=outSize, blockSize=BlockSize_i, dtype=dtype, zarrSubSize=zarrSubSize);
end

% set up parallel computing 
numBatch = size(batchBBoxes, 1);
taskSize = max(5, min(20, round(numBatch / 5000))); % the number of batches a job should process
numTasks = ceil(numBatch / taskSize);

maxJobNum = inf;
taskBatchNum = 1;

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    batchBBoxes_i = batchBBoxes(batchInds, :);
    
    zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;
    MIPZarrTmppaths_str = sprintf('{''%s''}', strjoin(MIPZarrTmppaths, ''','''));
    
    funcStrs{i} = sprintf(['MIP_block([%s],''%s'',%s,''%s'',%s,%s,''uuid'',''%s'',', ...
        '''debug'',%s)'], strrep(num2str(batchInds, '%d,'), ' ', ''), zarrFullpath, ...
        MIPZarrTmppaths_str, zarrFlagFullpath, strrep(mat2str(batchBBoxes_i), ' ', ','), ...
        strrep(mat2str([poolSize, poolSize_1]), ' ', ','), uuid, string(debug));
end

% submit jobs 
inputFullpaths = repmat({zarrFullpath}, numTasks, 1);
if parseCluster || ~parseParfor
    memAllocate = prod(BatchSize) * byteNum / 2^30 * (2 + (2 * (~mccMode)));
    minTaskJobNum = max(min(numTasks, 10), round(numTasks / 5));
    is_done_flag = false;
    for i = 1 : 3
        if all(is_done_flag)
            break;
        end
        is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'memAllocate', memAllocate * 2^(i-1), 'maxJobNum', maxJobNum, ...
            'taskBatchNum', taskBatchNum, 'minTaskJobNum', minTaskJobNum, 'masterCompute', masterCompute, ...
            'parseCluster', parseCluster, 'jobLogDir', jobLogDir, 'mccMode', mccMode, 'ConfigFile', ConfigFile);
    end
elseif parseParfor
    GPUJob = false;
    nworker = 12;
    is_done_flag= matlab_parfor_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'taskBatchNum', taskBatchNum, 'GPUJob', GPUJob, 'nworker', nworker, 'uuid', uuid);
end

if ~all(is_done_flag)
    error('Some blocks are not finished!')
end

for i = 1 : 3
    if exist(MIPZarrpaths{i}, 'dir') && exist(MIPZarrTmppaths{i}, 'dir')
        rmdir(MIPZarrpaths{i}, 's');
    end
    if exist(MIPZarrTmppaths{i}, 'dir')
        movefile(MIPZarrTmppaths{i}, MIPZarrpaths{i});
    end
end

if mipSlab
    return;
end

% collect results and generate MIPs 
for i = 1 : 3
    if axis(i) == 0
        continue;
    end
    saveMIP_zarr(MIPZarrpaths{i}, MIPFullpaths{i}, dtype, (1 : 3) == i);
end

end

