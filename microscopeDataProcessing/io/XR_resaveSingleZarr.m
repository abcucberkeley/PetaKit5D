function [] = XR_resaveSingleZarr(zarrFullpath, resultFullpath, varargin)
% 
% Author: Xiongtao Ruan (07/31/2025)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath', @ischar); 
ip.addRequired('resultFullpath', @ischar); 
ip.addParameter('inputBbox', [], @isnumeric); % inputBbox for input
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @isnumeric); % size to process in one batch 
ip.addParameter('shardSize', [512, 512, 512], @isnumeric); % size in one shard 
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('cpusPerTask', 1, @isscalar);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(zarrFullpath, resultFullpath, varargin{:});

pr = ip.Results;
inputBbox = pr.inputBbox;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if isempty(uuid)
    uuid = get_uuid();
end

% check if the group folder exist
if exist(resultFullpath, 'dir')
    disp('The output result exists, skip it!');
    return;
end

[resultPath, resultfsname] = fileparts(resultFullpath);
if strcmp(resultFullpath(end), '/')
    resultTmpPath = [resultFullpath(1 : end - 1), '_', uuid];
else
    resultTmpPath = [resultFullpath, '_', uuid];
end    

if ~exist(zarrFullpath, 'dir')
    error('The input file %s does not exist, pleast double check!', zarrFullpath);
end

zInfo = getZarrInfo(zarrFullpath);
dtype = zInfo.dtype;
sz = zInfo.size;

outSize = sz;
if ~isempty(inputBbox)
    outSize = inputBbox(4 : 6) - inputBbox(1 : 3) + 1;
end

% framework
blockSize = min(outSize, blockSize);
batchSize = min(outSize, batchSize);
batchSize = ceil(batchSize ./ blockSize) .* blockSize;
% without a cluster, neighboring batches are read and written together, so
% that their chunks are compressed in parallel (a batch of one chunk kept a
% few cores busy). Their flags name their size, so that the flags of a run
% with other batches and the same uuid are not taken for them.
flagName = 'batches';
if ~parseCluster
    localBatchSize = local_batch_size(outSize, batchSize, blockSize, dataTypeToByteNumber(dtype));
    if any(localBatchSize ~= batchSize)
        flagName = sprintf('batches_%d_%d_%d', localBatchSize);
    end
    batchSize = localBatchSize;
end
bSubSz = ceil(outSize ./ batchSize);
numBatch = prod(bSubSz);

sameBatchSize = false;
BorderSize = 0;

[batchBBoxes, regionBBoxes, localBBoxes] = XR_zarrChunkCoordinatesExtraction(outSize, ...
    'batchSize', batchSize, 'blockSize', blockSize, 'sameBatchSize',sameBatchSize, ...
    'BorderSize', BorderSize);

% bbox for the actual input data
if ~isempty(inputBbox)
    batchBBoxes = batchBBoxes + [inputBbox(1 : 3), inputBbox(1 : 3)] - 1;
end

% process for each file using distributed computing
fprintf('Process resaving zarr data...\n')

zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', resultPath, resultfsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end

% initialize zarr file
if exist(resultTmpPath, 'dir')
    tmp_zInfo = getZarrInfo(resultTmpPath);
    if any(tmp_zInfo.blockSize ~= blockSize) || any(tmp_zInfo.size ~= outSize)
        rmdir(resultTmpPath, 's');
        rmdir(zarrFlagPath, 's');
        mkdir(zarrFlagPath);
    end
end
if ~exist(resultTmpPath, 'dir')
    dimSeparator = '.';
    if prod(ceil(outSize ./ blockSize)) > 10000
        dimSeparator = '/';
    end
    createzarr(resultTmpPath, dataSize=outSize, blockSize=blockSize, dtype=dtype, ...
        dimSeparator=dimSeparator);
end

taskSize = max(5, ceil(numBatch / 5000)); % the number of batches a job should process
numTasks = ceil(numBatch / taskSize);

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    batchBBoxes_i = batchBBoxes(batchInds, :);
    regionBBoxes_i = regionBBoxes(batchInds, :);

    zarrFlagFullpath = sprintf('%s/%s_%d_%d.mat', zarrFlagPath, flagName, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;

    funcStrs{i} = sprintf(['resaveZarrBlock([%s],''%s'',''%s'',''%s'',%s,%s)'], ...
        mat2str_comma(batchInds), zarrFullpath, resultTmpPath, zarrFlagFullpath, ...
        mat2str_comma(batchBBoxes_i), mat2str_comma(regionBBoxes_i));
end

inputFullpaths = repmat({zarrFullpath}, numTasks, 1);

dtype = getImageDataType(zarrFullpath);
byteNum = dataTypeToByteNumber(dtype);
memAllocate = prod(batchSize) * byteNum / 2^30 * 1.5;
is_done_flag= generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
    funcStrs, finalOutFullpath=resultFullpath, cpusPerTask=cpusPerTask, memAllocate=memAllocate, ...
    parseCluster=parseCluster, mccMode=mccMode, configFile=configFile);

if ~all(is_done_flag)
    is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, finalOutFullpath=resultFullpath, cpusPerTask=cpusPerTask, memAllocate=memAllocate, ...
        parseCluster=parseCluster, masterCompute=masterCompute, mccMode=mccMode, configFile=configFile);
end

if ~all(is_done_flag)
    error('%d / %d tasks failed!', sum(~is_done_flag), numel(is_done_flag));
end

% another copy of this task may have published the result: keep it. The
% job framework then reports all blocks as done, so this copy's temporary
% result can be incomplete (a copy from another run has its own uuid)
if exist(resultFullpath, 'dir')
    if exist(resultTmpPath, 'dir')
        rmdir(resultTmpPath, 's');
    end
elseif exist(resultTmpPath, 'dir')
    movefile(resultTmpPath, resultFullpath);
end
if exist(zarrFlagPath, 'dir')
    rmdir(zarrFlagPath, 's');
end

end


function [batchSize] = local_batch_size(outSize, batchSize, blockSize, byteNum)
% batch size for processing in this MATLAB: neighboring batches together
% (along y, then x, then z), up to one chunk per computational thread, with
% their data within 1/16 of this process's share of the memory: the memory
% (the job's in a Slurm job, at most what is available now) times its share
% of the CPUs (a local process pool divides the threads among its workers),
% so that processes that resave at once stay within it together. With the
% group, cpp-zarr's chunk buffers take about twice its data.

mem = getSystemMemory();
if isunix && ~ismac
    % free and reclaimable memory now
    tok = regexp(fileread('/proc/meminfo'), 'MemAvailable:\s+(\d+) kB', 'tokens', 'once');
    if ~isempty(tok)
        mem = min(mem, str2double(tok{1}) / 2^20, 'includenan');
    end
end
nCores = str2double(getenv('SLURM_CPUS_ON_NODE'));
if isnan(nCores)
    nCores = feature('numcores');
end
share = min(1, maxNumCompThreads / nCores);
maxBatchNum = min(floor(maxNumCompThreads / prod(batchSize ./ blockSize)), ...
    floor(mem * share * 2^30 / 16 / (prod(batchSize) * byteNum)), 'includenan');
if ~(maxBatchNum > 1)
    return;
end

bSubSz = ceil(outSize ./ batchSize);
n = ones(1, 3);
for d = 1 : 3
    n(d) = min(bSubSz(d), floor(maxBatchNum / prod(n)));
end
batchSize = batchSize .* n;

end
