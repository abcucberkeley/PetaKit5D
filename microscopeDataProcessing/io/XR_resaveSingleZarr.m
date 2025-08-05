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

    zarrFlagFullpath = sprintf('%s/batches_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
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

if exist(resultFullpath, 'dir')
    rmdir(resultFullpath, 's');
end
movefile(resultTmpPath, resultFullpath);
rmdir(zarrFlagPath, 's');

end
