function [] = N5ToZarr(n5Fullpath, varargin)
% convert N5 to zarr
% 
% Author: Xiongtao Ruan (09/17/2023)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('n5Fullpath'); 
ip.addParameter('resultDirStr', 'zarr/', @ischar);
ip.addParameter('bbox', [], @isnumeric); % bbox for input
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @isnumeric); % size to process in one batch 
ip.addParameter('zarrSubSize', [20, 20, 20], @isnumeric);
ip.addParameter('save16bit', true , @islogical); % save result data as 16 bit
ip.addParameter('flipEmptyValue', false , @islogical); % 65535 to 0 for empty region
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(n5Fullpath, varargin{:});

pr = ip.Results;
resultDirStr = pr.resultDirStr;
bbox = pr.bbox;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
zarrSubSize = pr.zarrSubSize;
flipEmptyValue = pr.flipEmptyValue;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if isempty(uuid)
    uuid = get_uuid();
end

% check if the folder exist
[dataPath, fsname, ext] = fileparts(n5Fullpath);
outPath = [dataPath, '/', resultDirStr];
if ~exist(outPath, 'dir')
    mkdir(outPath);
    if ~ispc
        fileattrib(outPath, '+w', 'g');
    end
end
resultFullpath = sprintf('%s/%s.zarr', outPath, fsname);
if exist(resultFullpath, 'dir')
    disp('The output result exists, skip it!');
    return;
end

if strcmp(resultFullpath(end), '/')
    outTmpPath = [resultFullpath(1 : end - 1), '_', uuid];
else
    outTmpPath = [resultFullpath, '_', uuid];
end    

if ~exist(n5Fullpath, 'dir')
    error('The input file %s does not exist, pleast double check!', n5Fullpath);
end
bim = blockedImage(n5Fullpath, 'Adapter', N5Adapter);
dtype = bim.ClassUnderlying;
sz = bim.Size;

if ~isempty(bbox)
    wdStart = bbox(1 : 3);
    imSize = bbox(4 : 6) - wdStart + 1;
else
    wdStart = [1, 1, 1];    
    imSize = sz;
end
outSize = imSize;

% framework
blockSize = min(outSize, blockSize);
batchSize = min(outSize, batchSize);
batchSize = ceil(batchSize ./ blockSize) .* blockSize;
bSubSz = ceil(outSize ./ batchSize);
numBatch = prod(bSubSz);

BorderSize = [0, 0, 0];

[outBatchBBoxes, ~, ~] = XR_zarrChunkCoordinatesExtraction(outSize, batchSize=batchSize, ...
    BlockSize=blockSize, sameBatchSize=false, BorderSize=BorderSize);

inBatchBboxes = outBatchBBoxes + [wdStart, wdStart] - 1;

% process for each block based on all BlockInfo use distributed computing
fprintf('Convert n5 file %s to zarr ...\n', n5Fullpath);

zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', outPath, fsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end

% initialize zarr file
if exist(outTmpPath, 'dir')
    bim = blockedImage(outTmpPath, 'Adapter', CZarrAdapter);
    if any(bim.BlockSize ~= blockSize) || any(bim.Size ~= outSize)
        rmdir(outTmpPath, 's');
        rmdir(zarrFlagPath, 's');
        mkdir(zarrFlagPath);
    end
end
if ~exist(outTmpPath, 'dir')
    createzarr(outTmpPath, dataSize=outSize, blockSize=blockSize, dtype=dtype, zarrSubSize=zarrSubSize);
end

taskSize = max(10, ceil(numBatch / 5000)); % the number of batches a job should process
numTasks = ceil(numBatch / taskSize);

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    inBatchBBoxes_i = inBatchBboxes(batchInds, :);
    outBatchBBoxes_i = outBatchBBoxes(batchInds, :);
    
    zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;
    Overwrite = false;
    
    funcStrs{i} = sprintf(['N5ToZarrBlock([%s],''%s'',''%s'',''%s'',[%s],[%s],', ...
        '''Overwrite'',%s,''flipEmptyValue'',%s)'], strrep(num2str(batchInds, '%d,'), ' ', ''), ...
        n5Fullpath, outTmpPath, zarrFlagFullpath, strrep(mat2str(inBatchBBoxes_i), ' ', ','), ...
        strrep(mat2str(outBatchBBoxes_i), ' ', ','), string(Overwrite), string(flipEmptyValue));
end

inputFullpaths = repmat({n5Fullpath}, numTasks, 1);

memAllocate = prod(batchSize) * 4 / 2^30 * 2.5;
is_done_flag= generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
    funcStrs, parseCluster=parseCluster, masterCompute=masterCompute, cpusPerTask=cpusPerTask, ...
    memAllocate=memAllocate, mccMode=mccMode, configFile=configFile);

if ~all(is_done_flag)
    generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, funcStrs, ...
        parseCluster=parseCluster, masterCompute=masterCompute, cpusPerTask=cpusPerTask, ...
        memAllocate=memAllocate * 2, mccMode=mccMode, configFile=configFile);
end

if exist(resultFullpath, 'dir')
    rmdir(resultFullpath, 's');
end
movefile(outTmpPath, resultFullpath);
rmdir(zarrFlagPath, 's');

fprintf('Done!\n')
toc

end

