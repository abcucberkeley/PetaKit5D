function [is_done_flag] = XR_crop_zarr(zarrFullpath, cropFullpath, bbox, varargin)
% inplace cropping of zarr array (large)
%
% Author: Xiongtao Ruan (06/03/2022)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath', @(x) ischar(x) || iscell(x));
ip.addRequired('cropFullpath', @(x) ischar(x) || iscell(x));
ip.addRequired('bbox', @isnumeric);
ip.addParameter('pad', false, @islogical); % pad region that is outside the bbox
ip.addParameter('batchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256] , @isnumeric); % save as zarr
ip.addParameter('saveMIP', true, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(zarrFullpath, cropFullpath, bbox, varargin{:});

pr = ip.Results;
pad = pr.pad;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
saveMIP = pr.saveMIP;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
configFile = pr.configFile;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;

cropFullpath = strip(cropFullpath, 'right', filesep);
[cropPath, fsname, ext] = fileparts(cropFullpath);

if exist(cropFullpath, 'dir')
    disp('The output result exists, skip it!');
    return;
end

fprintf('Start Large-file Cropping for %s...\n', fsname);

if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(cropPath, jobLogDir);
end

tic
zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', cropPath, fsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end 

bim = blockedImage(zarrFullpath, 'Adapter', CZarrAdapter);
imSize = bim.Size;
dtype = bim.ClassUnderlying;

% crop for each block
if ~pad && any(bbox(4 : 6) > imSize)
    error('The bounding box is out of bound! Please check the bbox setting or turn on pad option.')
end
outSize = bbox(4 : 6) - bbox(1 : 3) + 1;
batchSize = min(outSize, batchSize);
blockSize = min(outSize, blockSize);
sameBatchSize = false;
BorderSize = 0;

[batchBBoxes, regionBBoxes, localBBoxes] = XR_zarrChunkCoordinatesExtraction(outSize, ...
    'batchSize', batchSize, 'blockSize', blockSize, 'sameBatchSize',sameBatchSize, ...
    'BorderSize', BorderSize);

% bbox for the actual input data
batchBBoxes = batchBBoxes + [bbox(1 : 3), bbox(1 : 3)] - 1;

% initialize zarr file
cropTempPath = sprintf('%s/%s_%s.zarr', cropPath, fsname, uuid);
if ~exist(cropTempPath, 'dir')
    dimSeparator = '.';
    if prod(ceil(outSize ./ blockSize)) > 10000
        dimSeparator = '/';
    end
    createzarr(cropTempPath, dataSize=outSize, blockSize=blockSize, dtype=dtype, dimSeparator=dimSeparator);
end

% set up parallel computing 
numBatch = size(batchBBoxes, 1);
taskSize = max(5, round(numBatch / 5000)); % the number of batches a job should process
numTasks = ceil(numBatch / taskSize);

maxJobNum = inf;
taskBatchNum = 1;

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    batchBBoxes_i = batchBBoxes(batchInds, :);
    regionBBoxes_i = regionBBoxes(batchInds, :);
    
    zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;
        
    funcStrs{i} = sprintf(['XR_crop_block([%s],''%s'',''%s'',''%s'',%s,%s,''uuid'',''%s'',', ...
        '''debug'',%s)'], strrep(num2str(batchInds, '%d,'), ' ', ''), zarrFullpath, ...
        cropTempPath, zarrFlagFullpath, strrep(mat2str(batchBBoxes_i), ' ', ','), ...
        strrep(mat2str(regionBBoxes_i), ' ', ','), uuid, string(debug));
end

% submit jobs 
inputFullpaths = repmat({zarrFullpath}, numTasks, 1);
if parseCluster || ~parseParfor
    memAllocate = prod(batchSize) * 4 * 2 / 1024^3;
    is_done_flag= generic_computing_frameworks_wrapper(inputFullpaths, ...
        outputFullpaths, funcStrs, 'cpusPerTask', cpusPerTask, 'maxJobNum', maxJobNum, ...
        'memAllocate', memAllocate, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, ...
        'parseCluster', parseCluster, mccMode=mccMode, configFile=configFile);

    if ~all(is_done_flag)
        is_done_flag= generic_computing_frameworks_wrapper(inputFullpaths, ...
            outputFullpaths, funcStrs, 'cpusPerTask', cpusPerTask, 'maxJobNum', maxJobNum, ...
            'memAllocate', memAllocate * 2, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, ...
            'parseCluster', parseCluster, mccMode=mccMode, configFile=configFile);
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

if exist(cropFullpath, 'dir') && exist(cropTempPath, 'dir')
    rmdir(cropFullpath, 's');
end
if exist(cropTempPath, 'dir')
    movefile(cropTempPath, cropFullpath);
end
if exist(zarrFlagPath, 'dir')
    rmdir(zarrFlagPath, 's');
end

% generate MIPs
if saveMIP
    XR_MIP_zarr(cropFullpath, 'axis', [1, 1, 1], parseCluster=parseCluster, mccMode=mccMode, configFile=configFile);
end

end

