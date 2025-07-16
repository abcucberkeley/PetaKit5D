function XR_unmix_channels_zarr(zarrFullpaths, unmixFactors, varargin)
% Unmixing a channel by subtracting another channel with a factor
%
% Required inputs:
%           zarrFullpaths : char or cell array. Directory paths for the zarr files for different channels.
%            unmixFactors : 1x#Channel numerial array. Unmixing factors for different channels.
%
% Parameters (as 'specifier'-value pairs):
%                mode : char (default: 'linear'). Mode of unmixing. linear: linear unmixing; gaussian: gaussian unmixing
%         unmixSigmas : empty or 1x#Channel (default: []). Sigmas for gaussian unmixing.
%       resultDirName : char (default: 'Unmixed'). Result directory under data paths.
%          channelInd : a number from 1 - #Channel (default: 1). The channel to unmix.
%           batchSize : 1x3 vector (default: [1024, 1024, 1024]). Batch size per unmixing task.
%           blockSize : 1x3 vector (default: [256, 256, 256]). Block/chunk size for zarr output.
%          borderSize : 1x3 vector (default: [0, 0, 0]). Padded border for each batch.
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%         parseParfor : true|false (default: false). Use matlab parfor for paralle processing.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%           jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%         cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%             mccMode : true|false (default: false). Use mcc mode.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%               debug : true|false (default: false). Debug mode reserved for future use.
%
% Author: Xiongtao Ruan (07/24/2022)
%
% xruan (02/07/2025): update to current way for job submission
% xruan (02/21/2025): add Gaussian smoothing method for unmixing

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('unmixFactors', @(x) isnumeric(x));
ip.addParameter('mode', 'linear', @ischar); % linear vs gaussian
ip.addParameter('unmixSigmas', [], @isnumeric); 
ip.addParameter('resultDirName', 'Unmixed', @ischar); 
ip.addParameter('channelInd', 1, @isnumeric); % unmix for which channel
ip.addParameter('batchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256] , @isvector); % in y, x, z
ip.addParameter('borderSize', [0, 0, 0] , @isvector); % in y, x, z
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 3, @isnumeric);
ip.addParameter('configFile', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(zarrFullpaths, unmixFactors, varargin{:});

pr = ip.Results;
mode = pr.mode;
unmixSigmas = pr.unmixSigmas;
pthstr = pr.pthstr;
channelInd = pr.channelInd;
BatchSize = pr.batchSize;
BlockSize = pr.blockSize;
borderSize = pr.borderSize;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
configFile = pr.configFile;
mccMode = pr.mccMode;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;

zarrFullpaths = strip(zarrFullpaths, 'right', filesep);
[dataPath, fsname, ext] = fileparts(zarrFullpaths{channelInd});

unmixPath = [dataPath, '/', pthstr, '/'];
if ~exist(unmixPath, 'dir')
    mkdir(unmixPath);
end

unmixFullpath = sprintf('%s/%s.zarr', unmixPath, fsname);
unmixTmppath = sprintf('%s/%s_%s.zarr', unmixPath, fsname, uuid);

done_flag = exist(unmixFullpath, 'file');
if all(done_flag)
    disp('The output results exist, skip it!');
    return;
end

fprintf('Start Large-file channel unmixing for %s...\n', fsname);

if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str] = checkSlurmCluster(dataPath, jobLogDir);
end

tic
zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', unmixPath, fsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end 

imSize = getImageSize(zarrFullpaths{channelInd});
dtype = getImageDataType(zarrFullpaths{channelInd});

% MIPs for each block
BatchSize = min(imSize, BatchSize);
BlockSize = min(BatchSize, BlockSize);
SameBatchSize = ~true;
[batchBBoxes, regionBBoxes, localBBoxes] = XR_zarrChunkCoordinatesExtraction(imSize, 'BatchSize', BatchSize, ...
    'BlockSize', BlockSize, 'SameBatchSize', SameBatchSize, 'BorderSize', borderSize);

% initialize zarr files
if ~exist(unmixTmppath, 'dir')
    outSize = imSize;
    dimSeparator = '.';
    if prod(ceil(outSize ./ BlockSize)) > 10000
        dimSeparator = '/';
    end    
    createzarr(unmixTmppath, dataSize=outSize, BlockSize=BlockSize, dtype=dtype, dimSeparator=dimSeparator);
end

% set up parallel computing 
numBatch = size(batchBBoxes, 1);
taskSize = max(1, min(ceil(numBatch / 10000), 20)); % the number of batches a job should process
numTasks = ceil(numBatch / taskSize);

taskBatchNum = 1;

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    batchBBoxes_i = batchBBoxes(batchInds, :);
    regionBBoxes_i = regionBBoxes(batchInds, :);
    localBBoxes_i = localBBoxes(batchInds, :);
    
    zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;
    zarrFullpaths_str = sprintf('{''%s''}', strjoin(zarrFullpaths, ''','''));
    
    switch mode
        case 'linear'
            funcStrs{i} = sprintf(['unmix_channels_block(%s,%s,''%s'',%s,''%s'',%s,''uuid'',''%s'',', ...
                '''debug'',%s)'], strrep(mat2str(batchInds), ' ', ','), zarrFullpaths_str, ...
                unmixTmppath, strrep(mat2str(unmixFactors), ' ', ','), zarrFlagFullpath, ...
                strrep(mat2str(batchBBoxes_i), ' ', ','), uuid, string(debug));
            memFactor = numel(zarrFullpaths) + 1;
        case 'gaussian'
            funcStrs{i} = sprintf(['unmix_channels_gaussian_block(%s,%s,''%s'',', ...
                '%s,%s,''%s'',%s,%s,%s,''uuid'',''%s'',''debug'',%s)'], strrep(mat2str(batchInds), ' ', ','), ...
                zarrFullpaths_str, unmixTmppath, strrep(mat2str(unmixFactors), ' ', ','), ...
                strrep(mat2str(unmixSigmas), ' ', ','), zarrFlagFullpath, strrep(mat2str(batchBBoxes_i), ' ', ','), ...
                strrep(mat2str(regionBBoxes_i), ' ', ','), strrep(mat2str(localBBoxes_i), ' ', ','), ...
                uuid, string(debug));
            memFactor = (numel(zarrFullpaths) + 1) * 1.5;
    end
end

% submit jobs
memAllocate = prod(BatchSize) * 4 / 1024^3 * memFactor;
inputFullpaths = repmat({zarrFullpaths{channelInd}}, numTasks, 1);
if parseCluster || ~parseParfor 
    is_done_flag= generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, finalOutFullpath=unmixFullpath, parseCluster=parseCluster, ...
        masterCompute=masterCompute, taskBatchNum=taskBatchNum, cpusPerTask=cpusPerTask, ...
        memAllocate=memAllocate, mccMode=mccMode, ConfigFile=configFile);

    if ~all(is_done_flag)
        is_done_flag= generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, finalOutFullpath=unmixFullpath, parseCluster=parseCluster, ...
            masterCompute=masterCompute, taskBatchNum=taskBatchNum, cpusPerTask=cpusPerTask, ...
            memAllocate=memAllocate, mccMode=mccMode, ConfigFile=configFile);
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

if exist(unmixFullpath, 'dir') && exist(unmixTmppath, 'dir')
    rmdir(unmixFullpath, 's');
end
if exist(unmixTmppath, 'dir')
    movefile(unmixTmppath, unmixFullpath);
end


end


