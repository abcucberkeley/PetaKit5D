function [] = rescale_zarr_file(frameFullpath, rsPath, varargin)
% RL decon for big data that can not be loaded to memory, the idea is to let
% workers load certain chunks and decon the chunk and save the results to
% the right coordinates in the disk.
% 
% xruan (11/13/2021): add support for parfor based parallel computing


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('outFullpath', @(x) ischar(x) || iscell(x));
ip.addParameter('save16bit', true , @islogical);
ip.addParameter('rsFactor', 60000 , @isnumeric);
ip.addParameter('rsRange', [0, 65535] , @isnumeric);
ip.addParameter('batchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256] , @isvector); % in y, x, z
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('parseCluster', true, @islogical); %  
ip.addParameter('parseParfor', false, @islogical); %  
ip.addParameter('GPUJob', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('masterCompute', true, @islogical); % 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(frameFullpath, rsPath, varargin{:});

pr = ip.Results;    

% parameters
save16bit = pr.save16bit;
rsFactor = pr.rsFactor;
rsRange = pr.rsRange;
GPUJob = pr.GPUJob;
useGPU = true;

% check if background information available, if not, estimate background
% info. Currently use 99. 
% simplified version related options

batchSize = pr.batchSize;
BlockSize = pr.BlockSize;

tic
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
uuid = pr.uuid;
if isempty(uuid)
    uuid = get_uuid();
end

[dataPath, fsname, ext] = fileparts(frameFullpath);
rsFullpath = sprintf('%s/%s.zarr', rsPath, fsname);

if exist(rsFullpath, 'dir')
    disp('The output result exists, skip it!');
    return;
end

fprintf('Start Large-file RL Decon for %s...\n', fsname);

if strcmp(rsFullpath(end), '/')
    rsTmppath = [rsFullpath(1 : end-1), '_', uuid];
else
    rsTmppath = [rsFullpath, '_', uuid];
end


tic
fprintf(['reading ' fsname '...\n'])
switch ext
    case {'.tif', '.tiff'}
        bim = blockedImage(frameFullpath, 'Adapter', MPageTiffAdapter);
    case '.zarr'
        bim = blockedImage(frameFullpath, 'Adapter', ZarrAdapter);
end
toc

% not consider edge erosion for now

% dtype = class(im);
if save16bit
    dtype = 'uint16';
else
    dtype = 'single';
end

if parseCluster
    if  ~exist(jobLogDir, 'dir')
        warning('The job log directory does not exist, use %s/job_logs as job log directory.', deconPath)
        jobLogDir = sprintf('%s/job_logs', deconPath);
        mkdir(jobLogDir);
    end
end

tic

zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', rsPath, fsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end

imSize = bim.Size;

% BlockSize = nv_bim.BlockSize;
sameBatchSize = false;
BorderSize = 0;
[batchBBoxes, regionBBoxes] = XR_zarrChunkCoordinatesExtraction(imSize, 'batchSize', batchSize, ...
    'BlockSize', BlockSize, 'sameBatchSize', sameBatchSize, 'BorderSize', BorderSize);


% initialize zarr file
init_val = zeros(1, dtype);
if ~exist(rsTmppath, 'dir')
    rs_bim = blockedImage(rsTmppath, imSize, BlockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');
    rs_bim.Adapter.close();
end

taskSize = 5; % the number of batches a job should process
numBatch = size(batchBBoxes, 1);
numTasks = ceil(numBatch / taskSize);

if GPUJob
    maxJobNum = inf;
    cpusPerTask = 5;
    cpuOnlyNodes = false;
    taskBatchNum = 5;
    SlurmParam = '-p abc --qos abc_normal -n1 --mem=167G --gres=gpu:1';
else
    maxJobNum = inf;
    cpusPerTask = 1;
    cpuOnlyNodes = true;
    taskBatchNum = 1;
    SlurmParam = '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M';
end

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    batchBBoxes_i = batchBBoxes(batchInds, :);
    regionBBoxes_i = regionBBoxes(batchInds, :);
    
    zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;
    
    funcStrs{i} = sprintf(['rescale_zarr_block([%s],''%s'',''%s'',''%s'',%0.20f,%s,', ...
        '%s,%s)'], ...
        strrep(num2str(batchInds, '%d,'), ' ', ''), frameFullpath, rsTmppath, zarrFlagFullpath, ...
        rsFactor, strrep(mat2str(rsRange), ' ', ','), strrep(mat2str(batchBBoxes_i), ' ', ','), ...
        strrep(mat2str(regionBBoxes_i), ' ', ','));
end

inputFullpaths = repmat({frameFullpath}, numTasks, 1);
% submit jobs
if parseCluster
    is_done_flag= slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
        'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);

    if ~all(is_done_flag)
        slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
            'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);
    end
elseif parseParfor
    is_done_flag= matlab_parfor_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'GPUJob', GPUJob);
end

if exist(rsFullpath, 'dir')
    rmdir(rsFullpath, 's');
end
movefile(rsTmppath, rsFullpath);

% generate MIP z file
deconMIPPath = sprintf('%s/MIPs/', deconPath);
if ~exist(deconMIPPath, 'dir')
    mkdir(deconMIPPath);
    fileattrib(deconMIPPath, '+w', 'g');
end
deconMIPname = sprintf('%s%s_MIP_z.tif', deconMIPPath, fsname);
saveMIP_zarr(rsFullpath, deconMIPname);
toc

end

