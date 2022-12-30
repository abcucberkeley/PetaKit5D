function [] = RLdecon_large_in_place(frameFullpath, xyPixelSize, dz, deconPath, PSF, varargin)
% RL decon for big data that can not be loaded to memory, the idea is to let
% workers load certain chunks and decon the chunk and save the results to
% the right coordinates in the disk.
% 
% xruan (11/13/2021): add support for parfor based parallel computing


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('pixelSize', @isnumeric);
ip.addRequired('dz', @isnumeric);
ip.addRequired('deconPath', @(x) ischar(x) || isempty(x));
ip.addRequired('PSF', @(x) ischar(x) || isempty(x));
ip.addParameter('Save16bit', false , @islogical);
ip.addParameter('Deskew', false , @islogical);
ip.addParameter('SkewAngle', -32.45 , @isnumeric);
ip.addParameter('flipZstack', false, @islogical); 
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('dzPSF', 0.1 , @isnumeric); %in um
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @isnumeric); 
ip.addParameter('OTFCumThresh', 0.9, @isnumeric); % OTF cumutative sum threshold
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x)); % decon in skewed space
ip.addParameter('fixIter', false, @islogical); % CPU Memory in Gb
ip.addParameter('BatchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256] , @isvector); % in y, x, z
ip.addParameter('deconMaskFns', {}, @iscell); % 2d masks to filter regions to decon, in xy, xz, yz orders
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('cpusPerTask', 4, @isnumeric);
ip.addParameter('GPUJob', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('psfGen', true, @islogical); % psf generation

ip.parse(frameFullpath, xyPixelSize, dz, deconPath, PSF, varargin{:});

pr = ip.Results;    

% parameters
dzPSF = pr.dzPSF;
DeconIter = pr.DeconIter;
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
skewed = pr.skewed;
Deskew = pr.Deskew;
SkewAngle = pr.SkewAngle;
flipZstack = pr.flipZstack;
Save16bit = pr.Save16bit;
GPUJob = pr.GPUJob;
useGPU = GPUJob;
psfGen = pr.psfGen;

% check if background information available, if not, estimate background
% info. Currently use 99. 
Background = pr.Background;
if isempty(Background)
    Background = 100;
end
% simplified version related options
fixIter = pr.fixIter;
debug = pr.debug;

BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
deconMaskFns = pr.deconMaskFns;

tic
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes =  pr.cpuOnlyNodes;
uuid = pr.uuid;
if isempty(uuid)
    uuid = get_uuid();
end

[dataPath, fsname, ext] = fileparts(frameFullpath);
deconFullpath = sprintf('%s/%s.zarr', deconPath, fsname);

if exist(deconFullpath, 'dir')
    disp('The output result exists, skip it!');
    return;
end

fprintf('Start Large-file RL Decon for %s...\n', fsname);

if strcmp(deconFullpath(end), '/')
    deconTmppath = [deconFullpath(1 : end-1), '_', uuid];
else
    deconTmppath = [deconFullpath, '_', uuid];
end

% generate psf and use the cropped psf to decide the border size
pp = readtiff(PSF);
if psfGen
    medFactor = 1.5;
    PSFGenMethod = 'masked';
    if Deskew
        dz_data_ratio = sind(SkewAngle);
    else
        dz_data_ratio = 1;
    end
    psf = psf_gen_new(pp, dzPSF, dz*dz_data_ratio, medFactor, PSFGenMethod);
    py = find(squeeze(sum(psf, [2, 3])));
    px = find(squeeze(sum(psf, [1, 3])));
    pz = find(squeeze(sum(psf, [1, 2])));
    cropSz = [min(py(1) - 1, size(psf, 1) - py(end)), min(px(1) - 1, size(psf, 2) - px(end)), min(pz(1) - 1, size(psf, 3) - pz(end))] - 1;
    cropSz = max(0, cropSz);
    bbox = [cropSz + 1, size(psf) - cropSz];
    psf = psf(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
    
    % write generated psf to disk in case it will be deleted in the c
    psfgen_folder = sprintf('%s/psfgen/', deconPath);
    mkdir(psfgen_folder);
    [~, psf_fsn] = fileparts(PSF);
    psfgen_filename = sprintf('%s/%s.tif', psfgen_folder, psf_fsn);
    if ~exist(psfgen_filename, 'file')
        writetiff(psf, psfgen_filename);
    end
else
    psf = pp;
end


tic
fprintf(['reading ' fsname '...\n'])

% not consider edge erosion for now

% dtype = class(im);
if Save16bit
    dtype = 'uint16';
else
    dtype = 'single';
end

if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
end

zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', deconPath, fsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end

imSize = getImageSize(frameFullpath);
toc

tic
% BlockSize = nv_bim.BlockSize;
SameBatchSize = true;
BorderSize = round((size(psf) + 10) / 2);
[batchBBoxes, regionBBoxes] = XR_zarrChunkCoordinatesExtraction(imSize, 'BatchSize', BatchSize, ...
    'BlockSize', BlockSize, 'SameBatchSize', SameBatchSize, 'BorderSize', BorderSize);
scaleFactor = 1.0;

% initialize zarr file
init_val = zeros(1, dtype);
if ~exist(deconTmppath, 'dir')
    try
        decon_bim = blockedImage(deconTmppath, imSize, BlockSize, init_val, "Adapter", CZarrAdapter, 'Mode', 'w');
    catch ME
        disp(ME)
        decon_bim = blockedImage(deconTmppath, imSize, BlockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');
    end        
    decon_bim.Adapter.close();
end

taskSize = 20; % the number of batches a job should process
numBatch = size(batchBBoxes, 1);

if GPUJob
    if parseCluster
        taskSize = max(100, round(numBatch / 5000));
    end

    maxJobNum = inf;
    cpusPerTask = 4;
    cpuOnlyNodes = false;
    taskBatchNum = 1;
    SlurmParam = '-p abc_a100 --qos abc_normal -n1 --mem=125G --gres=gpu:1';
    masterCompute = false;
else
    maxJobNum = inf;
    % cpusPerTask = 24;
    cpuOnlyNodes = false;
    taskBatchNum = 1;
    SlurmParam = '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M';
end

numTasks = ceil(numBatch / taskSize);

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    batchBBoxes_i = batchBBoxes(batchInds, :);
    regionBBoxes_i = regionBBoxes(batchInds, :);
    
    zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;
    Overwrite = false;
    deconMaskFns_str = sprintf('{''%s''}', strjoin(deconMaskFns, ''','''));

    funcStrs{i} = sprintf(['RLdecon_for_zarr_block([%s],''%s'',''%s'',''%s'',''%s'',%s,%s,', ...
        '%0.20d,%0.20d,''Save16bit'',%s,''Overwrite'',%s,''SkewAngle'',%0.20d,''flipZstack'',%s,', ...
        '''Background'',%0.20d,''dzPSF'',%0.20d,''DeconIter'',%d,''RLMethod'',''%s'',', ...
        '''wienerAlpha'',%0.20f,''OTFCumThresh'',%0.20f,''skewed'',[%s],''fixIter'',%s,', ...
        '''scaleFactor'',%d,''useGPU'',%s,''deconMaskFns'',%s,''uuid'',''%s'',', ...
        '''debug'',%s,''psfGen'',%s)'], ...
        strrep(num2str(batchInds, '%d,'), ' ', ''), frameFullpath, PSF, deconTmppath, zarrFlagFullpath, ...
        strrep(mat2str(batchBBoxes_i), ' ', ','), strrep(mat2str(regionBBoxes_i), ' ', ','), xyPixelSize, ...
        dz, string(Save16bit), string(Overwrite), SkewAngle, string(flipZstack), Background, dzPSF, ...
        DeconIter, RLMethod, wienerAlpha, OTFCumThresh, string(skewed), string(fixIter), scaleFactor, ...
        string(useGPU), deconMaskFns_str, uuid, string(debug), string(psfGen));
end

inputFullpaths = repmat({frameFullpath}, numTasks, 1);
% submit jobs
if ~parseParfor
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
        funcStrs, 'taskBatchNum', taskBatchNum, 'GPUJob', GPUJob, 'uuid', uuid);
end

if ~all(is_done_flag, 'all')
    error('The decon of some chunks are missing!');
end

if exist(deconFullpath, 'dir') && exist(deconTmppath, 'dir')
    rmdir(deconFullpath, 's');
end
if exist(deconTmppath, 'dir')
    movefile(deconTmppath, deconFullpath);
end

% generate MIP z file
deconMIPPath = sprintf('%s/MIPs/', deconPath);
if ~exist(deconMIPPath, 'dir')
    mkdir(deconMIPPath);
    fileattrib(deconMIPPath, '+w', 'g');
end
XR_MIP_zarr(deconFullpath, 'axis', [1, 1, 1]);
toc

end

