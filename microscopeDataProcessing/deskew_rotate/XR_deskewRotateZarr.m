function [is_done_flag] = XR_deskewRotateZarr(frameFullpath, xyPixelSize, dz, varargin)
% Deskew and/or rotate data for a single zarr file that cannot be fitted to
% memory
% 
%
% Author: Xiongtao Ruan (02/16/2022)
%
% Based on XR_deskewRotateFrame.m
%
% xruan (12/14/2022): add support for input bbox, that is, crop the data before dsr


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('xyPixelSize'); 
ip.addRequired('dz'); 
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('Crop', false, @islogical);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('Save16bit', false , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('SaveMIP', true , @islogical); % save MIP-z for ds and dsr. 
ip.addParameter('saveZarr', false , @islogical); % save as zarr
ip.addParameter('BatchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @isvector); % in y, x, z
ip.addParameter('inputBbox', [], @(x) isempty(x) || isvector(x));
ip.addParameter('taskSize', [], @isnumeric);
ip.addParameter('DSRCombined', true, @islogical); % combined processing 
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x)); % resampling after rotation 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.addParameter('surffix', '', @ischar); % suffix for the folder
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpuOnlyNodes', false, @islogical);
ip.addParameter('cpusPerTask', 8, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(frameFullpath, xyPixelSize, dz, varargin{:});

pr = ip.Results;
Crop = pr.Crop;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
ObjectiveScan = pr.ObjectiveScan;
flipZstack = pr.flipZstack;
Save16bit = pr.Save16bit;
SaveMIP = pr.SaveMIP;
DSRCombined = pr.DSRCombined;
resample = pr.resample;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
inputBbox = pr.inputBbox;
taskSize = pr.taskSize;
Interp = pr.Interp;
surffix = pr.surffix;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpuOnlyNodes = pr.cpuOnlyNodes;
cpusPerTask = pr.cpusPerTask;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;

% decide zAniso
if ObjectiveScan
    zAniso = dz / xyPixelSize;
else
    theta = SkewAngle * pi / 180;
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
end

if ~exist(frameFullpath, 'dir')
    error('The input zarr file %s does not exist!', frameFullpath);
end

[dataPath, fsname, ext] = fileparts(frameFullpath);
dsrPath = [dataPath, '/DSR/'];
if ~exist(dsrPath, 'dir')
    mkdir(dsrPath);
end

dsrFullpath = sprintf('%s/%s.zarr', dsrPath, fsname);

if exist(dsrFullpath, 'dir')
    disp('The output result exists, skip it!');
    return;
end

fprintf('Start Large-file Deskew and Rotate for %s...\n', fsname);

if strcmp(dsrFullpath(end), '/')
    dsrTmppath = [dsrFullpath(1 : end-1), '_', uuid];
else
    dsrTmppath = [dsrFullpath, '_', uuid];
end

tic
fprintf(['reading ' fsname '...\n'])
switch ext
    case {'.tif', '.tiff'}
        bim = blockedImage(frameFullpath, 'Adapter', MPageTiffAdapter);
    case '.zarr'
        bim = blockedImage(frameFullpath, 'Adapter', CZarrAdapter);
end
toc

if Save16bit
    dtype = 'uint16';
else
    dtype = 'single';
end

if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
end

tic
zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', dsrPath, fsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end

% map input and output for xz
bimSize = bim.Size;
if ~isempty(inputBbox)
    wdStart = inputBbox(1 : 3);
    imSize = inputBbox(4 : 6) - wdStart + 1;
else
    wdStart = [1, 1, 1];    
    imSize = bimSize;
end

ny = imSize(1);
nx = imSize(2);
nz = imSize(3);

if ~ObjectiveScan
    % outSize = round([ny nxDs/cos(theta) h]);
    % calculate height; first & last 2 frames have interpolation artifacts
    outSize = round([ny, (nx-1)*cos(theta)+(nz-1)*zAniso/sin(abs(theta)), (nx-1)*sin(abs(theta))-4]);
else
    % exact proportions of rotated box
    outSize = round([ny, nx*cos(theta)+nz*zAniso*sin(abs(theta)), nz*zAniso*cos(theta)+nx*sin(abs(theta))]);
end
% change border size to +/-2 in y. 
BorderSize = [2, 0, 0, 2, 0, 0];

% set batches along y axis
BatchSize = min(imSize, BatchSize);
BatchSize(2 : 3) = imSize(2 : 3);
BlockSize = min(imSize, BlockSize);
BlockSize = min(BatchSize, BlockSize);

bSubSz = ceil(imSize ./ BatchSize);
numBatch = prod(bSubSz);

[Y, X, Z] = ndgrid(1 : bSubSz(1), 1 : bSubSz(2), 1 : bSubSz(3));
bSubs = [Y(:), X(:), Z(:)];
clear Y X Z

batchBBoxes = zeros(numBatch, 6);
regionBBoxes = zeros(numBatch, 6);

batchBBoxes(:, 1 : 3) = (bSubs - 1) .* BatchSize + wdStart; 
batchBBoxes(:, 4 : 6) = min(batchBBoxes(:, 1 : 3) + BatchSize - 1, imSize + wdStart - 1);

borderSizes(:, 1 : 3) = batchBBoxes(:, 1 : 3) - max(1, batchBBoxes(:, 1 : 3) - BorderSize(1 : 3));
borderSizes(:, 4 : 6) = min(bimSize, batchBBoxes(:, 4 : 6) + BorderSize(4 : 6)) - batchBBoxes(:, 4 : 6);

batchBBoxes(:, 1 : 3) = batchBBoxes(:, 1 : 3) - borderSizes(:, 1 : 3);
batchBBoxes(:, 4 : 6) = batchBBoxes(:, 4 : 6) + borderSizes(:, 4 : 6);

regionBBoxes(:, 1) = (bSubs(:, 1) - 1) .* BatchSize(1) + 1; 
regionBBoxes(:, 2 : 3) = 1;
regionBBoxes(:, 4 : 6) = min(regionBBoxes(:, 1 : 3) + [BatchSize(1), outSize(2 : 3)] - 1, outSize);

% initialize zarr file
init_val = zeros(1, dtype);
if ~exist(dsrTmppath, 'dir')
    try
        dsr_bim = blockedImage(dsrTmppath, outSize, BlockSize, init_val, "Adapter", CZarrAdapter, 'Mode', 'w');
    catch ME
        disp(ME)
        dsr_bim = blockedImage(dsrTmppath, outSize, BlockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');
    end
    dsr_bim.Adapter.close();
end

% set up parallel computing 
numBatch = size(batchBBoxes, 1);
if isempty(taskSize)
    taskSize = max(1, min(10, round(numBatch / 5000))); % the number of batches a job should process
end
numTasks = ceil(numBatch / taskSize);

maxJobNum = inf;
taskBatchNum = 1;
SlurmParam = '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M';

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    batchBBoxes_i = batchBBoxes(batchInds, :);
    borderSizes_i = borderSizes(batchInds, :);
    regionBBoxes_i = regionBBoxes(batchInds, :);
    
    zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;
    Overwrite = false;
    
    funcStrs{i} = sprintf(['XR_deskewRotateBlock([%s],''%s'',''%s'',''%s'',%s,%s,%s,', ...
        '%0.20d,%0.20d,''Overwrite'',%s,''SkewAngle'',%0.20d,''Reverse'',%s,', ...
        '''flipZstack'',%s,''Interp'',''%s'',''uuid'',''%s'',''debug'',%s)'], ...
        strrep(num2str(batchInds, '%d,'), ' ', ''), frameFullpath, ...
        dsrTmppath, zarrFlagFullpath, strrep(mat2str(batchBBoxes_i), ' ', ','), ...
        strrep(mat2str(regionBBoxes_i), ' ', ','), strrep(mat2str(borderSizes_i), ' ', ','), ...
        xyPixelSize, dz, string(Overwrite), SkewAngle, string(Reverse), string(flipZstack), ...
        Interp, uuid, string(debug));
end

inputFullpaths = repmat({frameFullpath}, numTasks, 1);

% submit jobs
if parseCluster || ~parseParfor
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
        funcStrs, 'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'GPUJob', GPUJob, 'uuid', uuid);
end

if exist(dsrFullpath, 'dir') && exist(dsrTmppath, 'dir')
    rmdir(dsrFullpath, 's');
end
if exist(dsrTmppath, 'dir')
    movefile(dsrTmppath, dsrFullpath);
end

% generate MIP z file
if SaveMIP
    dsrMIPPath = sprintf('%s/MIPs/', dsrPath);
    if ~exist(dsrMIPPath, 'dir')
        mkdir(dsrMIPPath);
        fileattrib(dsrMIPPath, '+w', 'g');
    end
    % dsrMIPname = sprintf('%s%s_MIP_z.tif', dsrMIPPath, fsname);
    % saveMIP_zarr(dsrFullpath, dsrMIPname);
    XR_MIP_zarr(dsrFullpath);
end
toc

end

