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
% xruan (06/03/2023): add support for resampling
% xruan (08/02/2023): add support for MIP mask based cropping of input and output result


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpath', @(x) ischar(x) || iscell(x));
ip.addRequired('xyPixelSize'); 
ip.addRequired('dz'); 
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('Crop', false, @islogical);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('DSRCombined', true, @(x) islogical(x)); % combined processing 
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('Save16bit', false , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('SaveMIP', true , @islogical); % save MIP-z for ds and dsr. 
ip.addParameter('saveZarr', false , @islogical); % save as zarr
ip.addParameter('BatchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @isvector); % in y, x, z
ip.addParameter('zarrSubSize', [20, 20, 20], @isnumeric); % zarr subfolder size
ip.addParameter('inputBbox', [], @(x) isempty(x) || isvector(x));
ip.addParameter('taskSize', [], @isnumeric);
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x)); % resampling after rotation 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.addParameter('maskFns', {}, @iscell); % 2d masks to filter regions to deskew and rotate, in xy, xz, yz order
ip.addParameter('surffix', '', @ischar); % suffix for the folder
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 8, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(frameFullpath, xyPixelSize, dz, varargin{:});

pr = ip.Results;
Crop = pr.Crop;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
ObjectiveScan = pr.ObjectiveScan;
flipZstack = pr.flipZstack;
Save16bit = pr.Save16bit;
SaveMIP = pr.SaveMIP;
resample = pr.resample;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
zarrSubSize = pr.zarrSubSize;
inputBbox = pr.inputBbox;
taskSize = pr.taskSize;
Interp = pr.Interp;
maskFns = pr.maskFns;
surffix = pr.surffix;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

% decide zAniso
if ObjectiveScan
    zAniso = dz / xyPixelSize;
else
    theta = SkewAngle * pi / 180;
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
end
if iscell(frameFullpath)
    frameFullpath = frameFullpath{1};
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

if Save16bit
    dtype = 'uint16';
else
    dtype = 'single';
end

if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath, jobLogDir);
end

tic
zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', dsrPath, fsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end

rs = [1, 1, 1];
if ~isempty(resample)
    rs = resample;
end

% map input and output for xz
bimSize = getImageSize(frameFullpath);

% use MIP masks to decide the input and output boudning box
outputBbox = [];
if ~isempty(maskFns)
    fprintf('Compute input and out bounding boxes with MIP masks...\n')
    disp(maskFns(:));

    BorderSize = [100, 100, 100];
    BorderSize = max(10, round(BorderSize ./ rs));
    [maskInBbox, maskOutBbox] = XR_getDeskeRotateBoxesFromMasks(maskFns, BorderSize, ...
        xyPixelSize, dz, SkewAngle, Reverse, rs, inputBbox);

    inputBbox = maskInBbox;
    outputBbox = maskOutBbox;

    fprintf('Input bounding box: %s\n', mat2str(inputBbox));
    fprintf('Output bounding box: %s\n', mat2str(outputBbox));    
end

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
    dsrOutSize = round([ny, (nx-1)*cos(theta)+(nz-1)*zAniso/sin(abs(theta)), (nx-1)*sin(abs(theta))-4]);
else
    % exact proportions of rotated box
    dsrOutSize = round([ny, nx*cos(theta)+nz*zAniso*sin(abs(theta)), nz*zAniso*cos(theta)+nx*sin(abs(theta))]);
end
dsrOutSize = round(dsrOutSize ./ rs([1,2,3]));

outSize = dsrOutSize;
if isempty(outputBbox)
    outputBbox = [1, 1, 1, dsrOutSize];
end
outSize(2 : 3) = outputBbox(5 : 6) - outputBbox(2 : 3) + 1;

% change border size to +/-1 in y. 
% BorderSize = [2, 0, 0, 2, 0, 0];
BorderSize = [1, 0, 0, 1, 0, 0];

% set batches along y axis
outBatchSize = min(outSize, BatchSize);
outBatchSize(2 : 3) = outSize(2 : 3);
outBlockSize = min(imSize, BlockSize);
outBlockSize = min(outBatchSize, outBlockSize);

[outBatchBBoxes, outRegionBBoxes, outLocalBboxes] = XR_zarrChunkCoordinatesExtraction(dsrOutSize, ...
    BatchSize=outBatchSize, BlockSize=outBlockSize, bbox=outputBbox, SameBatchSize=false, BorderSize=BorderSize(1 : 3));

inBatchBboxes = outBatchBBoxes;
inBatchBboxes(:, 1) = round((inBatchBboxes(:, 1) - 1) .* rs(1)) + wdStart(1);
inBatchBboxes(:, 2) = wdStart(2);
inBatchBboxes(:, 3) = wdStart(3);
inBatchBboxes(:, 4) = min(inBatchBboxes(:, 1) + round((outBatchBBoxes(:, 4) - outBatchBBoxes(:, 1) + 1) .* rs(1) - 1), bimSize(1));
inBatchBboxes(:, 5) = inBatchBboxes(:, 2) + imSize(2) - 1;
inBatchBboxes(:, 6) = inBatchBboxes(:, 3) + imSize(3) - 1;

% initialize zarr file
if exist(dsrTmppath, 'dir')
    bim = blockedImage(dsrTmppath, 'Adapter', CZarrAdapter);
    if any(bim.BlockSize ~= outBlockSize) || any(bim.Size ~= outSize)
        rmdir(dsrTmppath, 's');
        rmdir(zarrFlagPath, 's');
        mkdir(zarrFlagPath);
    end
end
if ~exist(dsrTmppath, 'dir')
    createzarr(dsrTmppath, dataSize=outSize, blockSize=outBlockSize, dtype=dtype, zarrSubSize=zarrSubSize);                
end

% set up parallel computing 
numBatch = size(inBatchBboxes, 1);
if isempty(taskSize)
    taskSize = max(2, min(10, round(numBatch / 5000))); % the number of batches a job should process
end
numTasks = ceil(numBatch / taskSize);

maxJobNum = inf;
taskBatchNum = 1;

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    inBatchBBoxes_i = inBatchBboxes(batchInds, :);
    outBatchBBoxes_i = outBatchBBoxes(batchInds, :);
    outRegionBBoxes_i = outRegionBBoxes(batchInds, :);
    outLocalBBoxes_i = outLocalBboxes(batchInds, :);
    
    zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;
    Overwrite = false;
    
    funcStrs{i} = sprintf(['XR_deskewRotateBlock([%s],''%s'',''%s'',''%s'',%s,%s,%s,%s,', ...
        '%0.20d,%0.20d,''Overwrite'',%s,''SkewAngle'',%0.20d,''Reverse'',%s,', ...
        '''flipZstack'',%s,''Interp'',''%s'',''resample'',%s,''uuid'',''%s'',''debug'',%s)'], ...
        strrep(num2str(batchInds, '%d,'), ' ', ''), frameFullpath, dsrTmppath, ...
        zarrFlagFullpath, strrep(mat2str(inBatchBBoxes_i), ' ', ','), strrep(mat2str(outBatchBBoxes_i), ' ', ','), ...
        strrep(mat2str(outRegionBBoxes_i), ' ', ','), strrep(mat2str(outLocalBBoxes_i), ' ', ','), ...
        xyPixelSize, dz, string(Overwrite), SkewAngle, string(Reverse), string(flipZstack), ...
        Interp, strrep(mat2str(resample), ' ', ','), uuid, string(debug));
end

inputFullpaths = repmat({frameFullpath}, numTasks, 1);

% submit jobs
if parseCluster || ~parseParfor
    inSize = inBatchBboxes(1, 4 : 6) - inBatchBboxes(1, 1 : 3) + 1;
    outSize = outBatchBBoxes(1, 4 : 6) - outBatchBBoxes(1, 1 : 3) + 1;
    outSize(2 : 3) = max(outSize(2 : 3), dsrOutSize(2 : 3));
    memAllocate = prod(inSize) * 4 / 1024^3 * 2 + prod(outSize) * 4 / 1024^3 .* (prod(resample) * 0.6 + 1.5);

    is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', cpusPerTask, 'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, ...
        'masterCompute', masterCompute, 'parseCluster', parseCluster, 'memAllocate', memAllocate, ...
        'mccMode', mccMode, 'ConfigFile', ConfigFile);

    if ~all(is_done_flag)
        is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'cpusPerTask', cpusPerTask, 'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, ...
            'masterCompute', masterCompute, 'parseCluster', parseCluster, 'memAllocate', memAllocate, ...
            'mccMode', mccMode, 'ConfigFile', ConfigFile);
    end
elseif parseParfor
    is_done_flag= generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'clusterType', 'parfor', 'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, ...
        'GPUJob', GPUJob, 'uuid', uuid);
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
    XR_MIP_zarr(dsrFullpath, mccMode=mccMode, ConfigFile=ConfigFile);
end
toc

end

