function XR_ExM_PunctaRemoval_zarr(zarrFullpath, varargin)
% puncta removal for a zarr file based on GU_ExM_calcPunctaDensityVol.m
%
%
% Required inputs:
%        zarrFullpath : char or cell array. Directory path(s) to the input Zarr file(s).
%
% Parameters (as 'specifier'-value pairs):
%                mode : char or cell (default: 'updated'). Processing mode. 
%                       Options include: 'updated': updated method, 'bg_interpolated': based on interpolated background, and 'bg_estimation': background estimation.
%          bgFullpath : char (default: ''). File path to user-provided background data.
%            bgFactor : numeric (default: 1.0). Scaling factor applied to user-provided background data.
%               Sigma : numeric (default: 2.5). Sigma value for 3D Gaussian filtering.
%          OTSUMaxPer : numeric (default: 99.9). Maximum percentile of data used for Otsu thresholding.
%        MinThreshold : numeric array (default: [0, 0]). Minimum threshold values for each channel.
%                       If empty, thresholds are determined using Otsu's method.
%        MaxThreshold : numeric array (default: []). Maximum threshold values per channel.
%                       If empty, thresholds are determined using Otsu's method.
%       BaseThreshold : numeric array (default: []). Base threshold values.
%                       If empty, thresholds are determined using Otsu's method.
%            volThrsh : numeric (default: 1000). Volume threshold to remove small isolated objects.
%              offset : numeric (default: 0). Offset value added to cleaned image to avoid zero background.
%        localWinSize : 1x3 vector (default: [15, 15, 23]). Local window size for regional cropping (y, x, z).
%          SigmaThrsh : numeric (default: 4). Sigma threshold for point detection.
%            intThrsh : numeric (default: 10000). Intensity threshold for removing bright peaks.
%          initDetect : char (default: '2d'). Initial detection method. Options: '2d' (MIP), '3d' (volume).
%         detVolThrsh : numeric (default: 5000). Volume threshold for object detection.
%                bbox : vector or empty (default: []). Bounding box specifying region of interest.
%           batchSize : 1x3 vector (default: [1024, 1024, 1024]). Batch size for processing (y, x, z).
%           blockSize : 1x3 vector (default: [256, 256, 256]). Block size for Zarr output (y, x, z).
%           Save16bit : true|false (default: false). Save output as 16-bit.
%              suffix : char (default: ''). Optional suffix for output directory.
%                uuid : char (default: ''). UUID string to label intermediate or final results.
%               debug : true|false (default: false). Enable debug mode.
%        parseCluster : true|false (default: true). Use SLURM cluster for parallel processing.
%       masterCompute : true|false (default: true). Whether the master job node also runs computations.
%         cpusPerTask : numeric (default: 1). Number of CPU cores allocated per task.
%             mccMode : true|false (default: false). Enable compiled MATLAB (MCC) mode.
%          configFile : char (default: ''). Path to config file for job submission.
%
% 
% Author: Xiongtao Ruan (02/18/2022)
% xruan (02/22/2022) add support for upper bound, and weighted threshold.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath', @(x) ischar(x) || iscell(x));
ip.addParameter('mode', 'updated', @(x) ischar(x) || iscell(x)); % updated, bg_interpolated, bg_estimation
ip.addParameter('bgFullpath', '', @(x) ischar(x)); % user provided background path
ip.addParameter('bgFactor', 1.0, @(x) isnumeric(x)); % scaling factor for user provided background
ip.addParameter('Sigma', 2.5, @isnumeric); % for 3D Gauss filtering
ip.addParameter('OTSUMaxPer', 99.9, @isnumeric); % Max percentile of data for OTSU calculation
ip.addParameter('MinThreshold', [0,0], @isnumeric); % if left empty, will use OTSU to calculate intensity to threshold; [ch1, ch2]
ip.addParameter('MaxThreshold', [], @isnumeric); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('BaseThreshold', [], @isnumeric); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('volThrsh', 1000, @isnumeric); % volume threshold for removing small isolated objects
ip.addParameter('offset', 0, @isnumeric); % offset to add to the cleaned image to make the background non-zero
ip.addParameter('localWinSize', [15, 15, 23], @isvector); % local window size for cropping local region
ip.addParameter('SigmaThrsh', 4, @isnumeric); % sigma threshold for point detection
ip.addParameter('intThrsh', 10000, @isnumeric); % intensity threshold for the peak to be removed
ip.addParameter('initDetect', '2d', @ischar); % initial detection method, 2d mip or 3d stack
ip.addParameter('detVolThrsh', 5000, @isnumeric); % volume threshold for the detetion
ip.addParameter('bbox', [], @(x) isempty(x) || isvector(x)); % bbox for the region to be processed. 
ip.addParameter('batchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256] , @isvector); % in y, x, z
ip.addParameter('Save16bit', false , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(zarrFullpath, varargin{:});
pr = ip.Results;
mode = pr.mode;
bgFullpath = pr.bgFullpath;
bgFactor = pr.bgFactor;
Sigma = pr.Sigma;
OTSUMaxPer = pr.OTSUMaxPer;
MinThreshold = pr.MinThreshold;
MaxThreshold = pr.MaxThreshold;
BaseThreshold = pr.BaseThreshold;
volThrsh = pr.volThrsh;
offset = pr.offset;
localWinSize = pr.localWinSize;
SigmaThrsh = pr.SigmaThrsh;
intThrsh = pr.intThrsh;
initDetect = pr.initDetect;
detVolThrsh = pr.detVolThrsh;
bbox = pr.bbox;
BatchSize = pr.batchSize;
BlockSize = pr.blockSize;
Save16bit = pr.Save16bit;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;

parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
configFile = pr.configFile;

bgEstimate = strcmp(mode, 'bg_estimation');

if ~exist(zarrFullpath, 'dir')
    error('The input zarr file %s does not exist!', zarrFullpath);
end

[dataPath, fsname, ext] = fileparts(zarrFullpath);
if bgEstimate
    outRtPath = [dataPath, '/Puncta_background_estimation/'];
    mkdir(outRtPath);
    outPath = sprintf('%s/%s_batch_%d_%d_%d/', outRtPath, fsname, BatchSize(1), BatchSize(2), BatchSize(3));
    outFullpath = sprintf('%s/%s_batch_%d_%d_%d_background_values.zarr', outRtPath, fsname, BatchSize(1), BatchSize(2), BatchSize(3));
else
    outPath = [dataPath, '/Puncta_Removed/'];
    outFullpath = sprintf('%s/%s.zarr', outPath, fsname);
    
    fprintf('Start Large-file puncta removal for %s...\n', fsname);
end
if exist(outFullpath, 'dir')
    disp('The output result exists, skip it!');
    return;
end
outTmppath = [outFullpath, '_', uuid];

if ~exist(outPath, 'dir')
    mkdir(outPath);
end

% save parameters
save('-v7.3', [outPath, '/parameters.mat'], 'pr');
s = jsonencode(pr, PrettyPrint=true);
fid = fopen([outPath, '/parameters.json'], 'w');
fprintf(fid, s);
fclose(fid);

fprintf(['reading ' fsname '...\n'])

if Save16bit
    dtype = 'uint16';
else
    dtype = 'single';
end

if parseCluster
    [parseCluster] = checkSlurmCluster(dataPath);
end

tic

% map input and output for xz
imSize = getImageSize(zarrFullpath);

SameBatchSize = false;
BorderSize = [100, 100, 100];
if bgEstimate
    BorderSize = min(round(BatchSize / 2), BorderSize);
end
% batchBBoxes: input batch bbox
% regionBBoxes: output batch bbox
% localBBoxes: bbox for the region from processed input to output batch.
[batchBBoxes, regionBBoxes, localBBoxes] = XR_zarrChunkCoordinatesExtraction(imSize, 'BatchSize', BatchSize, ...
    'BlockSize', BlockSize, 'SameBatchSize', SameBatchSize, 'BorderSize', BorderSize, 'bbox', bbox);

if isempty(bbox)
    outSize = imSize;
else
    outSize = bbox(4 : 6) - bbox(1 : 3) + 1;
end    

% check if intermediate result exists
fresh_run = true & ~bgEstimate;
if ~bgEstimate && exist(outTmppath, 'dir')
    bim = blockedImage(outTmppath, 'Adapter', CZarrAdapter);
    out_tmp_size = bim.Size;
    out_tmp_block_size = bim.BlockSize;
    if numel(outSize) == numel(out_tmp_size) && all(outSize == out_tmp_size) && ...
            all(BlockSize == out_tmp_block_size)
        fresh_run = false;
    else
        rmdir(outTmppath, 's');
    end
end

% do not overwrite the intermediate data if it exists
if fresh_run
    % initialize zarr file
    dimSeparator = '.';
    if prod(ceil(outSize ./ BlockSize)) > 10000
        dimSeparator = '/';
    end
    createzarr(outTmppath, dataSize=outSize, BlockSize=BlockSize, dtype=dtype, dimSeparator=dimSeparator);
end

if ~bgEstimate
    zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', outPath, fsname, uuid);
    if fresh_run && exist(zarrFlagPath, 'dir')
        rmdir(zarrFlagPath, 's');
    end
    mkdir_recursive(zarrFlagPath);
end

% set up parallel computing 
taskSize = 5; % the number of batches a job should process
% taskSize = min(200, taskSize * ceil(prod(1024 ./ BatchSize)));
taskSize = taskSize * ceil(prod(1024 ./ BatchSize));
numBatch = size(batchBBoxes, 1);
numTasks = ceil(numBatch / taskSize);
if numTasks > 20000
    taskSize = min(100000, ceil(numBatch / 20000));
    numTasks = ceil(numBatch / taskSize);
end

taskBatchNum = 1;

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    batchBBoxes_i = batchBBoxes(batchInds, :);
    
    if bgEstimate
        outputFullpath = sprintf('%s/blocks_%d_%d.mat', outPath, batchInds(1), batchInds(end));
        outputFullpaths{i} = outputFullpath;
        Overwrite = false;
        
        funcStrs{i} = sprintf(['XR_ExM_PunctaRemoval_background_estimation_block([%s],''%s'',''%s'',%s,', ...
            '''Overwrite'',%s,''Sigma'',%0.20d,''OTSUMaxPer'',%0.20d,''MinThreshold'',%0.20d,', ...
            '''MaxThreshold'',[%0.20d],''BaseThreshold'',[%0.20d],''uuid'',''%s'',''debug'',%s)'], ...
            strrep(num2str(batchInds, '%d,'), ' ', ''), zarrFullpath, outputFullpath, ...
            strrep(mat2str(batchBBoxes_i), ' ', ','), string(Overwrite), Sigma, OTSUMaxPer, ...
            MinThreshold, MaxThreshold, BaseThreshold, uuid, string(debug));
    else
        regionBBoxes_i = regionBBoxes(batchInds, :);
        localBBoxes_i = localBBoxes(batchInds, :);        
        zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
        outputFullpaths{i} = zarrFlagFullpath;
        Overwrite = false;
        switch mode
            case 'bg_interpolated'
                funcStrs{i} = sprintf(['XR_ExM_PunctaRemoval_interpolated_background_block([%s],', ...
                    '''%s'',''%s'',''%s'',''%s'',%s,%s,%s,''Overwrite'',%s,''bgFactor'',%0.20d,''Sigma'',%0.20d,', ...
                    '''OTSUMaxPer'',%0.20d,''volThrsh'',%0.20d,''offset'',%d,''localWinSize'',%s,', ...
                    '''SigmaThrsh'',%s,''intThrsh'',%s,''initDetect'',''%s'',''detVolThrsh'',%s,', ...
                    '''uuid'',''%s'',''debug'',%s)'], strrep(num2str(batchInds, '%d,'), ' ', ''), ...
                    zarrFullpath, outTmppath, bgFullpath, zarrFlagFullpath, strrep(mat2str(batchBBoxes_i), ' ', ','), ...
                    strrep(mat2str(regionBBoxes_i), ' ', ','), strrep(mat2str(localBBoxes_i), ' ', ','), ...
                    string(Overwrite), bgFactor, Sigma, OTSUMaxPer, volThrsh, offset, ...
                    strrep(mat2str(localWinSize), ' ', ','), strrep(mat2str(SigmaThrsh), ' ', ','), ...
                    strrep(mat2str(intThrsh), ' ', ','), initDetect, strrep(mat2str(detVolThrsh), ' ', ','), ...
                    string(debug));
            case 'updated'
                funcStrs{i} = sprintf(['XR_ExM_PunctaRemoval_updated_block([%s],''%s'',''%s'',''%s'',%s,%s,', ...
                    '%s,''Overwrite'',%s,''Sigma'',%0.20d,''OTSUMaxPer'',%0.20d,''MinThreshold'',%0.20d,', ...
                    '''MaxThreshold'',[%0.20d],''BaseThreshold'',[%0.20d],''volThrsh'',%0.20d,''offset'',%d,', ...
                    '''localWinSize'',%s,''SigmaThrsh'',%s,''intThrsh'',%s,''initDetect'',''%s'',''detVolThrsh'',%s,', ...
                    '''uuid'',''%s'',''debug'',%s)'], strrep(num2str(batchInds, '%d,'), ' ', ''), ...
                    zarrFullpath, outTmppath, zarrFlagFullpath, strrep(mat2str(batchBBoxes_i), ' ', ','), ...
                    strrep(mat2str(regionBBoxes_i), ' ', ','), strrep(mat2str(localBBoxes_i), ' ', ','), ...
                    string(Overwrite), Sigma, OTSUMaxPer, MinThreshold, MaxThreshold, BaseThreshold, ...
                    volThrsh, offset, strrep(mat2str(localWinSize), ' ', ','), strrep(mat2str(SigmaThrsh), ' ', ','), ...
                    strrep(mat2str(intThrsh), ' ', ','), initDetect, strrep(mat2str(detVolThrsh), ' ', ','), ...
                    uuid, string(debug));
        end
    end
end

inputFullpaths = repmat({zarrFullpath}, numTasks, 1);

memAllocate = prod(BatchSize + BorderSize * 2) * 4 / 1024^3 * 5;
% submit jobs
is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, funcStrs, ...
    finalOutFullpath=outFullpath, parseCluster=parseCluster, masterCompute=masterCompute, taskBatchNum=taskBatchNum, ...
    cpusPerTask=cpusPerTask, memAllocate=memAllocate, mccMode=mccMode, ConfigFile=configFile);


if ~all(is_done_flag)
    error('The computing for some batchs is missing!');
end

if bgEstimate
    % collect results
    fprintf('Collect local background values...\n')
    bgSize = ceil(outSize ./ BatchSize);
    bg_mat = zeros(bgSize, 'single');
    for i = 1 : numTasks
        batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
        regionBBoxes_i = regionBBoxes(batchInds, :);

        blockFn = sprintf('%s/blocks_%d_%d.mat', outPath, batchInds(1), batchInds(end));
        a = load(blockFn, 'bg_mat');

        subs_i = (regionBBoxes_i(:, 1 : 3) - 1) ./ BatchSize + 1;
        inds_i = sub2ind(bgSize, subs_i(:, 1), subs_i(:, 2), subs_i(:, 3));
        bg_mat(inds_i) = a.bg_mat;
    end

    BlockSize = ceil(1024 ./ BatchSize);
    dtype = 'single';
    dimSeparator = '.';
    createzarr(outTmppath, dataSize=bgSize, BlockSize=BlockSize, dtype=dtype, dimSeparator=dimSeparator);
    writezarr(bg_mat, outTmppath);
    fprintf('Done!\n')    
end

if exist(outFullpath, 'dir') && exist(outTmppath, 'dir')
    rmdir(outFullpath, 's');
end
if exist(outTmppath, 'dir')
    movefile(outTmppath, outFullpath);
end
    
if ~bgEstimate
    if exist(zarrFlagPath, 'dir')
        rmdir(zarrFlagPath, 's');
    end    
    % generate MIP z file
    XR_MIP_zarr(outFullpath, parseCluster=parseCluster, mccMode=mccMode, ConfigFile=ConfigFile);
end
toc

end

