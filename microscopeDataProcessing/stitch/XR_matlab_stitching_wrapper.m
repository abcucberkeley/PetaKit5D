function [] = XR_matlab_stitching_wrapper(dataPath, imageListFileName, varargin)
% wrapper for java stitching pipeline. 
% 
% 
% Inputs :   
%        imageDirName : full path for the directory of the images to be stitched
%   imageListFileName : full path for the coordinate information csv file
%
% Options (as 'specifier'-value pairs): 
%
%         'axisOrder' : Axis mapping of the coordinate system. Default: 'xyz'.
%         'ffcorrect' : flat field correction. Default: false.
%        'Resolution' : Image resolution in um, 2D (same for xy) or 3D vector.  
%         'resultDir' : name of stitching result directory. The directory is child directory of imageDirname
%       'blendMethod' : Blend method for overlap regions. Available: none, mean, max, median. Default: mean
%           'padSize' : Pad or crop the stitched image, empty (default) or 
%                       a 1X3 vector of integers (y, x, z). Postive for pad and negative for crop. 
%      'outBbox' : Crop the stitched image by some bounding box, empty (default, no crop) 
%                       or a 3X2 vector for start and end indices of the bounding box (y, x, z).      
%        'zNormalize' : normalize background along z-axis by the median. Default: false.
%      'parseCluster' : Use slurm-based cluster computing. Default: true. 
%         'jobLogDir' : Log directory for the slurm jobs. 
%       'cpusPerTask' : Number of cpus for a job. Default: 12
%         'save16bit' : true|{false}. Save final results as 16bit or single. 
%              'uuid' : unique string for a job for saving files. 
%       'maxTrialNum' : Max number of times to rerun failure cases. 
%      'unitWaitTime' : For computing without cluster, the wait time per      
%                       file in minutes, in order to check whether the computing is done. 
%
%
% Author: Xiongtao Ruan (02/18/2020)
%
% xruan: add xcorr based stitching
% xruan (04/06/2020): add option for computing of only first time point and
% save parameter mat files for the running. 
% xruan (06/19/2020): add option for applying shifts from primary channel for other channels. 
% xruan (06/19/2020): add option for applying shifts from primary channel of the first time point. 
% xruan (07/11/2020): For primary first option, if masterCPU is true, use
%                     master script to do the computing to save computing resources. 
% xruan (07/11/2020): add option for cpu only nodes
% xruan (07/13/2020): add option for streaming mode (where the data is coming 
%                     in real time). Also, check if the frame
% xruan (07/21/2020): add support for new name format (Tile_XXX before Scan)
% xruan (07/21/2020): fix issue for default primary channel in case it is not exist. 
% xruan (07/26/2020): add option to use existing stitch info to stitch for all images (including first time point) 
% xruan (08/01/2020): add support for file names without CamA/B
% xruan (08/02/2020): add support for increase cpuPerTasks if the result is
%                     expect to be large
% xruan (08/02/2020): add support for primary channel for no xcorrshift
% xruan (08/17/2020): add support for stitching of DSR decon (only for
% existing DSR decon files).
% xruan (08/20/2020): add support for objective scan
% xruan (08/23/2020): add option for overlap type (full)
% xruan (09/23/2020): add prefix for filename in case string before Scan
% xruan (10/05/2020): add zarr-based stitching pipeline as an option; use
%                     tile number in each dimension for pad size computing
% xruan (10/06/2020): filter out partial file records in image list
% xruan (10/18/2020): remove incomplete time point for primary or primaryfirst options
% xruan (10/23/2020): add support for chosen time points and allow empty
%                     folder for raw data if using DSR or DSR decon
% xruan (10/24/2020): add support for user-defined processing on tiles
% xruan (12/06/2020): add support for subIter/fullIter for Iter, and
%                     flipped tiles
% xruan (12/09/2020): add support for using primary coordinates for secondary channels/tps
% xruan (02/24/2021): add support for user defined xy, z max offsets for xcorr registration
% xruan (03/24/2021): add support for processing of given channels
% xruan (07/04/2021): ignore tif files recorded in Image List files that
%   do not match the file pattern. 
% xruan (07/05/2021): add support for user defined resample (arbitary factor)
% xruan (07/15/2021): extend subiteration to unlimited number of sets
% xruan (08/25/2021): add support for channel-specific user functions
% xruan (10/13/2021): add support for cropping tiles; add support for
%                     skewed space stitching with reference (for decon data)
% xruan (10/28/2021): add support for IO scan
% xruan (12/17/2021): add support for 2d stitching (e.g., MIPs), only with dsr for now
% xruan (01/03/2022): add support for online stitch (partial stitch) for each layer
% xruan (01/25/2022): add support for loading tileFullpaths and coordinates from file.
% xruan (01/31/2022): add support for negative tile numbers in xyz. 
% xruan (06/20/2022): add input variable axisWeight for user defined weights for optimization
% xruan (08/25/2022): change CropToSize to tileOutBbox (more generic)
% xruan (12/13/2022): change xcorr thresh as user defined parameter with default 0.25
% xruan (02/15/2023): add support for multi-location stitching (tiles are
%   in different subfolders and/or subgroups). 
% xruan (04/29/2023): add support for zarr files directly as inputs, 
%                     and also for stitching of large zarr files


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @(x) ischar(x) || iscell(x));
ip.addRequired('imageListFileName', @(x) ischar(x) || iscell(x));
% ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('streaming', false, @islogical);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('multiLoc', false, @islogical); % use subregions from different folders
ip.addParameter('processedDirStr', '', @ischar); % path for using existing exist processed data (i.e., DSR, decon)
ip.addParameter('stitchInfoFullpath', '', @ischar); % use exist stitch info for stitching
ip.addParameter('DS', false, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.5, @isnumeric);
ip.addParameter('skewAngle', 32.45, @isnumeric);
ip.addParameter('reverse', false, @islogical);
ip.addParameter('parseSettingFile', false, @islogical); % use setting file to decide whether filp Z stack or not.
ip.addParameter('axisOrder', 'x,y,z', @ischar); % stage coordinates axis order
ip.addParameter('dataOrder', 'y,x,z', @ischar); % data axis order, in case data in zyx order
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('IOScan', false, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('largeZarr', false, @islogical);
ip.addParameter('poolSize', [], @isnumeric);
ip.addParameter('batchSize', [500, 500, 500], @isnumeric);
ip.addParameter('blockSize', [500, 500, 500], @isnumeric);
ip.addParameter('shardSize', [], @isnumeric);
ip.addParameter('resampleType', 'xy_isotropic', @ischar); % by default use xy isotropic
ip.addParameter('resampleFactor', [], @isnumeric); % user-defined resample factor
ip.addParameter('inputBbox', [], @isnumeric); % crop input tile before processing
ip.addParameter('tileOutBbox', [], @isnumeric); % crop tile after processing 
ip.addParameter('tileOffset', 0, @isnumeric); % offset added to tile
ip.addParameter('resultDirName', 'matlab_stitch', @ischar);
ip.addParameter('blendMethod', 'none', @ischar);
ip.addParameter('overlapType', '', @ischar); % '', 'none', 'half', or 'full'
ip.addParameter('xcorrShift', true, @islogical);
ip.addParameter('xyMaxOffset', 300, @(x) isvector(x) && numel(x) <= 2); % max offsets in xy axes
ip.addParameter('zMaxOffset', 50, @isnumeric); % max offsets in z axis
ip.addParameter('xcorrDownsample', [2, 2, 1], @isnumeric); % max offsets in z axis
ip.addParameter('xcorrThresh', 0.25, @isnumeric); % threshold of of xcorr, ignore shift if xcorr below this threshold.
ip.addParameter('padSize', [], @(x) isnumeric(x) && (isempty(x) || numel(x) == 3));
ip.addParameter('outBbox', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6));
ip.addParameter('xcorrMode', 'primaryFirst', @(x) ismember(x, {'primary', 'primaryFirst', 'all'}));  % 'primary': choose one channel as primary channel, 
        % 'all': xcorr shift for each channel;  % 'primaryFirst': the primary channel of first time point
ip.addParameter('shiftMethod', 'grid', @ischar); % {'local', 'global', 'grid', 'group', 'test'}
ip.addParameter('axisWeight', [1, 0.1, 10], @isnumeric); % axis weight for optimization, y, x, z
ip.addParameter('groupFile', '', @ischar); % file to define tile groups
ip.addParameter('primaryCh', '', @(x) isempty(x) || ischar(x)); % format: CamA_ch0. If it is empty, use the first channel as primary channel
ip.addParameter('usePrimaryCoords', false, @islogical); 
ip.addParameter('save16bit', false, @islogical);
ip.addParameter('edgeArtifacts', 0, @isnumeric);
ip.addParameter('distBboxes', [], @isnumeric); % bounding boxes for distance transform
ip.addParameter('saveMIP', true, @islogical);
ip.addParameter('stitchMIP', [], @(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3))); % 1x3 vector or vector, by default, stitch MIP-z
ip.addParameter('onlineStitch', false, @(x) islogical(x)); % support for online stitch (with partial number of tiles). 
ip.addParameter('processFunPath', '', @(x) isempty(x) || ischar(x) || iscell(x)); % path of user-defined process function handle
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 8, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 0.1, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPath, imageListFileName, varargin{:});

pr = ip.Results;
% Overwrite = pr.Overwrite;
streaming = pr.streaming;
channelPatterns = pr.channelPatterns;
multiLoc = pr.multiLoc;
processedDirStr = pr.processedDirStr;
stitchInfoFullpath = pr.stitchInfoFullpath;
DS = pr.DS;
DSR = pr.DSR;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
skewAngle = pr.skewAngle;
reverse = pr.reverse;
parseSettingFile = pr.parseSettingFile;
axisOrder = pr.axisOrder;
dataOrder = pr.dataOrder;
objectiveScan = pr.objectiveScan;
IOScan =  pr.IOScan;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
poolSize = pr.poolSize;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
shardSize = pr.shardSize;
resampleType = pr.resampleType;
resampleFactor = pr.resampleFactor;
inputBbox = pr.inputBbox;
tileOutBbox = pr.tileOutBbox;
tileOffset = pr.tileOffset;
resultDirName = pr.resultDirName;
blendMethod = pr.blendMethod;
overlapType = pr.overlapType;
xcorrShift = pr.xcorrShift;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
xcorrDownsample = pr.xcorrDownsample;
xcorrThresh = pr.xcorrThresh;
padSize = pr.padSize;
outBbox = pr.outBbox;
xcorrMode = pr.xcorrMode;
shiftMethod = pr.shiftMethod;
axisWeight = pr.axisWeight;
groupFile = pr.groupFile;
primaryCh = pr.primaryCh;
usePrimaryCoords = pr.usePrimaryCoords;
save16bit = pr.save16bit;
edgeArtifacts = pr.edgeArtifacts;
distBboxes = pr.distBboxes;
saveMIP = pr.saveMIP;
stitchMIP = pr.stitchMIP;
onlineStitch = pr.onlineStitch;
processFunPath = pr.processFunPath;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
mccMode = pr.mccMode;
configFile = pr.configFile;

% make root directory
if multiLoc
    stitching_rt = [dataPath{1}, '/', resultDirName];
else
    if iscell(dataPath)
        dataPath = dataPath{1};
    end
    if iscell(imageListFileName)
        imageListFileName = imageListFileName{1};
    end
    stitching_rt = [dataPath, '/', resultDirName];
end
if ~exist(stitching_rt, 'dir')
    mkdir(stitching_rt);
    fileattrib(stitching_rt, '+w', 'g');
end

% save parameters 
save('-v7.3', [stitching_rt, '/parameters.mat'], 'pr');
writetable(struct2table(pr, 'AsArray', true), [stitching_rt, '/parameters.txt'])

% check if axis order is valid
axisOrder = strrep(axisOrder, ' ', '');
pattern = '^(-?x,?-?y,?-?z|-?y,?-?x,?-?z|-?z,?-?y,?-?x|-?x,?-?z,?-?y|-?x,?-?z,?-?y|-?y,?-?z,?-?x)$';
if ~regexpi(axisOrder, pattern)
    error("The axisOrder is not right, it must has the form like 'y,x,z' or '-x,y,z' (flipped in x-axis)!");
end

% % save xcorr info
stitchInfoDir = 'stitchInfo';
stitch_info_path = [stitching_rt, '/', stitchInfoDir];
if ~exist(stitch_info_path, 'dir')
    mkdir(stitch_info_path);
    fileattrib(stitch_info_path, '+w', 'g');            
end

% temporary directory for intermediate results
stitching_tmp = [stitching_rt, '/', 'tmp'];
if ~exist(stitching_tmp, 'dir')
    mkdir(stitching_tmp);
    fileattrib(stitching_tmp, '+w', 'g');            
end

% check if a slurm-based computing cluster exists
if parseCluster
    if multiLoc
        [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath{1}, jobLogDir);
    else
        [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath, jobLogDir);
    end
end

% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end

if any(stitchMIP)
    cpusPerTask = 1;
end

%% parse image list information

useProcessedData = ~isempty(processedDirStr);

if multiLoc
    [tab, primary_tab, fullIter, Ch, Cam, stackn, nz, specifyCam, prefix, zlayerStitch, xcorrMode, stitchInfoFullpath] = ...
        stitch_parse_multi_loc_image_list_information(dataPath, imageListFileName, streaming=streaming, ...
        onlineStitch=onlineStitch, stitchInfoFullpath=stitchInfoFullpath, stitchInfoPath=stitch_info_path, ...
        zarrFile=zarrFile, channelPatterns=channelPatterns, useProcessedData=useProcessedData, ...
        processedDirStr=processedDirStr, xcorrMode=xcorrMode, primaryCh=primaryCh);
else
    [tab, primary_tab, fullIter, Ch, Cam, stackn, nz, specifyCam, prefix, zlayerStitch, xcorrMode, stitchInfoFullpath] = ...
        stitch_parse_image_list_information(dataPath, imageListFileName, streaming=streaming, ...
        onlineStitch=onlineStitch, stitchInfoFullpath=stitchInfoFullpath, stitchInfoPath=stitch_info_path, ...
        zarrFile=zarrFile, channelPatterns=channelPatterns, useProcessedData=useProcessedData, ...
        processedDirStr=processedDirStr, xcorrMode=xcorrMode, primaryCh=primaryCh);
end

if ischar(processFunPath)
    processFunPath = repmat({processFunPath}, numel(channelPatterns), 1);
elseif iscell(processFunPath)
    for i = 1 : numel(processFunPath)
        if isempty(processFunPath{i})
            processFunPath{i} = '';
        end
    end
end

%% do stitching computing
row_exist_flag = true(numel(fullIter), numel(Cam), numel(stackn), numel(Ch), numel(nz)); % flag for whether the run exists.
is_done_flag = false(numel(fullIter), numel(Cam), numel(stackn), numel(Ch), numel(nz));
trial_counter = zeros(numel(fullIter), numel(Cam), numel(stackn), numel(Ch), numel(nz));
max_trial_num = maxTrialNum;

if parseCluster
    job_ids = -ones(numel(fullIter), numel(Cam), numel(stackn), numel(Ch), numel(nz));
    job_status_flag = false(numel(fullIter), numel(Cam), numel(stackn), numel(Ch), numel(nz));
end

% set wait counter for streaming option
if streaming
    stream_counter = 0;
    if onlineStitch
        stream_max_counter = 5;
    else
        stream_max_counter = 100 * size(tab, 1);
    end
end

if onlineStitch
    tileNums = zeros(numel(fullIter), numel(Cam), numel(stackn), numel(Ch), numel(nz));
end

bigStitchData = false;
if ~onlineStitch && ~any(stitchMIP) && ~streaming && size(tab, 1) > 20
    tile_fnames = tab.Filename;
    if useProcessedData
        if multiLoc
            tile_fullpath_1 =  [dataPath{tab.did(1)}, '/', processedDirStr, '/', tile_fnames{1}];
        else
            tile_fullpath_1 =  [dataPath, '/', processedDirStr, '/', tile_fnames{1}];
        end
    else
        if multiLoc
            tile_fullpath_1 =  [dataPath{tab.did(1)}, '/', tile_fnames{1}];
        else
            tile_fullpath_1 =  [dataPath, '/', tile_fnames{1}];
        end
    end
    sz = getImageSize(tile_fullpath_1);
    if prod(sz) * 4 * numel(tile_fnames) > 100 * 1024^3
        bigStitchData = true;
    end
end

while ~all(is_done_flag | trial_counter >= max_trial_num, 'all')
    % exit the job if no new images are transferred.
    if streaming && stream_counter > stream_max_counter
        break;
    end
    
    lastF = find(~is_done_flag & trial_counter < maxTrialNum, 1, 'last');
    for n = 1:numel(fullIter)
        for ncam = 1:numel(Cam)
            for s = 1:numel(stackn)
                for c = 1:numel(Ch)
                    for z = 1:numel(nz)
                        if is_done_flag(n, ncam, s, c, z) || trial_counter(n, ncam, s, c, z) >= max_trial_num
                            continue;
                        end
                        if ~row_exist_flag(n, ncam, s, c, z)
                            is_done_flag(n, ncam, s, c, z) = true;
                            continue;
                        end

                        f = sub2ind([numel(fullIter), numel(Cam), numel(stackn), numel(Ch), numel(nz)], n, ncam, s, c, z);
                        if zlayerStitch
                            cur_tab = tab(tab.ch == Ch(c) & tab.camera == Cam(ncam) & strcmp(tab.fullIter, fullIter{n}) & tab.stack == stackn(s) & tab.z == nz(z), :);
                        else
                            cur_tab = tab(tab.ch == Ch(c) & tab.camera == Cam(ncam) & strcmp(tab.fullIter, fullIter{n}) & tab.stack == stackn(s), :);                            
                        end
                        
                        % obtain filenames                    
                        if isempty(cur_tab)
                            row_exist_flag(n, ncam, s, c, z) = false;
                            is_done_flag(n, ncam, s, c, z) = true;

                            % if the primary channel is missing, then skip
                            % other channels if it is primary mode
                            if (ncam == 1 && c == 1) && strcmp(xcorrMode, 'primary')
                                row_exist_flag(1, 1 : end, 1, 1 : end, 1) = false;
                                is_done_flag(1, 1 : end, 1, 1 : end, 1) = true;
                            end
                            continue;
                        end

                        % for primaryfirst option, check whether the number of
                        % tiles are the same as primary channel
                        if strcmpi(xcorrMode, 'primaryFirst') && size(cur_tab, 1) ~= size(primary_tab, 1) && ~onlineStitch
                            fprintf('The number of tiles is different from the number of tiles in the primary channel, skip it!\n');
                            disp(cur_tab.Filename);
                            is_done_flag(n, ncam, s, c, z) = true;
                            continue;
                        end
                        
                        laser = unique(cur_tab.laser);
                        abstime = unique(cur_tab.abstime);
                        % reltime = unique(cur_t.t);
                        fpgatime = cur_tab.fpgatime(1);
                        xyz = [cur_tab.StageX_um_, cur_tab.StageY_um_, cur_tab.StageZ_um_];
                        % also add data path id as part of the tileIdx
                        tileIdx = [cur_tab.x, cur_tab.y, cur_tab.z, cur_tab.did];

                        if numel(laser) > 1
                            laser = laser(1);
                        end

                        tile_fnames = cur_tab.Filename;
                        if multiLoc
                            nF = numel(tile_fnames);
                            tile_fullpaths = cell(nF, 1);
                            for f = 1 : nF
                                tile_fullpaths{f} = sprintf('%s/%s', dataPath{cur_tab.did(f)}, tile_fnames{f});
                            end
                        else
                            tile_fullpaths = cellfun(@(x) [dataPath, '/', x], tile_fnames, 'unif', 0);
                        end

                        % parse setting file 
                        flippedTile = [];
                        if parseSettingFile
                            settingInfo = XR_parseSettingFiles_wrapper(tile_fullpaths);
                            flippedTile = [settingInfo.StageInterval] < 0;
                        end

                        % check if files exist for streaming option, and also
                        % if useProcessedData is true, check the processed files exist.
                        if streaming
                            is_tile_exist = cellfun(@(x) exist(x, 'file'), tile_fullpaths);
                            if ~all(is_tile_exist)
                                stream_counter = stream_counter + 1;
                                continue;
                            else 
                                stream_counter = 0;
                            end
                        end
                        if useProcessedData
                            if multiLoc
                                nF = numel(tile_fnames);
                                tile_processed_fullpaths = cell(nF, 1);
                                for f = 1 : nF
                                    tile_processed_fullpaths{f} = sprintf('%s/%s/%s', dataPath{cur_tab.did(f)}, processedDirStr, tile_fnames{f});
                                end
                            else
                                tile_processed_fullpaths = cellfun(@(x) [dataPath, '/', processedDirStr, '/', x], tile_fnames, 'unif', 0);
                            end
                            is_processed_tile_exist = cellfun(@(x) exist(x, 'file'), tile_processed_fullpaths);
                            if streaming
                                if ~all(is_processed_tile_exist)
                                    stream_counter = stream_counter + 1;
                                    continue;
                                else
                                    stream_counter = 0;
                                end
                            else
                                if ~all(is_processed_tile_exist) 
                                    % is_done_flag(n, ncam, s, c) = true;
                                    trial_counter(n, ncam, s, c, z) = trial_counter(n, ncam, s, c, z) + 1;
                                    continue; 
                                end
                            end
                        else
                            processedDirStr = '';
                        end
                        
                        if zlayerStitch
                            z_str = sprintf('_%03dz', nz(z));
                        else
                            z_str = '';
                        end
                        if onlineStitch
                             z_str = sprintf('%s_%03dntile', z_str, size(cur_tab, 1));
                        end

                        % for secondary channels, first check whether the
                        % stitching info file exist
                        isPrimaryCh = true;
                        % the stitch Info full path for options except 'primaryFirst'
                        switch xcorrMode
                            case 'primaryFirst'
                                if ~(n == 1 && ncam == 1 && s == 1 && c == 1 && z == 1)
                                    isPrimaryCh = false;
                                end 
                            case {'stitchInfo'}
                                isPrimaryCh = false;
                            case {'primary', 'all'}
                                if specifyCam
                                    stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s.mat', ...
                                        stitch_info_path, prefix, fullIter{n}, Cam(ncam), Ch(c), stackn(s), laser, abstime, fpgatime, z_str);
                                else
                                    stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s.mat', ...
                                        stitch_info_path, prefix, fullIter{n}, Ch(c), stackn(s), laser, abstime, fpgatime, z_str);                            
                                end
                        end

                        if true || xcorrShift 
                            % get the primary info path for secondary channels for 'primary' option
                            if strcmpi(xcorrMode, 'primary') && ~(ncam == 1 && c == 1)
                                isPrimaryCh = false;
                                primary_tab = tab(tab.ch == Ch(1) & tab.camera == Cam(1) & strcmp(tab.fullIter, fullIter{n}) & tab.stack == stackn(s), :);
                                p_laser = unique(primary_tab.laser);
                                p_abstime = unique(primary_tab.abstime);
                                p_fpgatime = primary_tab.fpgatime(1);
                                p_z = unique(primary_tab.z);
                                p_z = p_z(1);

                                if numel(p_laser) > 1
                                    p_laser = p_laser(1);
                                end
                                if zlayerStitch
                                    pz_str = sprintf('_%03dz', p_z);
                                else
                                    pz_str = '';
                                end
                                if specifyCam
                                    stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s.mat', ...
                                        stitch_info_path, prefix, fullIter{n}, Cam(1), Ch(1), stackn(s), p_laser, p_abstime, p_fpgatime, pz_str);
                                else
                                    stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s.mat', ...
                                        stitch_info_path, prefix, fullIter{n}, Ch(1), stackn(1), p_laser, p_abstime, p_fpgatime, pz_str);        
                                end
                            end

                            % for secondary channels, if the stitchInfo file  not exist, wait the stitching for the primary channel
                            if ~isPrimaryCh && ~exist(stitchInfoFullpath, 'file')
                                continue;
                            end
                        end

                        % also use flag based check of completion, to support distributed computing with same submission
                        if specifyCam
                            stitch_save_fsname = sprintf('%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s', ...
                                stitching_rt, prefix, fullIter{n}, Cam(ncam), Ch(c), stackn(s), laser, abstime, fpgatime, z_str);
                            cur_tmp_fname = sprintf('%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s.tmp', ...
                                stitching_rt, prefix, fullIter{n}, Cam(ncam), Ch(c), stackn(s), laser, abstime, fpgatime, z_str);
                        else
                            stitch_save_fsname = sprintf('%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s', ...
                                stitching_rt, prefix, fullIter{n}, Ch(c), stackn(s), laser, abstime, fpgatime, z_str);
                            cur_tmp_fname = sprintf('%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s.tmp', ...
                                stitching_rt, prefix, fullIter{n}, Ch(c), stackn(s), laser, abstime, fpgatime, z_str);
                        end

                        stitch_save_fname = [stitch_save_fsname, '.zarr'];
                        ftype = 'dir';
                        [~, fsname] = fileparts(stitch_save_fname);

                        % for online stitch check if old results with fewer tiles exist, if so, delete them. 
                        if onlineStitch
                            % xruan: for bidirectional scan, the fpga time may change
                            % dir_info = dir([stitch_save_fsname(1:end-8), '*ntile.zarr']);
                            dir_info = dir([stitch_save_fsname(1 : end - 32), '_*msecAbs_', z_str(2:5), '*ntile.zarr']);
                            stitch_save_files = {dir_info.name}';

                            stitch_tile_nums = regexp(stitch_save_files, '_(\d+)ntile', 'tokens');
                            stitch_tile_nums = cellfun(@(x) str2double(x{1}{1}), stitch_tile_nums);
                            
                            for st = 1 : numel(stitch_tile_nums)
                                cur_stitch_save_fname = [stitching_rt, '/', stitch_save_files{st}];
                                [~, cur_stitch_save_fsname] = fileparts(cur_stitch_save_fname);
                                cur_stitch_save_mip_fname = sprintf('%s/MIPs/%s_MIP_z.tif', stitching_rt, cur_stitch_save_fsname);
                                if stitch_tile_nums(st) < size(cur_tab, 1) && exist(cur_stitch_save_fname, ftype)
                                    switch ftype
                                        case 'file'
                                            delete(cur_stitch_save_fname); 
                                        case 'dir'
                                            rmdir(cur_stitch_save_fname, 's');
                                    end
                                    delete(cur_stitch_save_mip_fname);
                                    is_done_flag(n, ncam, s, c, z) = false;
                                    tileNums(f) = size(cur_tab, 1);
                                end
                            end
                        end
                            
                        % first check if the computing is done
                        if exist(stitch_save_fname, ftype)
                            is_done_flag(n, ncam, s, c, z) = true; 
                            if exist(cur_tmp_fname, 'file')
                                delete(cur_tmp_fname);
                            end
                            continue;
                        end

                        tile_fullpaths_str = sprintf('{''%s''}', strjoin(tile_fullpaths, ''','''));
                        xyz_str = strrep(mat2str(xyz), ' ', ',');  
                        tileIdx_str = strrep(mat2str(tileIdx), ' ', ',');
                        tileInfoFullpath = '';
                        flippedTile_str = strrep(num2str(flippedTile, '%d,'), ' ', '');

                        cind = cellfun(@(x) contains(tile_fullpaths{1}, x), channelPatterns);                        
                        processFunPath_str = sprintf('{''%s''}', strjoin(processFunPath(cind, :), ''','''));

                        % for tile number greater than 10, save the info to the disk and load it for the function
                        if numel(tile_fullpaths) > 10 && parseCluster && job_ids(n, ncam, s, c, z) == -1
                            fprintf('Save tile paths and coordinates to disk...\n');
                            tileInfoFullpath = sprintf('%s/stitchInfo/%s_tile_info.mat', stitching_rt, fsname);
                            tileInfoTmppath = sprintf('%s/stitchInfo/%s_tile_info_%s.mat', stitching_rt, fsname, uuid);
                            save('-v7.3', tileInfoTmppath, 'tile_fullpaths', 'xyz', 'tileIdx', 'flippedTile');
                            movefile(tileInfoTmppath, tileInfoFullpath);
                            
                            tile_fullpaths_str = '{}';
                            xyz_str = '[]';
                            tileIdx_str = '[]';
                            flippedTile_str = '';
                        end
                        
                        func_str = sprintf(['XR_stitching_frame_zarr_dev_v1(%s,%s,''axisOrder'',''%s'',''dataOrder'',''%s'',', ...
                            '''xyPixelSize'',%0.10f,''dz'',%0.10f,''skewAngle'',%0.10f,''reverse'',%s,''objectiveScan'',%s,', ...
                            '''IOScan'',%s,''resultPath'',''%s'',''tileInfoFullpath'',''%s'',''stitchInfoDir'',''%s'',', ...
                            '''stitchInfoFullpath'',''%s'',''processedDirStr'',''%s'',''DS'',%s,''DSR'',%s,''zarrFile'',%s,', ...
                            '''largeZarr'',%s,''poolSize'',%s,''batchSize'',%s,''blockSize'',%s,''shardSize'',%s,''resampleType'',''%s'',', ...
                            '''resample'',[%s],''inputBbox'',%s,''tileOutBbox'',%s,''tileOffset'',%d,''blendMethod'',''%s'',', ...
                            '''overlapType'',''%s'',''xcorrShift'',%s,''xyMaxOffset'',%0.10f,''zMaxOffset'',%0.10f,''xcorrDownsample'',%s,', ...
                            '''xcorrThresh'',%0.10f,''shiftMethod'',''%s'',''axisWeight'',[%s],''groupFile'',''%s'',''isPrimaryCh'',%s,', ...
                            '''usePrimaryCoords'',%s,''padSize'',[%s],''outBbox'',[%s],''save16bit'',%s,', ...
                            '''tileIdx'',%s,''flippedTile'',[%s],''processFunPath'',%s,''stitchMIP'',%s,''bigStitchData'',%s,', ...
                            '''edgeArtifacts'',%0.10f,''distBboxes'',%s,''saveMIP'',%s,''parseCluster'',%s,''uuid'',''%s'',', ...
                            '''mccMode'',%s,''configFile'',''%s'')'], ...
                            tile_fullpaths_str, xyz_str, axisOrder, dataOrder, xyPixelSize, dz, skewAngle, string(reverse), ...
                            string(objectiveScan), string(IOScan), stitch_save_fname, tileInfoFullpath, stitchInfoDir, ...
                            stitchInfoFullpath, processedDirStr, string(DS), string(DSR), string(zarrFile), string(largeZarr), ...
                            strrep(mat2str(poolSize), ' ', ','), strrep(mat2str(batchSize), ' ', ','), strrep(mat2str(blockSize), ' ', ','), ...
                            strrep(mat2str(shardSize), ' ', ','), resampleType, strrep(num2str(resampleFactor, '%.10d,'), ' ', ''),  ...
                            strrep(mat2str(inputBbox), ' ', ','), strrep(mat2str(tileOutBbox), ' ', ','), tileOffset, blendMethod, overlapType, ...
                            string(xcorrShift), xyMaxOffset, zMaxOffset, strrep(mat2str(xcorrDownsample), ' ', ','), xcorrThresh, shiftMethod, ...
                            strrep(mat2str(axisWeight), ' ', ','), groupFile, string(isPrimaryCh), string(usePrimaryCoords), ...
                            num2str(padSize, '%d,'), strrep(num2str(outBbox, '%d,'), ' ', ''), ...
                            string(save16bit), tileIdx_str, flippedTile_str, processFunPath_str, strrep(mat2str(stitchMIP), ' ', ','), ...
                            string(bigStitchData), edgeArtifacts, strrep(mat2str(distBboxes), ' ', ','), string(saveMIP), ...
                            string(parseCluster), uuid, string(mccMode), configFile);

                        if exist(cur_tmp_fname, 'file') || (parseCluster && ~(masterCompute && xcorrShift && strcmpi(xcorrMode, 'primaryFirst') && isPrimaryCh))
                            % for cluster computing with master, check whether the job still alive. Otherwise, use waiting time
                            % for the check
                            if parseCluster
                                if ~(masterCompute && f == lastF && job_ids(f) == 0)
                                    if ~exist('imSize', 'var')
                                        if useProcessedData
                                            imSize = getImageSize(tile_processed_fullpaths{1});
                                        else
                                            imSize = getImageSize(tile_fullpaths{1});
                                        end
                                    end
                                    datasize = prod(imSize);
                                    totalSize = datasize * numel(tile_fullpaths);
                                    % assume for double
                                    totalDsize = totalSize * 4 / 1024^3;
                                    mem_factor = 2;

                                    lastFile = masterCompute && f==lastF;
    
                                    % allocate 5 time of the size
                                    memAllocate = totalDsize * mem_factor;
                                    job_id = job_ids(n, ncam, s, c, z);
                                    task_id = rem(f, 5000);
                                    
                                    [job_id, ~, submit_status] = generic_single_job_submit_wrapper(func_str, job_id, task_id, ...
                                        'jobLogFname', job_log_fname, 'jobErrorFname', job_log_error_fname, ...
                                        masterCompute=masterCompute, cpusPerTask=cpusPerTask, lastFile=lastFile, ...
                                        memAllocate=memAllocate, mccMode=mccMode, configFile=configFile);
    
                                    job_ids(n, ncam, s, c, z) = job_id;
                                    trial_counter(n, ncam, s, c, z) = trial_counter(n, ncam, s, c, z) + submit_status;
                                end
                            else
                                per_file_wait_time = unitWaitTime; % minite
                                temp_file_info = dir(cur_tmp_fname);
                                if minutes(datetime('now') - temp_file_info.date) < size(cur_tab, 1) * per_file_wait_time
                                    continue; 
                                else
                                    fclose(fopen(cur_tmp_fname, 'w'));
                                end
                            end
                        else
                            fclose(fopen(cur_tmp_fname, 'w'));
                        end

                        if ~parseCluster || (parseCluster && masterCompute && (f == lastF || ...
                                (xcorrShift && strcmpi(xcorrMode, 'primaryFirst') && isPrimaryCh ...
                                && job_ids(n, ncam, s, c, z) == -1)))
                            trial_counter(n, ncam, s, c, z) = trial_counter(n, ncam, s, c, z) + 1;
                            t0 = tic;
                            feval(str2func(['@()', func_str]));
                            fprintf('\nStitching result %s is generated. ', fsname)
                            toc(t0);
                            fprintf('\n');
                        end

                        % check if computing is done
                        if exist(stitch_save_fname, ftype)
                            is_done_flag(n, ncam, s, c, z) = true;
                            if exist(cur_tmp_fname, 'file')
                                delete(cur_tmp_fname);
                            end
                        end
                    end
                end
            end
        end
    end

    if streaming 
        stream_counter = stream_counter + 1;
    end

    % wait 30 seconds if some tasks are still computing
    if ~all(is_done_flag | trial_counter >= max_trial_num, 'all') || ...
            (parseCluster && any(job_status_flag & ~is_done_flag, 'all'))
        pause(5);
    end
end

if exist(stitching_tmp, 'dir')
    rmdir(stitching_tmp, 's');
end

end

