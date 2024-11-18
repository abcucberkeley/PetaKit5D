function [] = XR_matlab_stitching_wrapper(dataPaths, imageListFullpaths, varargin)
% Wrapper for the stitching pipeline.
% 
% 
% Required inputs :   
%          dataPaths : char or cell array. Directory paths for the tiles. Either a string for a single path or a cell array of paths for tiles within different directories.
% imageListFullpaths : char or cell array. Path(s) for image list. The order of imageListFullpaths need to map to the dataPaths. 
%
% Parameters (as 'specifier'-value pairs): 
%      resultDirName : char (default: 'matlab_stitch'). Result directory under data path.
%          streaming : true|false (default: false). True for real-time processing, and false for existing data
%    channelPatterns : a cell array (default: {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}).  Channel identifiers for included channels. 
%           multiLoc : true|{false}. Tiles are in multiple folders, with a group per folder.
%    processedDirStr : empty or char (default: ''). Processed directory name under dataPath.
% stitchInfoFullpath : empty or char (default: ''). Stitch info path. If provided, use the information to stitch the data.
%                 DS : true|false (default: false). Data is in deskewed space.
%                DSR : true|false (default: false). Data is in deskew/rotated space (with stage coordinates).
%        xyPixelSize : a number (default: 0.108). Pixel size in um.
%                 dz : a number (default: 0.5). Scan interval in um.
%          skewAngle : a number (default: 32.45). Skew angle (in degree) of the stage.
%            reverse : true|false (default: false). Inverse direction of z axis. 
%   parseSettingFile : true|false (default: false). Use the setting file to decide whether filp z stacks or not.
%          axisOrder : char (default: 'x,y,z'). Axis order mapping for coordinates in image list. With combinations of -, x, y, and z. '-yxz' means negative -y map to x, x maps to y, and z maps to z.
%          dataOrder : char (default: 'y,x,z'). Axis order mapping for data. 'y,x,z' means the first, second and third axes are y, x, and z, respectively.
%      objectiveScan : true|false (default: false). Objective scan. This is the scan as in the DS space.
%             IOScan : true|false (default: false). Inverted objective scan. This is the scan with the stage coordinates (DSR space). 
%           zarrFile : true|false (default: false). Use Zarr file as input.
%          largeFile : true|false (default: false). Use large-scale stitching strategy, i.e., use MIP slabs for registration and distance transforms.
%           poolSize : empty or 1x6 or 1x8 vectors (default: []). MIP slab generation downsampling factor for stitching large tiles. 
%          batchSize : 1x3 vector (default: [512, 512, 512]). Batch size per stitching task.
%          blockSize : 1x3 vector (default: [256, 256, 256]). Block/chunk size for zarr output.
%          shardSize : empty or 1x3 vector (default: []). If not empty, set shard size in zarr output.
%       resampleType : 'given'|'isotropic'|'xy_isotropic' (default: 'isotropic'). given: user-defined, xy_isotropic: xy isotropic, and z 
%     resampleFactor : empty or 1x1, 1x2 or 1x3 vector (default: []). Resampling factor. Empty: no resampling; axis order yxz.
%          inputBbox : empty or 1x6 vector (default: []). Input bounding box for crop for tiles. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%        tileOutBbox : empty or 1x6 vector (default: []). Crop tiles after preprocessing. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%         tileOffset : a number (default: 0). Offset added to tiles during preprocessing.
%        blendMethod : 'none'|'feather'|'mean'|'median'|'max' (default: 'feather'). Blending method for stitching. 'none' means taking half overlap from each tile; 'mean', 'median' and 'max' means use the mean, median or max intensity in the overlap region. 'feather' means feather blending.
%        overlapType : ''|'none'|'half'|'full' (default: ''). Method to handle overlap regions for blending (feather use 'full'). 'none', 'half', and 'full' means the overlap region to use for blending. If empty, use 'none' for none blending, and 'half' for others.
%         xcorrShift : true|false (default: false). Xcorr registration for stitching.
%        xyMaxOffset : a number (default: 300). Max offsets in voxel in xy axes in the registration.
%         zMaxOffset : a number (default: 50). Max offset in voxel in z axis in the registration.
%    xcorrDownsample : 1x3 vector (default: [2, 2, 1]). Downsampling factor for registration in yxz axis order.
%        xcorrThresh : a number (default: 0.25). Threshold for max xcorr value. The registration for a pair of tiles is not considered if belowing this threshold.
%            outBbox : empty or 1x6 vector (default: []). Only output the stitching result for the given output. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%          xcorrMode : 'primary'|'primaryFirst'|'all' (default: 'primaryFirst'). Xcorr registration mode. 'primary': use all time points in primary channel. 'primaryFirst': use the first time point in primary channel. 'all': registration within each time point and channel.
%        shiftMethod : 'local'|'global'|'grid'|'group' (default: 'grid'). Registration optimization method. 'local': use minimum spinning tree; 'global': use all pairs of neighbors; 'grid': use direct neighbors; 'group': use grid method within group, then across groups.
%         axisWeight : 1x3 vector (default: [1, 0.1, 10]). Weights for the axes to be considered during registration optimization in yxz orders.
%          groupFile : char (default: ''). File to define tile groups for 'group' option in shiftMethod. The format is the tile indices x, y, z, followed by gropu id.
%          primaryCh : empty or char (default: ''). The primary channel for registration. Must be one from channelPatterns.
%   usePrimaryCoords : true|false (default: false). Use the coordinates from the primary channel for stitching.
%          save16bit : true|false (default: true). Save 16bit result for stitch result; otherwise save as single. 
%      edgeArtifacts : a number (default: 2). The number of voxels from the border to erode to remove edge artifacts. 
%         distBboxes : empty or 1x6 or #Tilex6 vectors (default: []). Bounding boxes to define the distance transform for feather blending. Low weights are assigned outside the bounding boxes.
%            saveMIP : true|false (default: true). Save MIPs after stitching.
%          stitchMIP : empty or 1x1 or 1x3 bool vector (default: []). Stitching MIP for given axis. Order yxz. 
%       onlineStitch : true|{false}. Support for online stitch (with partial number of tiles).
%     processFunPath : empty or char (default: ''). Path for user-defined process function handle. Support .mat or .txt formats.
%       parseCluster : true|false (default: true). Use slurm cluster for the processing.
%      masterCompute : true|false (default: true). Master job node is involved in the processing.
%          jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%        cpusPerTask : a number (default: 8). The number of cpu cores per task for slurm job submission.
%               uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%        maxTrialNum : a number (default: 3). The max number of retries for a task.
%       unitWaitTime : a number (default: 0.1). The wait time per file in minutes to check whether the computing is done.
%            mccMode : true|false (default: false). Use mcc mode.
%         configFile : empty or char (default: ''). Path for the config file for job submission.
%
%
% Author: Xiongtao Ruan (02/18/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('imageListFullpaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'matlab_stitch', @ischar);
ip.addParameter('streaming', false, @islogical);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('multiLoc', false, @islogical); 
ip.addParameter('processedDirStr', '', @ischar);
ip.addParameter('stitchInfoFullpath', '', @ischar);
ip.addParameter('DS', false, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.5, @isnumeric);
ip.addParameter('skewAngle', 32.45, @isnumeric);
ip.addParameter('reverse', false, @islogical);
ip.addParameter('parseSettingFile', false, @islogical);
ip.addParameter('axisOrder', 'x,y,z', @ischar);
ip.addParameter('dataOrder', 'y,x,z', @ischar);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('IOScan', false, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('poolSize', [], @isnumeric);
ip.addParameter('batchSize', [512, 512, 512], @isnumeric);
ip.addParameter('blockSize', [256, 256, 256], @isnumeric);
ip.addParameter('shardSize', [], @isnumeric);
ip.addParameter('resampleType', 'xy_isotropic', @ischar);
ip.addParameter('resampleFactor', [], @isnumeric);
ip.addParameter('inputBbox', [], @isnumeric);
ip.addParameter('tileOutBbox', [], @isnumeric);
ip.addParameter('tileOffset', 0, @isnumeric);
ip.addParameter('blendMethod', 'feather', @ischar);
ip.addParameter('overlapType', '', @ischar);
ip.addParameter('xcorrShift', true, @islogical);
ip.addParameter('xyMaxOffset', 300, @(x) isvector(x) && numel(x) <= 2);
ip.addParameter('zMaxOffset', 50, @isnumeric);
ip.addParameter('xcorrDownsample', [2, 2, 1], @isnumeric);
ip.addParameter('xcorrThresh', 0.25, @isnumeric);
ip.addParameter('outBbox', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6));
ip.addParameter('xcorrMode', 'primaryFirst', @(x) ismember(x, {'primary', 'primaryFirst', 'all'}));
ip.addParameter('shiftMethod', 'grid', @ischar);
ip.addParameter('axisWeight', [1, 0.1, 10], @isnumeric);
ip.addParameter('groupFile', '', @ischar);
ip.addParameter('primaryCh', '', @(x) isempty(x) || ischar(x));
ip.addParameter('usePrimaryCoords', false, @islogical); 
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('edgeArtifacts', 0, @isnumeric);
ip.addParameter('distBboxes', [], @isnumeric);
ip.addParameter('saveMIP', true, @islogical);
ip.addParameter('stitchMIP', [], @(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3)));
ip.addParameter('onlineStitch', false, @(x) islogical(x));
ip.addParameter('processFunPath', '', @(x) isempty(x) || ischar(x) || iscell(x));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 8, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 0.1, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, imageListFullpaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
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
largeFile = pr.largeFile;
poolSize = pr.poolSize;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
shardSize = pr.shardSize;
resampleType = pr.resampleType;
resampleFactor = pr.resampleFactor;
inputBbox = pr.inputBbox;
tileOutBbox = pr.tileOutBbox;
tileOffset = pr.tileOffset;
blendMethod = pr.blendMethod;
overlapType = pr.overlapType;
xcorrShift = pr.xcorrShift;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
xcorrDownsample = pr.xcorrDownsample;
xcorrThresh = pr.xcorrThresh;
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
    stitching_rt = [dataPaths{1}, '/', resultDirName];
else
    if iscell(dataPaths)
        dataPaths = dataPaths{1};
    end
    if iscell(imageListFullpaths)
        imageListFullpaths = imageListFullpaths{1};
    end
    stitching_rt = [dataPaths, '/', resultDirName];
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
        [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPaths{1}, jobLogDir);
    else
        [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPaths, jobLogDir);
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
        stitch_parse_multi_loc_image_list_information(dataPaths, imageListFullpaths, streaming=streaming, ...
        onlineStitch=onlineStitch, stitchInfoFullpath=stitchInfoFullpath, stitchInfoPath=stitch_info_path, ...
        zarrFile=zarrFile, channelPatterns=channelPatterns, useProcessedData=useProcessedData, ...
        processedDirStr=processedDirStr, xcorrMode=xcorrMode, primaryCh=primaryCh);
else
    [tab, primary_tab, fullIter, Ch, Cam, stackn, nz, specifyCam, prefix, zlayerStitch, xcorrMode, stitchInfoFullpath] = ...
        stitch_parse_image_list_information(dataPaths, imageListFullpaths, streaming=streaming, ...
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
            tile_fullpath_1 =  [dataPaths{tab.did(1)}, '/', processedDirStr, '/', tile_fnames{1}];
        else
            tile_fullpath_1 =  [dataPaths, '/', processedDirStr, '/', tile_fnames{1}];
        end
    else
        if multiLoc
            tile_fullpath_1 =  [dataPaths{tab.did(1)}, '/', tile_fnames{1}];
        else
            tile_fullpath_1 =  [dataPaths, '/', tile_fnames{1}];
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
                                tile_fullpaths{f} = sprintf('%s/%s', dataPaths{cur_tab.did(f)}, tile_fnames{f});
                            end
                        else
                            tile_fullpaths = cellfun(@(x) [dataPaths, '/', x], tile_fnames, 'unif', 0);
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
                                    tile_processed_fullpaths{f} = sprintf('%s/%s/%s', dataPaths{cur_tab.did(f)}, processedDirStr, tile_fnames{f});
                                end
                            else
                                tile_processed_fullpaths = cellfun(@(x) [dataPaths, '/', processedDirStr, '/', x], tile_fnames, 'unif', 0);
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
                        if numel(tile_fullpaths) > 10 && parseCluster
                            tileInfoFullpath = sprintf('%s/stitchInfo/%s_tile_info.mat', stitching_rt, fsname);
                            if ~exist(tileInfoFullpath, 'file')
                                fprintf('Save tile paths and coordinates to disk...\n');
                                tileInfoTmppath = sprintf('%s/stitchInfo/%s_tile_info_%s.mat', stitching_rt, fsname, uuid);
                                save('-v7.3', tileInfoTmppath, 'tile_fullpaths', 'xyz', 'tileIdx', 'flippedTile');
                                movefile(tileInfoTmppath, tileInfoFullpath);
                            end
                            
                            tile_fullpaths_str = '{}';
                            xyz_str = '[]';
                            tileIdx_str = '[]';
                            flippedTile_str = '';
                        end

                        func_str = sprintf(['XR_stitching_frame_zarr_dev_v1(%s,%s,''axisOrder'',''%s'',''dataOrder'',''%s'',', ...
                            '''xyPixelSize'',%0.10f,''dz'',%0.10f,''skewAngle'',%0.10f,''reverse'',%s,''objectiveScan'',%s,', ...
                            '''IOScan'',%s,''resultPath'',''%s'',''tileInfoFullpath'',''%s'',''stitchInfoDir'',''%s'',', ...
                            '''stitchInfoFullpath'',''%s'',''processedDirStr'',''%s'',''DS'',%s,''DSR'',%s,''zarrFile'',%s,', ...
                            '''largeFile'',%s,''poolSize'',%s,''batchSize'',%s,''blockSize'',%s,''shardSize'',%s,''resampleType'',''%s'',', ...
                            '''resampleFactor'',[%s],''inputBbox'',%s,''tileOutBbox'',%s,''tileOffset'',%d,''blendMethod'',''%s'',', ...
                            '''overlapType'',''%s'',''xcorrShift'',%s,''xyMaxOffset'',%0.10f,''zMaxOffset'',%0.10f,''xcorrDownsample'',%s,', ...
                            '''xcorrThresh'',%0.10f,''shiftMethod'',''%s'',''axisWeight'',[%s],''groupFile'',''%s'',''isPrimaryCh'',%s,', ...
                            '''usePrimaryCoords'',%s,''outBbox'',[%s],''save16bit'',%s,', ...
                            '''tileIdx'',%s,''flippedTile'',[%s],''processFunPath'',%s,''stitchMIP'',%s,''bigStitchData'',%s,', ...
                            '''edgeArtifacts'',%0.10f,''distBboxes'',%s,''saveMIP'',%s,''parseCluster'',%s,''uuid'',''%s'',', ...
                            '''mccMode'',%s,''configFile'',''%s'')'], ...
                            tile_fullpaths_str, xyz_str, axisOrder, dataOrder, xyPixelSize, dz, skewAngle, string(reverse), ...
                            string(objectiveScan), string(IOScan), stitch_save_fname, tileInfoFullpath, stitchInfoDir, ...
                            stitchInfoFullpath, processedDirStr, string(DS), string(DSR), string(zarrFile), string(largeFile), ...
                            strrep(mat2str(poolSize), ' ', ','), strrep(mat2str(batchSize), ' ', ','), strrep(mat2str(blockSize), ' ', ','), ...
                            strrep(mat2str(shardSize), ' ', ','), resampleType, strrep(num2str(resampleFactor, '%.10d,'), ' ', ''),  ...
                            strrep(mat2str(inputBbox), ' ', ','), strrep(mat2str(tileOutBbox), ' ', ','), tileOffset, blendMethod, overlapType, ...
                            string(xcorrShift), xyMaxOffset, zMaxOffset, strrep(mat2str(xcorrDownsample), ' ', ','), xcorrThresh, shiftMethod, ...
                            strrep(mat2str(axisWeight), ' ', ','), groupFile, string(isPrimaryCh), string(usePrimaryCoords), ...
                            strrep(num2str(outBbox, '%d,'), ' ', ''), ...
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

