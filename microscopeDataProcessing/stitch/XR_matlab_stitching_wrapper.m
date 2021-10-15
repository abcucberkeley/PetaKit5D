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
%       'BlendMethod' : Blend method for overlap regions. Available: none, mean, max, median. Default: mean
%           'padSize' : Pad or crop the stitched image, empty (default) or 
%                       a 1X3 vector of integers (y, x, z). Postive for pad and negative for crop. 
%      'boundboxCrop' : Crop the stitched image by some bounding box, empty (default, no crop) 
%                       or a 3X2 vector for start and end indices of the bounding box (y, x, z).      
%        'zNormalize' : normalize background along z-axis by the median. Default: false.
%      'parseCluster' : Use slurm-based cluster computing. Default: true. 
%         'jobLogDir' : Log directory for the slurm jobs. 
%       'cpusPerTask' : Number of cpus for a job. Default: 12
%         'Save16bit' : true|{false}. Save final results as 16bit or single. 
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
% do not match the file pattern. 
% xruan (07/05/2021): add support for user defined resample (arbitary factor)
% xruan (07/15/2021): extend subiteration to unlimited number of sets
% xruan (08/25/2021): add support for channel-specific user functions
% xruan (10/13/2021): add support for cropping tiles; add support for
%   skewed space stitching with reference (for decon data)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @ischar);
ip.addRequired('imageListFileName', @ischar);
% ip.addParameter('Overwrite', true, @islogical);
ip.addParameter('Streaming', false, @islogical);
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('useExistDSR', false, @islogical); % use exist DSR for the processing
ip.addParameter('DSRDirstr', 'DSR', @ischar); % path for DSRDir str, if it is not true
ip.addParameter('useExistDecon', false, @islogical); % use exist decon for the processing (change to more generic decon, for not only DSR decon)
ip.addParameter('DeconDirstr', '', @ischar); % path for decon str, if it is not true
ip.addParameter('stitchInfoFullpath', '', @ischar); % use exist stitch info for stitching
ip.addParameter('DS', false, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('parseSettingFile', false, @islogical); % use setting file to decide whether filp Z stack or not.
ip.addParameter('axisOrder', 'x,y,z', @ischar);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('resampleType', 'xy_isotropic', @ischar); % by default use xy isotropic
ip.addParameter('resample', [], @isnumeric); % user-defined resample factor
ip.addParameter('CropToSize', [], @isnumeric); % size cropped to 
ip.addParameter('TileOffset', 0, @isnumeric); % offset added to tile
ip.addParameter('ffcorrect', false, @islogical);
ip.addParameter('Resolution', [0.108, 0.5], @isnumeric);
ip.addParameter('resultDir', 'matlab_stitch', @ischar);
ip.addParameter('BlendMethod', 'none', @ischar);
ip.addParameter('overlapType', '', @ischar); % '', 'none', 'half', or 'full'
ip.addParameter('xcorrShift', true, @islogical);
ip.addParameter('xyMaxOffset', 300, @isnumeric); % max offsets in xy axes
ip.addParameter('zMaxOffset', 50, @isnumeric); % max offsets in z axis
ip.addParameter('padSize', [], @(x) isnumeric(x) && (isempty(x) || numel(x) == 3));
ip.addParameter('boundboxCrop', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6));
ip.addParameter('zNormalize', false, @islogical);
ip.addParameter('onlyFirstTP', false, @islogical); % only compute first time point (for deciding cropping bouding box)
ip.addParameter('timepoints', [], @isnumeric); % stitch for given time points, nx1 
ip.addParameter('subtimepoints', [], @isnumeric); % stitch for given sub time points (subtimepoints), nx1
ip.addParameter('subsubtimepoints', [], @isnumeric); % stitch for given subsub time points (subtimepoints), nx1
ip.addParameter('xcorrMode', 'primaryFirst', @(x) strcmpi(x, 'primary') || strcmpi(x, 'primaryFirst') || strcmpi(x, 'all')); % 'primary': choose one channel as primary channel, 
                                                                                          % 'all': xcorr shift for each channel, 
                                                                                          % 'primaryFirst': the primary channel of first time point
ip.addParameter('primaryCh', '', @(x) isempty(x) || ischar(x)); % format: CamA_ch0. If it is empty, use the first channel as primary channel
ip.addParameter('usePrimaryCoords', false, @islogical); 
ip.addParameter('Save16bit', false, @islogical);
ip.addParameter('pipeline', 'matlab', @(x) strcmpi(x, 'matlab') || strcmpi(x, 'zarr'));
ip.addParameter('processFunPath', '', @(x) isempty(x) || ischar(x) || iscell(x)); % path of user-defined process function handle
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 8, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 0.1, @isnumeric);
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2020b; matlab -nodisplay -nosplash -nodesktop -r', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);

ip.parse(dataPath, imageListFileName, varargin{:});

pr = ip.Results;
% Overwrite = pr.Overwrite;
Streaming = pr.Streaming;
ChannelPatterns = pr.ChannelPatterns;
useExistDSR = pr.useExistDSR;
DSRDirstr = pr.DSRDirstr;
useExistDecon = pr.useExistDecon;
DeconDirstr = pr.DeconDirstr;
stitchInfoFullpath = pr.stitchInfoFullpath;
DS = pr.DS;
DSR = pr.DSR;
Reverse = pr.Reverse;
parseSettingFile = pr.parseSettingFile;
axisOrder = pr.axisOrder;
ObjectiveScan = pr.ObjectiveScan;
resampleType = pr.resampleType;
resample = pr.resample;
CropToSize = pr.CropToSize;
TileOffset = pr.TileOffset;
ffcorrect = pr.ffcorrect;
Resolution = pr.Resolution;
resultDir = pr.resultDir;
BlendMethod = pr.BlendMethod;
overlapType = pr.overlapType;
xcorrShift = pr.xcorrShift;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
padSize = pr.padSize;
boundboxCrop = pr.boundboxCrop;
zNormalize = pr.zNormalize;
onlyFirstTP = pr.onlyFirstTP;
timepoints = pr.timepoints;
subtimepoints = pr.subtimepoints;
subsubtimepoints = pr.subtimepoints;
xcorrMode = pr.xcorrMode;
primaryCh = pr.primaryCh;
usePrimaryCoords = pr.usePrimaryCoords;
Save16bit = pr.Save16bit;
pipeline = pr.pipeline;
processFunPath = pr.processFunPath;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes = pr.cpuOnlyNodes;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
SlurmParam = pr.SlurmParam;
MatlabLaunchStr = pr.MatlabLaunchStr;

switch pipeline
    case 'matlab'
        stitch_function_str = 'XR_stitching_frame_v1';
    case 'zarr'
        stitch_function_str = 'XR_stitching_frame_zarr_dev_v1';
end

px = Resolution(1);
dz = Resolution(end);

% make root directory
stitching_rt = [dataPath, '/', resultDir];
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

if useExistDecon
    useExistDSR = false;
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str, jobLogDir] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
end

% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end

%% parse image list information
% read image list csv file
t = readtable(imageListFileName, 'Delimiter','comma');
% t_column_name = t.Properties.VariableNames;
fn = t.Filename;
% remove records for partial volume files
partialvol = ~cellfun(@isempty, regexpi(fn, '_part\d+.tif'));
if any(partialvol)
    t(partialvol, :) = [];
    fn = t.Filename;
end   

specifyCam = true;
if all(~cellfun(@isempty, regexp(fn, '_Cam\w_ch', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>\d+)x_(?<y>\d+)y_(?<z>\d+)z_(?<t>\d+)t.tif';
elseif all(~cellfun(@isempty, regexp(fn, '_ch[0-9]_', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>\d+)x_(?<y>\d+)y_(?<z>\d+)z_(?<t>\d+)t.tif';
    specifyCam = false;
end

tmp = regexpi(fn, expression, 'names');

matched_inds = true(numel(tmp), 1);

for f = 1:numel(tmp)
    if isempty(tmp{f})
        matched_inds(f) = false;
        continue;
    end
    
    t.prefix{f} = tmp{f}.prefix;
    t.Iter(f) = str2double(tmp{f}.Iter);
    t.subIter{f} = tmp{f}.subIter;
    t.fullIter{f} = [tmp{f}.Iter, tmp{f}.subIter];
    if specifyCam
        t.camera(f) = (tmp{f}.Cam);
    else
        % use N to represent cam if it is not contained in the filename
        t.camera(f) = 'N';
    end
    t.ch(f) = str2double(tmp{f}.ch);
    t.stack(f) = str2double(tmp{f}.stack);
    t.laser(f) = str2double(tmp{f}.laser);
    t.abstime(f) = str2double(tmp{f}.abstime);
    t.fpgatime(f) = str2double(tmp{f}.fpgatime);
    t.x(f) = str2double(tmp{f}.x);
    t.y(f) = str2double(tmp{f}.y);
    t.z(f) = str2double(tmp{f}.z);
    t.t(f) = str2double(tmp{f}.t);
end

t = t(matched_inds, :);

prefix = unique(t.prefix);
if ~isempty(prefix)
    prefix = prefix{1};
else
    prefix = '';
end
fullIter = unique(t.fullIter);
Ch = unique(t.ch);
Cam = unique(t.camera);
stackn = unique(t.stack);

% check whether the image files in the image list file exist 
dir_info = dir([dataPath, '/', '*.tif']);
imageFnames = {dir_info.name}';
if isempty(dir_info)
    if useExistDSR
        dir_info = dir([dataPath, DSRDirstr, '/', '*.tif']);
        imageFnames = {dir_info.name}';
    elseif useExistDecon
        dir_info = dir([dataPath, '/', DeconDirstr, '/', '*.tif']);
        imageFnames = cellfun(@(x) [x(1 : end - 10), '.tif'], {dir_info.name}', 'unif', 0);
    else
        error('The tiles do not exist!');
    end
end

% filter filenames by channel patterns
include_flag = false(numel(imageFnames), 1);
imageFullnames = cellfun(@(x) [dataPath, '/', x], imageFnames, 'unif', 0);
for c = 1 : numel(ChannelPatterns)
    include_flag = include_flag | contains(imageFullnames, ChannelPatterns{c}) | contains(imageFullnames, regexpPattern(ChannelPatterns{c}));
end
imageFnames = imageFnames(include_flag);

[~, tFsnames] = fileparts(t.Filename);
image_file_exist_flag = true(numel(tFsnames), 1);
for f = 1 : numel(t.Filename)
    if ~any(contains(imageFnames, tFsnames{f})) % || (useExistDecon && ~any(contains(imageFnames, tFsnames{f})))
        image_file_exist_flag(f) = false;
    end
end
if ~all(image_file_exist_flag)
    warning('Some files in the image list file do not exist! Ignore them in the stitching')
    disp(t.Filename(~image_file_exist_flag));
    if ~Streaming
        t(~image_file_exist_flag, :) = [];
        
        fullIter = unique(t.fullIter);
        Ch = unique(t.ch);
        Cam = unique(t.camera);
        % ntiles = numel(unique(t.x)) * numel(unique(t.y)) * numel(unique(t.z));
        stackn = unique(t.stack);
    end
end

if ~isempty(stitchInfoFullpath) 
    if exist(stitchInfoFullpath, 'file')
        % xcorrShift = true;
        xcorrMode = 'stitchInfo';
    else
        error('The user defined stitch info file %s does not exist!', stitchInfoFullpath);
    end
end

if (strcmpi(xcorrMode, 'primary') || strcmpi(xcorrMode, 'primaryFirst'))
    % first check whether the format is right
    if ~isempty(primaryCh)
        fprintf('The primary channel is %s.', primaryCh);
        if specifyCam && isempty(regexpi(primaryCh, 'Cam[a-z]_ch[0-9]', 'once'))
            error("primaryCh must be empty or with format 'Cam[a-z]_ch[0-9]'");
        elseif ~specifyCam && isempty(regexpi(primaryCh, 'ch[0-9]', 'once'))
            error("primaryCh must be empty or with format 'ch[0-9]'");            
        end
        pCam = primaryCh(4);
        pCh = str2double(primaryCh(end));
        if ~any(contains(fn, ['_', primaryCh, '_'], 'IgnoreCase', true))
            error('The given primary channel %s does not exist!', primaryCh);
        end
    else
        warning('The primary channel is not set, use the first available channel as primary channel...');
        pCam = [];
        for ncam = 1 : numel(Cam)
            for c = 1 : numel(Ch)
                if specifyCam
                    primaryCh = sprintf('Cam%s_ch%d', Cam(ncam), Ch(c));
                else
                    primaryCh = sprintf('ch%d', Ch(c));
                end
                if any(contains(fn, ['_', primaryCh, '_'], 'IgnoreCase', true))
                    pCam = Cam(ncam);
                    pCh =  Ch(c);
                    break;
                end
            end
            if ~isempty(pCam)
                break;
            end
        end       
        fprintf('Set primary channel as %s.', primaryCh);
    end

    % reorder Ch and Cam to make primary channel the first
    if Ch(1) ~= pCh
        Ch = [pCh; Ch(Ch ~= pCh)];
    end
    if ~strcmpi(Cam(1), pCam)
        ind = regexpi(Cam', pCam);
        Cam = Cam([ind, 1 : ind - 1, ind + 1 : end]);
    end
end

if onlyFirstTP
    fullIter = fullIter(1);
    timepoints = [];
end

% xruan (01/21/2021): modify to separate timepoints, subtimepoints, subsubtimepoints
if ~isempty(timepoints)
    Iter = cellfun(@(x) str2double(x(1 : 4)), fullIter);

    fullIter = fullIter(ismember(Iter, timepoints(:)));
    if ~isempty(subtimepoints)
        subIter = cellfun(@(x) str2double(x(6 : 9)), fullIter);
        fullIter = fullIter(ismember([Iter, subIter], [timepoints(:), subtimepoints(:)], 'rows'));
    end
    
    if ~isempty(subsubtimepoints)
        subIter = cellfun(@(x) str2double(x(6 : 9)), fullIter);
        subSubIter = cellfun(@(x) str2double(x(11 : 14)), fullIter);
        fullIter = fullIter(ismember([Iter, subIter, subSubIter], [timepoints(:), subtimepoints(:), subsubtimepoints(:)], 'rows'));
    end
end

if ischar(processFunPath)
    processFunPath = repmat({processFunPath}, numel(ChannelPatterns));
end


%% do stitching computing
row_exist_flag = true(numel(fullIter), numel(Cam), numel(stackn), numel(Ch)); % flag for whether the run exists.
is_done_flag = false(numel(fullIter), numel(Cam), numel(stackn), numel(Ch));
trial_counter = zeros(numel(fullIter), numel(Cam), numel(stackn), numel(Ch));
max_trial_num = maxTrialNum;

if parseCluster
    job_ids = -ones(numel(fullIter), numel(Cam), numel(stackn), numel(Ch));
    job_status_flag = false(numel(fullIter), numel(Cam), numel(stackn), numel(Ch));
end

% predefine stitchInfo when xcorrMode is 'primaryFirst'
if strcmp(xcorrMode, 'primaryFirst')
    primary_t = t(t.ch == Ch(1) & t.camera == Cam(1) & strcmp(t.fullIter, fullIter{1}) & t.stack == stackn(1), :);
    if isempty(primary_t) 
        error('The Image List Info for the primary channel for the first time point does not exist!');
    end
    
    p_laser = unique(primary_t.laser);
    p_abstime = unique(primary_t.abstime);
    p_fpgatime = primary_t.fpgatime(1);

    if numel(p_laser) > 1
        p_laser = p_laser(1);
    end
    if specifyCam
        stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs.mat', ...
            stitch_info_path, prefix, fullIter{1}, Cam(1), Ch(1), stackn(1), p_laser, p_abstime, p_fpgatime);
    else
        stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs.mat', ...
            stitch_info_path, prefix, fullIter{1}, Ch(1), stackn(1), p_laser, p_abstime, p_fpgatime);        
    end
end

% set wait counter for streaming option
if Streaming
    stream_counter = 0;
    stream_max_counter = 100 * size(t, 1);
end

while ~all(is_done_flag | trial_counter >= max_trial_num, 'all')
    lastF = find(~is_done_flag & trial_counter < maxTrialNum, 1, 'last');
    for n = 1:numel(fullIter)
        for ncam = 1:numel(Cam)
            for s = 1:numel(stackn)
                for c = 1:numel(Ch)                    
                    if is_done_flag(n, ncam, s, c) || trial_counter(n, ncam, s, c) >= max_trial_num
                        continue;
                    end
                    if ~row_exist_flag(n, ncam, s, c)
                        is_done_flag(n, ncam, s, c) = true;
                        continue;
                    end
                    
                    f = sub2ind([numel(fullIter), numel(Cam), numel(stackn), numel(Ch)], n, ncam, s, c);
                    task_id = f;
                    cur_t = t(t.ch == Ch(c) & t.camera == Cam(ncam) & strcmp(t.fullIter, fullIter{n}) & t.stack == stackn(s), :);

                    % obtain filenames                    
                    if isempty(cur_t)
                        row_exist_flag(n, ncam, s, c) = false;
                        is_done_flag(n, ncam, s, c) = true;

                        % if the primary channel is missing, then skip
                        % other channels if it is primary mode
                        if (ncam == 1 && c == 1) && strcmp(xcorrMode, 'primary')
                            row_exist_flag(1, 1 : end, 1, 1 : end) = false;
                            is_done_flag(1, 1 : end, 1, 1 : end) = true;
                        end
                        continue;
                    end
                    
                    % for primaryfirst option, check whether the number of
                    % tiles are the same as primary channel
                    if strcmpi(xcorrMode, 'primaryFirst') && size(cur_t, 1) ~= size(primary_t, 1)
                        fprintf('The number of tiles is different from the number of tiles in the primary channel, skip it!\n');
                        disp(cur_t.Filename);
                        is_done_flag(n, ncam, s, c) = true;
                        continue;
                    end
                    
                    laser = unique(cur_t.laser);
                    abstime = unique(cur_t.abstime);
                    % reltime = unique(cur_t.t);
                    fpgatime = cur_t.fpgatime(1);
                    xyz = [cur_t.StageX_um_, cur_t.StageY_um_, cur_t.StageZ_um_];
                    tileNum = [numel(unique(cur_t.x)), numel(unique(cur_t.y)), numel(unique(cur_t.z))];
                    
                    if numel(laser) > 1
                        laser = laser(1);
                    end
                    
                    tile_fnames = cur_t.Filename;
                    tile_fullpaths = cellfun(@(x) [dataPath, '/', x], tile_fnames, 'unif', 0);
                    
                    % parse setting file 
                    flippedTile = [];
                    if parseSettingFile
                        settingInfo = XR_parseSettingFiles_wrapper(tile_fullpaths);
                        flippedTile = [settingInfo.StageInterval] < 0;
                    end
                    
                    % check if files exist for streaming option, and also
                    % if useExistDSR is true, check the DSR files exist.
                    if Streaming
                        is_tile_exist = cellfun(@(x) exist(x, 'file'), tile_fullpaths);
                        if ~all(is_tile_exist)
                            stream_counter = stream_counter + 1;
                            continue;
                        else 
                            stream_counter = 0;
                        end
                    end
                    if useExistDSR
                        % DSRDirstr = 'DSR';
                        tile_dsr_fullpaths = cellfun(@(x) [dataPath, '/', DSRDirstr, '/', x], tile_fnames, 'unif', 0);
                        is_tile_dsr_exist = cellfun(@(x) exist(x, 'file'), tile_dsr_fullpaths);
                        if Streaming
                            if ~all(is_tile_dsr_exist)
                                stream_counter = stream_counter + 1;
                                continue;
                            else
                                stream_counter = 0;
                            end
                        else
                            if ~all(is_tile_dsr_exist) 
                                % is_done_flag(n, ncam, s, c) = true;
                                trial_counter(n, ncam, s, c) = trial_counter(n, ncam, s, c) + 1;
                                continue; 
                            end
                        end
                    else
                        DSRDirstr = '';
                    end
                    
                    if useExistDecon
                        tile_decon_fullpaths = cellfun(@(x) [dataPath, '/', DeconDirstr, '/', x(1 : end - 4), '_decon.tif'], tile_fnames, 'unif', 0);
                        is_tile_decon_exist = cellfun(@(x) exist(x, 'file'), tile_decon_fullpaths);
                        if Streaming
                            if ~all(is_tile_decon_exist)
                                stream_counter = stream_counter + 1;
                                continue;
                            else
                                stream_counter = 0;
                            end
                        else
                            if ~all(is_tile_decon_exist) 
                                is_done_flag(n, ncam, s, c) = true;
                                continue; 
                            end
                        end
                    else
                        DeconDirstr = '';
                    end

                    % for secondary channels, first check whether the
                    % stitching info file exist
                    isPrimaryCh = true;
                    % the stitch Info full path for options except 'primaryFirst'
                    if strcmpi(xcorrMode, 'primaryFirst')
                        if ~(n == 1 && ncam == 1 && s == 1 && c == 1)
                            isPrimaryCh = false;
                        end                        
                    elseif strcmpi(xcorrMode, 'stitchInfo')
                        isPrimaryCh = false;
                    else
                        if specifyCam
                            stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs.mat', ...
                                stitch_info_path, prefix, fullIter{n}, Cam(ncam), Ch(c), stackn(s), laser, abstime, fpgatime);
                        else
                            stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs.mat', ...
                                stitch_info_path, prefix, fullIter{n}, Ch(c), stackn(s), laser, abstime, fpgatime);                            
                        end
                    end
                    
                    if true || xcorrShift 
                        % get the primary info path for secondary channels
                        % for 'primary' option
                        if strcmpi(xcorrMode, 'primary') && ~(ncam == 1 && c == 1)
                            isPrimaryCh = false;
                            primary_t = t(t.ch == Ch(1) & t.camera == Cam(1) & strcmp(t.fullIter, fullIter{n}) & t.stack == stackn(s), :);
                            p_laser = unique(primary_t.laser);
                            p_abstime = unique(primary_t.abstime);
                            p_fpgatime = primary_t.fpgatime(1);

                            if numel(p_laser) > 1
                                p_laser = p_laser(1);
                            end
                            
                            if specifyCam
                                stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs.mat', ...
                                    stitch_info_path, prefix, fullIter{n}, Cam(1), Ch(1), stackn(s), p_laser, p_abstime, p_fpgatime);
                            else
                                stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs.mat', ...
                                    stitch_info_path, prefix, fullIter{n}, Ch(1), stackn(1), p_laser, p_abstime, p_fpgatime);        
                            end
                        end
                        
                        % for secondary channels, if the stitchInfo file
                        % not exist, wait the stitching for the primary
                        % channel
                        if ~isPrimaryCh && ~exist(stitchInfoFullpath, 'file')
                            continue;
                        end
                    end
                    
                    % also use flag based check of completion, to support
                    % distributed computing with same submission
                    if specifyCam
                        stitch_save_fsname = sprintf('%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs', ...
                            stitching_rt, prefix, fullIter{n}, Cam(ncam), Ch(c), stackn(s), laser, abstime, fpgatime);
                        cur_tmp_fname = sprintf('%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs.tmp', ...
                            stitching_rt, prefix, fullIter{n}, Cam(ncam), Ch(c), stackn(s), laser, abstime, fpgatime);
                    else
                        stitch_save_fsname = sprintf('%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs', ...
                            stitching_rt, prefix, fullIter{n}, Ch(c), stackn(s), laser, abstime, fpgatime);
                        cur_tmp_fname = sprintf('%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs.tmp', ...
                            stitching_rt, prefix, fullIter{n}, Ch(c), stackn(s), laser, abstime, fpgatime);
                    end
                    
                    switch pipeline
                        case 'matlab'
                            stitch_save_fname = [stitch_save_fsname, '.tif'];
                            ftype = 'file';
                        case 'zarr'
                            stitch_save_fname = [stitch_save_fsname, '.zarr'];
                            ftype = 'dir';
                    end
                                
                    % first check if the computing is done
                    if exist(stitch_save_fname, ftype)
                        is_done_flag(n, ncam, s, c) = true; 
                        if exist(cur_tmp_fname, 'file')
                            delete(cur_tmp_fname);
                        end
                        continue;
                    end
                    
                    tile_fullpaths_str = sprintf('{''%s''}', strjoin(tile_fullpaths, ''','''));
                    xyz_str = strrep(mat2str(xyz), ' ', ',');  
                    cind = cellfun(@(x) contains(tile_fullpaths{1}, x), ChannelPatterns);
                    func_str = sprintf(['%s(%s,%s,''axisOrder'',''%s'',''px'',%0.10f,''dz'',%0.10f,', ...
                        '''Reverse'',%s,''ObjectiveScan'',%s,''resultDir'',''%s'',''stitchInfoDir'',''%s'',''stitchInfoFullpath'',''%s'',', ...
                        '''DSRDirstr'',''%s'',''DeconDirstr'',''%s'',''DS'',%s,''DSR'',%s,''resampleType'',''%s'',''resample'',[%s],', ...
                        '''CropToSize'',%s,''TileOffset'',%d,''BlendMethod'',''%s'',''overlapType'',''%s'',''xcorrShift'',%s,''xyMaxOffset'',%0.10f,', ...
                        '''zMaxOffset'',%0.10f,''isPrimaryCh'',%s,''usePrimaryCoords'',%s,''padSize'',[%s],''boundboxCrop'',[%s],', ...
                        '''zNormalize'',%s,''Save16bit'',%s,''tileNum'',[%s],''flippedTile'',[%s],''processFunPath'',''%s'')'], ...
                        stitch_function_str, tile_fullpaths_str, xyz_str, axisOrder, px, dz, string(Reverse), string(ObjectiveScan), resultDir, stitchInfoDir, ...
                        stitchInfoFullpath, DSRDirstr, DeconDirstr, string(DS), string(DSR), resampleType, strrep(num2str(resample, '%.10d,'), ' ', ''), ...
                        strrep(mat2str(CropToSize), ' ', ','), TileOffset, BlendMethod, overlapType, string(xcorrShift), xyMaxOffset, zMaxOffset, string(isPrimaryCh), ...
                        string(usePrimaryCoords), num2str(padSize, '%d,'),  strrep(num2str(boundboxCrop, '%d,'), ' ', ''), string(zNormalize),  string(Save16bit), ...
                        strrep(num2str(tileNum, '%d,'), ' ', ''), strrep(num2str(flippedTile, '%d,'), ' ', ''), processFunPath{cind});

                    if exist(cur_tmp_fname, 'file') || (parseCluster && ~(masterCompute && strcmpi(xcorrMode, 'primaryFirst') && isPrimaryCh))
                        % for cluster computing with master, check whether
                        % the job still alive. Otherwise, use waiting time
                        % for the check
                        if parseCluster
                            job_status = check_slurm_job_status(job_ids(n, ncam, s, c), rem(task_id, 1000));

                            % kill the last pending job and use master node do the computing.
                            if job_status == 0.5 && (masterCompute && f == lastF)
                                system(sprintf('scancel %d_%d', job_ids(n, ncam, s, c), rem(task_id, 1000)), '-echo');
                            end

                            % if the job is still running, skip it. 
                            if job_status == 1 
                                continue;
                            end

                            % If there is no job, submit a job
                            if job_status == -1 && ~(masterCompute && f == lastF)
                                % check if memory is enough
                                if useExistDSR
                                    dir_info = dir(tile_dsr_fullpaths{1});
                                elseif useExistDecon
                                    dir_info = dir(tile_decon_fullpaths{1});                                    
                                else
                                    dir_info = dir(tile_fullpaths{1});                                    
                                end
                                datasize = dir_info.bytes;
                                totalSize = datasize * numel(tile_fullpaths);
                                % assume for double
                                totalDsize = totalSize * 8 / 1024^3;
                                if strcmp(BlendMethod, 'mean') || strcmp(BlendMethod, 'median')
                                    mem_factor = 10;
                                else
                                    mem_factor = 8;
                                end
                                if strcmp(pipeline, 'zarr')
                                    mem_factor = 1;
                                end
                                    
                                % allocate 5 time of the size
                                if cpusPerTask * 20 < totalDsize * mem_factor
                                    cpusPerTask = min(24, ceil(totalDsize * mem_factor / 20))
                                end
                                
                                matlab_setup_str = 'setup([],true)';
                                matlab_cmd = sprintf('%s;tic;%s;toc', matlab_setup_str, func_str);
                                stitch_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                                cmd = sprintf('sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s --wrap="echo Matlab command:  \\\"%s\\\"; %s"', ...
                                    rem(task_id, 5000), job_log_fname, job_log_error_fname, cpusPerTask, SlurmParam, ...
                                    slurm_constraint_str, matlab_cmd, stitch_cmd);
                                [status, cmdout] = system(cmd, '-echo');

                                job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                                job_id = str2double(job_id{1}{1});
                                job_ids(n, ncam, s, c) = job_id;
                            end
                        else
                            per_file_wait_time = unitWaitTime; % minite
                            temp_file_info = dir(cur_tmp_fname);
                            if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 < numel(cur_t) * per_file_wait_time
                                continue; 
                            else
                                fclose(fopen(cur_tmp_fname, 'w'));
                            end
                        end
                    else
                        fclose(fopen(cur_tmp_fname, 'w'));
                    end
                                        
                    if ~parseCluster || (parseCluster && masterCompute && (f == lastF || ...
                            (xcorrShift && strcmpi(xcorrMode, 'primaryFirst') && isPrimaryCh && job_ids(n, ncam, s, c) == -1)))
                        trial_counter(n, ncam, s, c) = trial_counter(n, ncam, s, c) + 1;
                        tic
                        feval(str2func(['@()', func_str]));
                        toc
                    end
                    
                    % check if computing is done
                    if exist(stitch_save_fname, ftype)
                        is_done_flag(n, ncam, s, c) = true;
                        if exist(cur_tmp_fname, 'file')
                            delete(cur_tmp_fname);
                        end
                    end
                end
            end
        end
    end
    
    % wait 30 seconds if some tasks are still computing
    if ~all(is_done_flag | trial_counter >= max_trial_num, 'all') || ...
            (parseCluster && any(job_status_flag & ~is_done_flag, 'all'))
        pause(30);
    end
        
    % exit the job if no new images are transferred.
    if Streaming && stream_counter > stream_max_counter
        break;
    end
end

rmdir(stitching_tmp, 's');

end

