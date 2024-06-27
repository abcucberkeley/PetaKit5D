function [] = XR_microscopeAutomaticProcessing(dataPaths, varargin)
% Automatic image processing pipeline for microscopy data. It can perform
% deskew w/o rotation, and stitching. It supports both existing data in the 
% directory or real-time transferring of data from other sources, i.e., microscope. 
%
%
% Required inputs :   
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings. 
%
% Parameters (as 'specifier'-value pairs): 
%           overwrite : true|false or 1x3 bool vector (default: false). Overwrite existing results
%           streaming : true|false (default: true). True for real-time processing, and false for existing data
%     channelPatterns : a cell array (default: {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}).  Channel identifiers for included channels. 
%           skewAngle : a number (default: 32.45). Skew angle (in degree) of the stage.
%                  dz : a number (default: 0.5). Scan interval in um.
%         xyPixelSize : a number (default: 0.108). Pixel size in um.
%             reverse : true|false (default: false). Inverse direction of z axis. 
%       objectiveScan : true|false (default: false). Objective scan.
%          zStageScan : true|false (default: false). Z stage scan (orthogonal to objective scan).
%           save16bit : 1x1 or 1x2 bool vector (default: [true, true]). Save 16bit result for deskew/rotate and stitch. 
%       dzFromEncoder : true|false (default: false). Estimate dz from encoder positions.
%            zarrFile : true|false (default: false). Use Zarr file as input.
%            saveZarr : true|false (default: false). Save results as Zarr files.
%         save3DStack : true|false (default: true). Save 3D stack for DS and DSR. 
%              deskew : true|false (default: true). Deskew the data.
%              rotate : true|false (default: true). Rotate deskewed data.
%              stitch : true|false (default: false). Stitch deskewed and rotated data.
%    parseSettingFile : true|false (default: false). Use the setting file to decide whether filp z stacks or not.
%          flipZstack : true|false (default: false). Flip z stacks.
%         DSRCombined : true|false (default: true). Use combined processing for deskew and rotation.
%        FFCorrection : true|false (default: false). Flat-field correction.
%           BKRemoval : true|false (default: false). REmove background during flat-field correction.
%          lowerLimit : a number between 0 and 1 (default: 0.4). Lower limit to cap the flat field image. 
%         constOffset : empty or a number (default: []). If empty, add the background; if not, add the const offset after flat field correction.
%        FFImagePaths : empty or a cell array of paths for corresponding channels (default: {'', '', ''}). Flat field image paths.
%     backgroundPaths : empty or a cell array of paths for corresponding channels (default: {'', '', ''}). Background image paths.
%        resampleType : 'given'|'isotropic'|'xy_isotropic' (default: 'isotropic'). given: user-defined, xy_isotropic: xy isotropic, and z 
%      resampleFactor : empty or 1x1, 1x2 or 1x3 vector (default: []). Resampling factor. Empty: no resampling; axis order yxz.
%           inputBbox : empty or 1x6 vector (default: []). Input bounding box for crop. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
% stitchResultDirName : empty or char (default: ''). Result directory name for stitching.
%  imageListFullpaths : empty or char or cell array (default: ''). Path(s) for image list.
%           axisOrder : char (default: 'xyz'). Axis order mapping for coordinates in image list. With combinations of -, x, y, and z. '-yxz' means negative -y map to x, x maps to y, and z maps to z.
%         blendMethod : 'none'|'feather'|'mean'|'median'|'max' (default: 'none'). Blending method for stitching.
%          xcorrShift : true|false (default: false). Xcorr registration for stitching.
%           xcorrMode : 'primary'|'primaryFirst'|'all' (default: 'primaryFirst'). Xcorr registration mode. 'primary': use all time points in primary channel. 'primaryFirst': use the first time point in primary channel. 'all': registration within each time point and channel.
%         xyMaxOffset : a number (default: 300). Max offsets in voxel in xy axes in the registration.
%          zMaxOffset : a number (default: 50). Max offset in voxel in z axis in the registration.
%       edgeArtifacts : a number (default: 2). The number of voxels from the border to erode to remove edge artifacts. 
%           primaryCh : empty or char (default: ''). The primary channel for registration. Must be one from channelPatterns.
%           stitchMIP : empty or 1x1 or 1x3 bool vector (default: []). Stitching MIP for given axis. Order yxz. 
%        onlineStitch : true|false (default: false). Support for online stitch (with partial number of tiles).
%   generateImageList : ''|'from_encoder'|'from_sqlite' (default: ''). Method to generate image list.
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%           jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%         cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%         maxTrialNum : a number (default: 3). The max number of retries for a task.
%        unitWaitTime : a number (default: 1). The wait time per file in minutes to check whether the computing is done.
%       minModifyTime : a number (default: 10). The minimum time in minutes for the latest modified file to decide whether it is fully transferred.
%       maxModifyTime : a number (default: 10). The maximum time in minutes to check whether there are coming new files.
%      maxWaitLoopNum : a number (default: 10). Number of maximum loops without any computing.
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%
%
% Author: Xiongtao Ruan (03/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 3) && islogical(x));
ip.addParameter('streaming', true,  @islogical);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('skewAngle', 32.45, @isscalar);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('reverse', true, @islogical);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('zStageScan', false, @islogical);
ip.addParameter('save16bit', [true, true], @(x) (numel(x) == 1 || numel(x) == 2) && islogical(x));
ip.addParameter('dzFromEncoder', false, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('saveZarr', false, @islogical);
ip.addParameter('save3DStack', true , @islogical);
% pipeline steps
ip.addParameter('deskew', true, @islogical);
ip.addParameter('rotate', true, @islogical);
ip.addParameter('stitch', false, @islogical);
% deskew and rotation options
ip.addParameter('parseSettingFile', false, @islogical);
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('DSRCombined', true, @islogical); 
ip.addParameter('FFCorrection', false, @islogical);
ip.addParameter('BKRemoval', false, @islogical);
ip.addParameter('lowerLimit', 0.4, @isnumeric);
ip.addParameter('constOffset', [], @(x) isnumeric(x));
ip.addParameter('FFImagePaths', {'','',''}, @iscell);
ip.addParameter('backgroundPaths', {'','',''}, @iscell);
ip.addParameter('resampleType', 'isotropic', @ischar);
ip.addParameter('resampleFactor', [], @isnumeric);
ip.addParameter('inputBbox', [], @isnumeric);
% stitch parameters
ip.addParameter('stitchResultDirName', '', @ischar);
ip.addParameter('imageListFullpaths', '', @(x) ischar(x) || iscell(x));
ip.addParameter('axisOrder', 'xyz', @(x) ischar(x));
ip.addParameter('blendMethod', 'none', @ischar);
ip.addParameter('xcorrShift', false, @islogical);
ip.addParameter('xcorrMode', 'primaryFirst', @(x) ismember(lower(x), {'primary', 'primaryfirst', 'all'}));
ip.addParameter('xyMaxOffset', 300, @isnumeric);
ip.addParameter('zMaxOffset', 50, @isnumeric);
ip.addParameter('edgeArtifacts', 2, @isnumeric);
ip.addParameter('primaryCh', '', @ischar);
ip.addParameter('stitchMIP', [], @(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3)));
ip.addParameter('onlineStitch', false, @(x) islogical(x));
ip.addParameter('generateImageList', '', @(x) ischar(x));
% job related parameters
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 1, @isnumeric);
ip.addParameter('minModifyTime', 1, @isnumeric);
ip.addParameter('maxModifyTime', 10, @isnumeric);
ip.addParameter('maxWaitLoopNum', 10, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root of XR_Repository. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../'];
cd(repo_rt);

pr = ip.Results;
overwrite = pr.overwrite;
streaming = pr.streaming;
skewAngle = pr.skewAngle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
objectiveScan = pr.objectiveScan;
zStageScan = pr.zStageScan;
reverse = pr.reverse;
channelPatterns = pr.channelPatterns;
save16bit = pr.save16bit;
resampleType = pr.resampleType;
resampleFactor = pr.resampleFactor;
dzFromEncoder = pr.dzFromEncoder;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
save3DStack = pr.save3DStack;
%deskew and rotate
deskew = pr.deskew;
rotate = pr.rotate;
parseSettingFile = pr.parseSettingFile;
flipZstack = pr.flipZstack;
DSRCombined = pr.DSRCombined;
FFCorrection = pr.FFCorrection;
BKRemoval = pr.BKRemoval;
lowerLimit = pr.lowerLimit;
constOffset = pr.constOffset;
FFImagePaths = pr.FFImagePaths;
backgroundPaths = pr.backgroundPaths;
inputBbox = pr.inputBbox;
% stitch parameters
stitch = pr.stitch;
stitchResultDirName = pr.stitchResultDirName;
imageListFullpaths = pr.imageListFullpaths;
axisOrder = pr.axisOrder;
blendMethod = pr.blendMethod;
xcorrShift = pr.xcorrShift;
xcorrMode = pr.xcorrMode;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
edgeArtifacts = pr.edgeArtifacts;
primaryCh = pr.primaryCh;
stitchMIP = pr.stitchMIP;
onlineStitch = pr.onlineStitch;
generateImageList = pr.generateImageList;
% job related
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
parseCluster = pr.parseCluster;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
minModifyTime = pr.minModifyTime;
maxModifyTime = pr.maxModifyTime;
maxWaitLoopNum = pr.maxWaitLoopNum;
mccMode = pr.mccMode;
configFile = pr.configFile;

% suppress directory exists warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');

if ischar(dataPaths)
    dataPaths = {dataPaths};
end
dataPaths = unique(dataPaths);

if ischar(imageListFullpaths)
    imageListFullpaths = {imageListFullpaths};
end

nd = numel(dataPaths);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), '/')
        dataPaths{d} = [dataPath, '/'];
    end
end
% add support for predefined folders for real time processing
if streaming
    for d = 1 : nd
        dataPath = dataPaths{d};
        if ~exist(dataPath, 'dir')
            mkdir(dataPath);
            fileattrib(dataPath, '+w', 'g');
        end
    end
end

% change separator backslash to slash
if ispc
    dataPaths = cellfun(@(x) strrep(x, '\', '/'), dataPaths, 'unif', 0);
    imageListFullpaths = cellfun(@(x) strrep(x, '\', '/'), imageListFullpaths, 'unif', 0);
end

if numel(overwrite) == 1
    overwrite = repmat(overwrite, 1, 3);
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath, jobLogDir);
end

% save zarr 
ext = '.tif';
if saveZarr
    ext = '.zarr';
end

% for stitching, enable DS and DSR
% For multiple datasets, if the image list is not provide for any dataset,
% set stitch as false. 
if stitch
    if ~streaming
        if numel(imageListFullpaths) ~= nd
            warning("The number of image list files does not match that of the data paths, please make sure image list files are provided for each dataset!")
            stitch = false;
        else
            for d = 1 : nd
                if ~exist(imageListFullpaths{d}, 'file')
                    warning('Image list filename %s does not exist, set stitch option as False', imageListFullpaths{d});
                    stitch = false;
                    break;
                end
            end
        end
    end
    
    deskew = true;
    rotate = true;
    stchPaths = cell(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};
        if isempty(stitchResultDirName)
            if strcmp(blendMethod, 'none')
                stitchResultDirName = 'matlab_stitch';
            else
                stitchResultDirName = sprintf('matlab_stitch_%s', blendMethod);                
            end
        end
        stchPath = [dataPath, '/', stitchResultDirName, '/'];
        if overwrite(3) && exist(stchPath, 'dir')
            rmdir(stchPath, 's');
        end
        if ~exist(stchPath, 'dir')
            mkdir(stchPath);
            fileattrib(stchPath, '+w', 'g');
        end
        stchPaths{d} = stchPath;
    end
    
    % check if axis order is valid
    axisOrder = strrep(axisOrder, ' ', '');
    pattern = '^(-?x,?-?y,?-?z|-?y,?-?x,?-?z|-?z,?-?y,?-?x|-?x,?-?z,?-?y|-?x,?-?z,?-?y|-?y,?-?z,?-?x)$';
    if ~regexpi(axisOrder, pattern)
        error("The axisOrder is not right, it must has the form like 'y,x,z' or '-x,y,z' (flipped in x-axis)!");
    end
end

% first check if DS and Deconvolution directories exist
if ~rotate
    DSRCombined = false;
end

if DSRCombined
    deskew = true;
    rotate = true;
end

if deskew
    if ~DSRCombined
        dsPaths = cell(nd, 1);
        for d = 1 : nd
            dataPath = dataPaths{d};
            dsPath = [dataPath, '/DS/'];
            if overwrite(1) && exist(dsPath, 'dir')
                rmdir(dsPath, 's');
            end
            if ~exist(dsPath, 'dir')
                mkdir(dsPath);
                fileattrib(dsPath, '+w', 'g');
            end
            dsPaths{d} = dsPath;
        end
    end
    
    % check LS and Background images for flat field correction
    if FFCorrection
        if numel(channelPatterns) ~= numel(FFImagePaths) || numel(channelPatterns) ~= numel(backgroundPaths) 
            error('The number of channels in channelPatterns does not match the number of LS files or Background Files!')
        end
        for c = 1 : numel(FFImagePaths)
            if ~exist(FFImagePaths{c}, 'file')
                error('LS Image file %s does not exist!', FFImagePaths{c});
            end
            if ~exist(backgroundPaths{c}, 'file')
                error('LS Image file %s does not exist!', backgroundPaths{c});
            end
        end
    end
else
    rotate = false;
    DS = false;
    DSR = false;
end
    
if rotate
    dsrPaths = cell(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};
        dsrPath = [dataPath, '/DSR/'];
        if overwrite(2) && exist(dsrPath, 'dir')
            rmdir(dsrPath, 's');
        end
        if ~exist(dsrPath, 'dir')
            mkdir(dsrPath);
            fileattrib(dsrPath, '+w', 'g');
        end
        dsrPaths{d} = dsrPath;
    end
end

% define resampleFactor based on resampleFactor type, resampleFactor parameter and the microscope setting
if ~isempty(resampleFactor)
    resampleType = 'given';
end
[resampleFactor, zAniso] = XR_checkResampleSetting(resampleType, resampleFactor, objectiveScan, skewAngle, xyPixelSize, dz);

% get actual dz from encoder positions
if dzFromEncoder && ~streaming
    dz_all = zeros(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};        
        dz_actual = XR_estimate_actual_step_size_from_encoder(dataPath, 'dz', dz);
        dz_all(d) = dz_actual;
    end   
else
    dz_all = repmat(dz, nd, 1);
end


%% check existing files and parse channels
Decon = false;
deconPaths = {};
[fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, latest_modify_times, ...
    FTP_inds, maskFullpaths] = XR_parseImageFilenames(dataPaths, channelPatterns, ...
    parseSettingFile, flipZstack, Decon, deconPaths, streaming, minModifyTime, zarrFile);

nF = numel(fnames);

% flags: for thee: deskew w/o rotate, decon w/o rotate, rotate
is_done_flag = false(nF, 2);
% ensure only one stitch wrapper is running of a dataset
stitch_running = false(nd, 1); 
trial_counter = zeros(nF, 2);
waitLoopCounter = 0;

% For no-streaming computing, first skip the disabled options. 
if ~streaming
    if ~stitch
        is_done_flag(:, 2) = true;
    end
end

if parseCluster
    job_ids = -ones(nF, 2);
    imSize_mat = zeros(nF, 3, 2);
    dataSize_mat = zeros(nF, 2);
    dataSize_mat(:, 1) = dataSizes;
end

% use while loop to perform computing for all images
ts = tic;
while ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') || ...
        (streaming && (nF == 0 || any(latest_modify_times < maxModifyTime) || waitLoopCounter < maxWaitLoopNum))
    for f = 1 : nF
        tic
        if all(is_done_flag(f, :))
            continue;
        end
        
        % first deskew and rotation
        fname = fnames{f};
        [~, fsname] = fileparts(fname);
        fdind = fdinds(f);
        partialvol = partialvols(f);
        gfname = gfnames{f};
        dataPath = dataPaths{fdind};
        dz_f = dz_all(fdind);
        
        frameFullpath = [dataPath, fname];
        % check wheter the file is deleted during the computing.
        if ~exist(frameFullpath, 'file')
            is_done_flag(f, :) = true;
            continue
        end
        
        task_id = rem(f, 5000);
        
        %% deskew w/o rotate
        if deskew
            if ~DSRCombined
                dsPath = dsPaths{fdind};
                dsFullpath = [dsPath, fsname, ext];
                tmpFullpath = sprintf('%s/%s.tmp', dsPath, fsname);
            end
            if rotate
                dsrPath = dsrPaths{fdind};
                dsrFullpath = [dsrPath, fsname, ext];
                tmpFullpath = sprintf('%s/%s.tmp', dsrPath, fsname);                
            end

            if (DSRCombined || exist(dsFullpath, 'file') || exist(dsFullpath, 'dir')) && (~rotate || exist(dsrFullpath, 'file') || exist(dsrFullpath, 'dir'))
                is_done_flag(f, 1) = true;
                if exist(tmpFullpath, 'file')
                    delete(tmpFullpath);
                end
            end
        else
            is_done_flag(f, 1) = true;
        end
        
        if ~is_done_flag(f, 1) 
            if FFCorrection
                % LLFFMapping =  ~cellfun(@isempty, regexpi(fname, channelPatterns));
                % change to contains.m to unify the matching
                LLFFMapping =  cellfun(@(x) contains(frameFullpath, x), channelPatterns);
                FFImage = FFImagePaths{LLFFMapping};
                BackgroundImage = backgroundPaths{LLFFMapping};
            else
                FFImage = '';
                BackgroundImage = '';
            end
            
            flipZstack = flipZstack_mat(f);
            
            % set up input file for either single volume file or a
            % group of files
            ds_input_path = {frameFullpath};
            if partialvol
                ds_input_path = cellfun(@(x) [dataPath, x], gfname, 'unif', 0);
            end
            ds_input_str = sprintf('{''%s''}', strjoin(ds_input_path, ''','''));

            func_str = sprintf(['XR_deskewRotateFrame(%s,%.20d,%.20d,''skewAngle'',%.20d,', ...
                '''objectiveScan'',%s,''zStageScan'',%s,''reverse'',%s,''FFCorrection'',%s,', ...
                '''BKRemoval'',%s,''lowerLimit'',%.20d,''constOffset'',[%s],''FFImage'',''%s'',', ...
                '''BackgroundImage'',''%s'',''rotate'',%s,''resampleFactor'',[%s],''DSRCombined'',%s,', ...
                '''inputBbox'',%s,''flipZstack'',%s,''save16bit'',%s,''save3DStack'',%s,''saveZarr'',%s,', ...
                '''uuid'',''%s'')'], ds_input_str, xyPixelSize, dz_f, skewAngle, string(objectiveScan), ...
                string(zStageScan), string(reverse),  string(FFCorrection), string(BKRemoval), ...
                lowerLimit, num2str(constOffset, '%0.10f'),  FFImage, BackgroundImage, string(rotate), ...
                strrep(num2str(resampleFactor, '%d,'), ' ', ''),  string(DSRCombined), strrep(mat2str(inputBbox), ' ', ','), ...
                string(flipZstack), string(save16bit(1)), string(save3DStack), string(saveZarr), uuid);
            
            if exist(tmpFullpath, 'file') || parseCluster
                if parseCluster
                    if ~DSRCombined
                        estRequiredMemory = XR_estimateComputingMemory(frameFullpath, 'steps', {'deskew'});
                    else
                        estRequiredMemory = dataSize_mat(f, 1) / 2^30 * 2 * (7 + 3 / prod(resampleFactor) + trial_counter(f, 1) * 8);
                    end
                    memAllocate = estRequiredMemory;

                    job_id = job_ids(f, 1);
                    [job_id, ~, submit_status] = generic_single_job_submit_wrapper(func_str, job_id, task_id, ...
                        'jobLogFname', job_log_fname, 'jobErrorFname', job_log_error_fname, ...
                        masterCompute=masterCompute, cpusPerTask=cpusPerTask, memAllocate=memAllocate, ...
                        mccMode=mccMode, configFile=configFile);

                    job_ids(f, 1) = job_id;
                    trial_counter(f, 1) = trial_counter(f, 1) + submit_status;
                else
                    temp_file_info = dir(tmpFullpath);
                    if minutes(datetime('now') - temp_file_info.date) < unitWaitTime
                        continue; 
                    else
                        fclose(fopen(tmpFullpath, 'w'));
                    end
                end
            else
                fclose(fopen(tmpFullpath, 'w'));
            end
            if ~parseCluster
                tic; feval(str2func(['@()', func_str])); toc;
                fprintf('\n');
                trial_counter(f, 1) = trial_counter(f, 1) + 1;
            end

            % check if computing is done
            if (DSRCombined || exist(dsFullpath, 'file') || exist(dsFullpath, 'dir')) && (~rotate || exist(dsrFullpath, 'file') || exist(dsrFullpath, 'dir'))
                is_done_flag(f, 1) = true;
                if exist(tmpFullpath, 'file')
                    delete(tmpFullpath);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% then stitching
        if stitch
            if is_done_flag(f, 2)
                continue;
            end
            stchPath = stchPaths{fdind};
            % dir_info = dir(sprintf('%s/%s*Abs.tif', stchPath, fname(regexp(fname, 'Scan.*') : end - 43)));
            stch_dir_info = dir(sprintf('%s/%s*Abs.zarr', stchPath, fsname(1 : end - 39)));

            if ~isempty(stch_dir_info) 
                stch_fname = stch_dir_info(1).name;
                stchFullpath = [stchPath, stch_fname];
                stchTmpFullpath = sprintf('%s.tmp', stchFullpath(1 : end - 4));
                if exist(stchFullpath, 'file') || exist(stchFullpath, 'dir')
                    is_done_flag(f, 2) = true;
                    if exist(stchTmpFullpath, 'file')
                        delete(stchTmpFullpath);
                    end
                end
            end
            if streaming && ~isempty(generateImageList) 
                stitch_flags_d = is_done_flag(:, 2) == 0 & fdinds == fdind;
                if any(stitch_flags_d) && (f == find(stitch_flags_d, 1, 'first'))
                    switch generateImageList
                        case 'from_encoder'
                            imageListFullpaths{fdind} = stitch_generate_imagelist_from_encoder(dataPath, dz_f, channelPatterns);
                        case 'from_sqlite'
                            imageListFullpaths{fdind} = stitch_generate_imagelist_from_sqlite(dataPath);                        
                    end
                end
            end
        else
            is_done_flag(f, 2) = true;
        end
        
        if ~is_done_flag(f, 2)
            % if the deskew result is not available, wait for the deskew data to be finished        
            if ~parseCluster && (~exist(frameFullpath, 'file') || ~exist(dsrFullpath, 'file'))
                % if DS or DSR, it means the deskew and rotation not done.
                continue;
            end
            
            stitch_DS = false;
            stitch_DSR = true;
            ProcessedDirStr = 'DSR';
            resampleType = 'isotropic';
            bbox = [];
            imageListFullpath = imageListFullpaths{fdind};
            channelPatterns_str = sprintf('{''%s''}', strjoin(channelPatterns, ''','''));
            if ~parseCluster
                mccMode = false;
                maxTrialNum = 1;
            end
            stitchZarrFile = saveZarr;

            func_str = sprintf(['XR_matlab_stitching_wrapper(''%s'',''%s'',''ResultDir'',''%s'',', ...
                '''streaming'',%s,''DS'',%s,''DSR'',%s,''channelPatterns'',%s,''ProcessedDirStr'',''%s'',', ...
                '''axisOrder'',''%s'',''resampleType'',''%s'',''resampleFactor'',[%s],''reverse'',%s,', ...
                '''parseSettingFile'',%s,''xcorrShift'',%s,''xcorrMode'',''%s'',''xyMaxOffset'',%.10f,', ...
                '''zMaxOffset'',%.10f,''edgeArtifacts'',%d,''blendMethod'',''%s'',', ...
                '''outBbox'',[%s],''save16bit'',%s,''primaryCh'',''%s'',''stitchMIP'',%s,''onlineStitch'',%s,', ...
                '''zarrFile'',%s,''maxTrialNum'',%d,''mccMode'',%s,''configFile'',''%s'',''uuid'',''%s'')'], ...
                dataPath, imageListFullpath, stitchResultDirName, string(parseCluster & streaming), string(stitch_DS), ...
                string(stitch_DSR), channelPatterns_str, ProcessedDirStr, axisOrder, resampleType, ...
                strrep(num2str(resampleFactor, '%.10d,'), ' ', ''), string(reverse), string(parseSettingFile), ...
                string(xcorrShift), xcorrMode, xyMaxOffset, zMaxOffset, edgeArtifacts, blendMethod, ...
                strrep(num2str(bbox, '%d,'), ' ', ''), string(save16bit(2)), primaryCh, ...
                strrep(mat2str(stitchMIP), ' ', ','), string(onlineStitch), ...
                string(stitchZarrFile), maxTrialNum, string(mccMode), configFile, uuid);

            if parseCluster
                dfirst_ind = find(fdinds == fdind, 1, 'first');
                job_id = job_ids(dfirst_ind, 2);
                task_id = rem(dfirst_ind, 5000);
                memFactor = 20;
                if any(stitchMIP)
                    memFactor = 2;
                end
                estRequiredMemory = dataSize_mat(f, 1) / 2^30 * 2 * memFactor;
                memAllocate = estRequiredMemory;

                [job_id, ~, submit_status] = generic_single_job_submit_wrapper(func_str, job_id, task_id, ...
                    'jobLogFname', job_log_fname, 'jobErrorFname', job_log_error_fname, ...
                    masterCompute=masterCompute, cpusPerTask=cpusPerTask, memAllocate=memAllocate, ...
                    mccMode=mccMode, configFile=configFile);

                job_ids(dfirst_ind, 2) = job_id;
                job_ids(f, 2) = job_id;
                trial_counter(dfirst_ind, 2) = trial_counter(dfirst_ind, 2) + submit_status;
                trial_counter(f, 2) = trial_counter(f, 2) + submit_status;
            else
                tic; feval(str2func(['@()', func_str])); toc;
                fprintf('\n');                
            end
            
            if ~isempty(stch_dir_info) && exist(stchFullpath, 'file')
                is_done_flag(f, 2) = true;
                if exist(stchTmpFullpath, 'file')
                    delete(stchTmpFullpath);
                end
            end
        end
    end
    
    %% wait for running jobs finishing and checking for new coming images
    nF_done = sum(all(is_done_flag, 2));
    sprintf('Time %d s: %d / %d (%0.3f) are finished!\n', toc(ts), nF_done, nF, nF_done / nF);
    
    if ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
        waitLoopCounter = 0;
        pause(30);
        % if not streaming, do not check for new files. 
        if ~streaming
            continue;
        end
    else
         % if not streaming, exit when all existing files are processed. 
        if streaming 
            if waitLoopCounter < maxWaitLoopNum
                waitLoopCounter = waitLoopCounter + 1;
                pause(30);
            end
        else
            break;
        end
    end
            
    % check whether there are new coming images (only for streaming option)
    [fnames_new, fdinds_new, gfnames_new, partialvols_new, dataSizes_new, flipZstack_mat_new, ...
        latest_modify_times_new, FTP_inds_new, maskFullpaths_new] = XR_parseImageFilenames(dataPaths, ...
        channelPatterns, parseSettingFile, flipZstack, Decon, deconPaths, streaming, minModifyTime, zarrFile);
    if isempty(fnames_new)
        continue;
    end
    
    new_inds = cellfun(@(x) ~any(contains(fnames, x)), fnames_new);
    if ~any(new_inds)
        continue;
    end
    
    % add new files and their computing flags
    cur_fnames = fnames_new(new_inds);
    cur_fdinds = fdinds_new(new_inds);
    cur_gfnames = gfnames_new(new_inds);
    cur_partialvols = partialvols_new(new_inds);
    cur_dataSizes = dataSizes_new(new_inds);
    cur_flipZstack_mat = flipZstack_mat_new(new_inds);
       
    nFnew = numel(cur_fnames);
    fnames(end + 1 : end + nFnew) = cur_fnames;
    fdinds = cat(1, fdinds, cur_fdinds);
    partialvols = cat(1, partialvols, cur_partialvols);
    gfnames(end + 1 : end + nFnew) = cur_gfnames;
    flipZstack_mat = cat(1, flipZstack_mat(:), cur_flipZstack_mat(:));
    
    latest_modify_times = max(latest_modify_times, latest_modify_times_new);
    
    nF = numel(fnames);
    is_done_flag = cat(1, is_done_flag, false(nFnew, 4));
    trial_counter = cat(1, trial_counter, zeros(nFnew, 4));
    
    if parseCluster
        job_ids = cat(1, job_ids, -ones(nFnew, 4));
        dataSize_mat = cat(1, dataSize_mat, [cur_dataSizes, zeros(nFnew, 1)]);
    end
end

end

