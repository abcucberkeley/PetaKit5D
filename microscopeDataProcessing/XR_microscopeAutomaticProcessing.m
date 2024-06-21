function [] = XR_microscopeAutomaticProcessing(dataPaths, varargin)
% automatic image processing pipeline for microscopy data. It can perform
% deskew w/o rotation, deconvolution, rotation afte deconvolution. It supports 
% both existing data in the directory or real-time transferring of data from
% other sources, i.e., microscope. Once all computing done, it will wait another 
% 5 minutes for potential new coming files. 
%
% Inputs :   
%            dataPath : directory paths for the data sets. Either a string
%                       for a single dataset or a cell array of paths for several
%                       datasets with same experiment settings. 
%
%
% Options (as 'specifier'-value pairs): 
%
%         'overwrite' : true|{false}, or a length 5 bool vector. overwrite existing results.
%         'streaming' : {true}|false. True for real-time processing, and false for existing data
%   'channelPatterns' : An cell array of channel identifies and orders for included channels.
%                       Supported formats: Cam[A/B]_ch[0-9] or ch[0-9].
%          'Channels' : Channel wavelength (currently not required).
%         'skewAngle' : skew angle of the stage. Default: 32.45.
%                'dz' : stage scan interval. Default: 0.5.
%       'xyPixelSize' : pixel size. Default: 0.108.
%           'reverse' : true|{false}. Inverse direction of z axis. 
%     'objectiveScan' : true|{false}. Objective scan. Default: 5.
%   'sCMOSCameraFlip' : true|{false}. sCMOS camera flip. 
%         'save16bit' : 1 X 4 bool vector. Save 16bit result for deskew/rotate, stitch, decon, and rotate after decon. 
%            'deskew' : {true}|false. deskew the data.
%            'rotate' : {true}|false. rotate deskewed data.
%            'stitch' : true|{false}. stitch deskewed and rotated data.
%             'Decon' : {true}|false. Deconvolution on data.
%  'rotateAfterDecon' : true|{false}. rotate deconvolution results using matlab rotation function.
%    'FFCorrection' : true|{false}. Flat-field correction.
%        'lowerLimit' : A number between 0 and 1. Lower limit of intensity range for binarization. 
%         'cudaDecon' : {true}|false. Use cudaDecon for deconvolution. If there is no GPU, false by default. 
%      'FFImagePaths' : An array of full paths of flag field, with channel orders defined in 'channelPatterns'.
%   'backgroundPaths' : An array of full paths of camera dark current images, with channel orders defined in 'channelPatterns'.
%   'stitchResultDirName' : matlab stitch directory name. 
%'imageListFullpaths' : A cell array of full paths. Image list csv file for tile coordinates for each dataset. 
%         'axisOrder' : Axis order for the coordinates. Format: 'xyz', '-x,y,z', 'y,x,z' etc.  
%       'blendMethod' : Method to handle overlap regions. Options: 'none' (default), 'mean', 'median' and 'max'. 
%        'xcorrShift' : {true}|false. Use cross-correlation based registration in stitching
%         'xcorrMode' : Reference registration mode. 'primary': primary channel; 
%                       'primaryfirst': first time point in primary channel; 'all': its own xcorr. 
%         'cudaDecon' : true|{false}. Use cuda decon method for deconvolution.
%          'cppDecon' : true|{false}. Use cpp decon method for deconvolution
%      'cppDeconPath' : Program path for cpp decon package. 
%       'loadModules' : Cammands to load dependency modules for cpp decon package (only for ABC cluster).
%     'cudaDeconPath' : Program path for cuda decon package. 
%        'OTFGENPath' : Program path for otf generation function in cuda decon package. 
%                'DS' : true|{false}. Use deskewed data for deconvolution.
%               'DSR' : {true}|false. Use deskewed rotated data for deconvolution. 
%        'Background' : Background intensity for deconvolution. Default: 99 if not provided.
%       'EdgeErosion' : Number of voxels from edges (defined in raw data) to 
%                       exclude in decon in order to remove edge artifacts. Default: 8. 
%        'ErodeByFTP' : {true}|false. Use the first time point and first
%                       channel to define an eroded mask to remove edge artifact. 
%       'deconrotate' : true|{false}. rotate deconvolution result within deconvolution steps.
%      'psfFullpaths' : Full paths of psf files. A file for a channel. The order is the same as channelPatterns.
%      'parseCluster' : Use slurm-based cluster computing.
%         'jobLogDir' : Log directory for the slurm jobs.
%       'cpusPerTask' : Number of cpus for a job. Default: 2
%      'cpuOnlyNodes' : {true}|false. Use CPU-only nodes in ABC cluster. 
%              'uuid' : unique string for a job for saving files. 
%       'maxTrialNum' : Max number of times to rerun failure cases. 
%      'unitWaitTime' : The wait time per file in minutes to check whether the computing is done.
%     'minModifyTime' : The minimum time in minutes for the latest modified file to decide whether it is fully transferred.
%     'maxModifyTime' : The maximum time in minutes to check whether there are coming new files.
%    'maxWaitLoopNum' : Number of maximum loops without any computing. Default: 10
%
%
% Author: Xiongtao Ruan (03/2020)
% 
% xruan (03/18/2020): add option of rotation after deconvolution; add
% option for skipping deskew.
% xruan (04/16/2020): add option for CPP decon
% xruan (07/13/2020): add option for matlab stitching, and update framework
%                     structure for each component
% xruan (07/17/2020): add overwrite option
% xruan (07/18/2020): add option for support of multiple datasets (with same settings)
% xruan (07/19/2020): add streaming option, set Channel pattern to filter channels to process
% xruan (07/23/2020): add CPU node only option, add mask file option for
% the the first time point for edge erosion in Decon
% xruan (07/28/2020): fix issue in decon, when ErodeByFTP true and the deconm
%                     result for the first time point exist but the mask file 
%                     doesn't exist, in which it will stuck. 
%                     Also check if the image size match the mask size, if
%                     not use its own mask for edge erosion.
% xruan (08/03/2020): add axis-order for the stitching
% xruan (08/07/2020): change file attributes for folders so that group users have write access
% xruan (08/08/2020): add more memory for per-cpu, fix issue for support of 
%                     multiple datasets in stitching; set xcorrShift as
%                     false by default
% xruan (08/21/2020): add support for Objective Scan
% xruan (08/26/2020): add support for user defined output resolution (dsr,
% stitch, decon) （to be done)
% xruan (10/05/2020): add support for zarr-stitching pipeline with option:
%                     'zarr' ('matlab' for original pipeline)
% xruan (10/06/2020): add support for input of multiple parts of same volume
% xruan (10/15/2020): add support for combined process for DS and rotate;
%                     also add support for only processing first timepoint
% xruan (10/20/2020): simplify image size check to data size
% xruan (11/02/2020): add option for rotated PSF
% xruan (12/05/2020): add option to remove background 
%                     add option to flip z stack in raw data (for negative X interval)
% xruan (12/18/2020): simplify code for processing functions; add support
%                     for not save 3d stack (for quick checking of results)
% xruan (01/13/2021): add support of resampleFactor for DSR and following analysis
% xruan (06/10/2021): add support for threshold and debug mode in matlab decon simplified version. 
% xruan (06/11/2021): add support for gpu computing for chuck decon in matlab decon wrapper
% xruan (07/05/2021): add support for user defined resampleFactor (arbitary factor)
% xruan (07/27/2021): add support for z-stage scan for deskew
% xruan (09/13/2021): add support for calcualating actual dz from encoder positions
% xruan (10/21/2021): add support for zarr file as input
% xruan (01/25/2022): add support for bbox crop before processing



ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x)); % data structure from loadConditionData
ip.addParameter('overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 5) && islogical(x));
ip.addParameter('streaming', true,  @islogical); % if true, check for new files. If false, assume all files transferred completely.
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('skewAngle', 32.45, @isscalar);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('reverse', true, @islogical);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('zStageScan', false, @islogical);
ip.addParameter('save16bit', [true, true, true, true], @(x) (numel(x) == 1 || numel(x) == 4) && islogical(x));
ip.addParameter('onlyFirstTP', false, @islogical);
ip.addParameter('dzFromEncoder', false, @islogical);
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('saveZarr', false, @islogical); % use zarr file as output
ip.addParameter('save3DStack', true , @islogical); % option to save 3D stack or not
% pipeline steps
ip.addParameter('deskew', true, @islogical);
ip.addParameter('rotate', true, @islogical);
ip.addParameter('stitch', false, @islogical);
% deskew and rotation options
ip.addParameter('parseSettingFile', false, @islogical); % use setting file to decide whether filp Z stack or not.
ip.addParameter('flipZstack', false, @islogical); % 
ip.addParameter('DSRCombined', true, @islogical); 
ip.addParameter('FFCorrection', false, @islogical);
ip.addParameter('BKRemoval', false, @islogical);
ip.addParameter('lowerLimit', 0.4, @isnumeric); % this value is the lowest
ip.addParameter('constOffset', [], @(x) isnumeric(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('FFImagePaths', {'','',''}, @iscell);
ip.addParameter('backgroundPaths', {'','',''}, @iscell);
ip.addParameter('resampleType', 'isotropic', @ischar); % resampleFactor type: given, isotropic, xy_isotropic
ip.addParameter('resampleFactor', [], @isnumeric); % resampleFactor
ip.addParameter('inputBbox', [], @isnumeric); % bbox for input in deskew and rotate
% stitch parameters
ip.addParameter('stitchPipeline', 'zarr', @ischar); % matlab or zarr
ip.addParameter('stitchResultDirName', '', @ischar);
ip.addParameter('imageListFullpaths', '', @(x) ischar(x) || iscell(x));
ip.addParameter('axisOrder', 'xyz', @(x) ischar(x));
ip.addParameter('blendMethod', 'none', @ischar);
ip.addParameter('xcorrShift', false, @islogical);
ip.addParameter('xcorrMode', 'primaryFirst', @(x) ismember(lower(x), {'primary', 'primaryfirst', 'all'})); % 'primary': choose one channel as primary channel, 
ip.addParameter('xyMaxOffset', 300, @isnumeric); % max offsets in xy axes
ip.addParameter('zMaxOffset', 50, @isnumeric); % max offsets in z axis
ip.addParameter('edgeArtifacts', 2, @isnumeric);
ip.addParameter('primaryCh', '', @ischar);
ip.addParameter('stitchMIP', [], @(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3))); % 1x3 vector or vector, by default, stitch MIP-z
ip.addParameter('onlineStitch', false, @(x) islogical(x)); % support for online stitch (with partial number of tiles). 
ip.addParameter('generateImageList', '', @(x) ischar(x)); % for real time processing, {'', 'from_encoder', 'from_sqlite'}
% job related parameters
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 1, @isnumeric);
ip.addParameter('minModifyTime', 1, @isnumeric); % the minimum duration of last modify time of a file, in minute.
ip.addParameter('maxModifyTime', 10, @isnumeric); % the maximum duration of last modify time of a file, in minute.
ip.addParameter('maxWaitLoopNum', 10, @isnumeric); % the max number of loops the loop waits with all existing files processed. 
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);
ip.addParameter('GPUConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root of XR_Repository. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../'];
cd(repo_rt);

pr = ip.Results;
overwrite = pr.overwrite;
streaming = pr.streaming;
% Resolution = pr.Resolution;
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
save3DStack = pr.save3DStack; % only for DS and DSR for now
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
onlyFirstTP = pr.onlyFirstTP;
primaryCh = pr.primaryCh;
stitchMIP = pr.stitchMIP;
onlineStitch = pr.onlineStitch;
generateImageList = pr.generateImageList;
% job related
largeFile = pr.largeFile;
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
GPUConfigFile = pr.GPUConfigFile;

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
    overwrite = repmat(overwrite, 1, 5);
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
        FTP_ind = FTP_inds(fdind);
        
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
                '''inputBbox'',%s,''flipZstack'',%s,''save16bit'',%s,''save3DStack'',%s,''saveZarr'',%s)'], ...
                ds_input_str, xyPixelSize, dz_f, skewAngle, string(objectiveScan), string(zStageScan), ...
                string(reverse),  string(FFCorrection), string(BKRemoval), lowerLimit, ...
                num2str(constOffset, '%0.10f'),  FFImage, BackgroundImage, string(rotate), ...
                strrep(num2str(resampleFactor, '%d,'), ' ', ''),  string(DSRCombined), strrep(mat2str(inputBbox), ' ', ','), ...
                string(flipZstack), string(save16bit(1)), string(save3DStack), string(saveZarr));
            
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
                        memAllocate=memAllocate, mccMode=mccMode, configFile=configFile);

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
                if exist(stchFullpath, 'file') || (strcmp(stitchPipeline, 'zarr') && exist(stchFullpath, 'dir'))
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
            % onlyFirstTP = false;
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
                '''onlyFirstTP'',%s,''outBbox'',[%s],''save16bit'',%s,', ...
                '''primaryCh'',''%s'',''stitchMIP'',%s,''onlineStitch'',%s,', ...
                '''zarrFile'',%s,''maxTrialNum'',%d,''mccMode'',%s,''configFile'',''%s'')'], ...
                dataPath, imageListFullpath, stitchResultDirName, string(parseCluster & streaming), string(stitch_DS), ...
                string(stitch_DSR), channelPatterns_str, ProcessedDirStr, axisOrder, resampleType, ...
                strrep(num2str(resampleFactor, '%.10d,'), ' ', ''), string(reverse), string(parseSettingFile), ...
                string(xcorrShift), xcorrMode, xyMaxOffset, zMaxOffset, edgeArtifacts, blendMethod, ...
                strrep(num2str(bbox, '%d,'), ' ', ''), string(save16bit(2)), primaryCh, ...
                strrep(mat2str(stitchMIP), ' ', ','), string(onlineStitch), ...
                string(stitchZarrFile), maxTrialNum, string(mccMode), configFile);

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
                    memAllocate=memAllocate, mccMode=mccMode, configFile=configFile);

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
    FTP_inds = FTP_inds_new;
    
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

