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
%         'Overwrite' : true|{false}, or a length 5 bool vector. Overwrite existing results.
%         'Streaming' : {true}|false. True for real-time processing, and false for existing data
%   'ChannelPatterns' : An cell array of channel identifies and orders for included channels.
%                       Supported formats: Cam[A/B]_ch[0-9] or ch[0-9].
%          'Channels' : Channel wavelength (currently not required).
%         'SkewAngle' : skew angle of the stage. Default: 32.45.
%                'dz' : stage scan interval. Default: 0.5.
%       'xyPixelSize' : pixel size. Default: 0.108.
%           'Reverse' : true|{false}. Inverse direction of z axis. 
%     'ObjectiveScan' : true|{false}. Objective scan. Default: 5.
%   'sCMOSCameraFlip' : true|{false}. sCMOS camera flip. 
%         'Save16bit' : 1 X 4 bool vector. Save 16bit result for deskew/rotate, stitch, decon, and rotate after decon. 
%            'Deskew' : {true}|false. Deskew the data.
%            'Rotate' : {true}|false. Rotate deskewed data.
%            'Stitch' : true|{false}. Stitch deskewed and rotated data.
%             'Decon' : {true}|false. Deconvolution on data.
%  'RotateAfterDecon' : true|{false}. Rotate deconvolution results using matlab rotation function.
%    'LLFFCorrection' : true|{false}. Flat-field correction.
%        'LowerLimit' : A number between 0 and 1. Lower limit of intensity range for binarization. 
%         'cudaDecon' : {true}|false. Use cudaDecon for deconvolution. If there is no GPU, false by default. 
%      'LSImagePaths' : An array of full paths of flag field, with channel orders defined in 'ChannelPatterns'.
%   'BackgroundPaths' : An array of full paths of camera dark current images, with channel orders defined in 'ChannelPatterns'.
%   'stitchResultDir' : matlab stitch directory name. 
%'imageListFullpaths' : A cell array of full paths. Image list csv file for tile coordinates for each dataset. 
%         'axisOrder' : Axis order for the coordinates. Format: 'xyz', '-x,y,z', 'y,x,z' etc.  
%       'BlendMethod' : Method to handle overlap regions. Options: 'none' (default), 'mean', 'median' and 'max'. 
%        'xcorrShift' : {true}|false. Use cross-correlation based registration in stitching
%         'xcorrMode' : Reference registration mode. 'primary': primary channel; 
%                       'primaryfirst': first time point in primary channel; 'all': its own xcorr. 
%      'boundboxCrop' : Bounding box to crop ROI after stitching. Empty or [y_min, x_min, z_min, y_max, x_max, z_max]. 
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
%       'deconRotate' : true|{false}. Rotate deconvolution result within deconvolution steps.
%      'psfFullpaths' : Full paths of psf files. A file for a channel. The order is the same as ChannelPatterns.
%      'parseCluster' : Use slurm-based cluster computing.
%         'jobLogDir' : Log directory for the slurm jobs.
%       'cpusPerTask' : Number of cpus for a job. Default: 2
%      'cpuOnlyNodes' : {true}|false. Use CPU-only nodes in ABC cluster. 
%              'uuid' : unique string for a job for saving files. 
%       'maxTrialNum' : Max number of times to rerun failure cases. 
%      'unitWaitTime' : The wait time per file in minutes to check whether the computing is done.
%     'minModifyTime' : The minimum time in minutes for the latest modified file to decide whether it is fully transferred.
%     'maxModifyTime' :  The maximum time in minutes to check whether there are coming new files.
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
% xruan (07/17/2020): add Overwrite option
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
% xruan (01/13/2021): add support of resample for DSR and following analysis
% xruan (06/10/2021): add support for threshold and debug mode in matlab decon simplified version. 
% xruan (06/11/2021): add support for gpu computing for chuck decon in matlab decon wrapper
% xruan (07/05/2021): add support for user defined resample (arbitary factor)
% xruan (07/27/2021): add support for z-stage scan for Deskew
% xruan (09/13/2021): add support for calcualating actual dz from encoder positions
% xruan (10/21/2021): add support for zarr file as input
% xruan (01/25/2022): add support for bbox crop before processing



ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths'); % data structure from loadConditionData
ip.addParameter('Overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 5) && islogical(x));
ip.addParameter('Streaming', true,  @islogical); % if true, check for new files. If false, assume all files transferred completely.
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('Channels', [488, 560, 642], @isnumeric);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('Reverse', true, @islogical);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('ZstageScan', false, @islogical);
ip.addParameter('sCMOSCameraFlip', false, @islogical);
ip.addParameter('Save16bit', [false, false, false, false], @(x) (numel(x) == 1 || numel(x) == 4) && islogical(x));
ip.addParameter('onlyFirstTP', false, @islogical);
ip.addParameter('dzFromEncoder', false, @islogical);
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('saveZarr', false, @islogical); % use zarr file as output
ip.addParameter('save3DStack', true , @islogical); % option to save 3D stack or not
% pipeline steps
ip.addParameter('Deskew', true, @islogical);
ip.addParameter('Rotate', true, @islogical);
ip.addParameter('Stitch', false, @islogical);
ip.addParameter('Decon', ~false, @islogical);
ip.addParameter('RotateAfterDecon', false, @islogical);
% deskew and rotation options
ip.addParameter('parseSettingFile', false, @islogical); % use setting file to decide whether filp Z stack or not.
ip.addParameter('flipZstack', false, @islogical); % 
ip.addParameter('DSRCombined', true, @islogical); 
ip.addParameter('LLFFCorrection', false, @islogical);
ip.addParameter('BKRemoval', false, @islogical);
ip.addParameter('LowerLimit', 0.4, @isnumeric); % this value is the lowest
ip.addParameter('constOffset', [], @(x) isempty(x) || isnumeric(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('LSImagePaths', {'','',''}, @iscell);
ip.addParameter('BackgroundPaths', {'','',''}, @iscell);
ip.addParameter('resampleType', 'isotropic', @ischar); % resample type: given, isotropic, xy_isotropic
ip.addParameter('resample', [], @isnumeric); % resample
ip.addParameter('InputBbox', [], @isnumeric); % bbox for input in deskew and rotate
% stitch parameters
ip.addParameter('stitchPipeline', 'matlab', @ischar); % matlab or zarr
ip.addParameter('stitchResultDir', '', @ischar);
ip.addParameter('imageListFullpaths', '', @(x) ischar(x) || iscell(x));
ip.addParameter('axisOrder', 'xyz', @(x) ischar(x));
ip.addParameter('BlendMethod', 'none', @ischar);
ip.addParameter('xcorrShift', false, @islogical);
ip.addParameter('xcorrMode', 'primaryFirst', @(x) strcmpi(x, 'primary') || strcmpi(x, 'primaryFirst') || strcmpi(x, 'all')); % 'primary': choose one channel as primary channel, 
ip.addParameter('xyMaxOffset', 300, @isnumeric); % max offsets in xy axes
ip.addParameter('zMaxOffset', 50, @isnumeric); % max offsets in z axis
ip.addParameter('timepoints', [], @isnumeric); % stitch for given time points
ip.addParameter('boundboxCrop', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6));
ip.addParameter('primaryCh', '', @ischar);
ip.addParameter('stitchMIP', [], @(x) islogical(x) && (numel(x) == 1 || numel(x) == 3)); % 1x3 vector or vector, byt default, stitch MIP-z
ip.addParameter('onlineStitch', false, @(x) islogical(x)); % support for online stitch (with partial number of tiles). 
ip.addParameter('generateImageList', '', @(x) ischar(x)); % for real time processing, {'', 'from_encoder', 'from_sqlite'}
% decon parameters
ip.addParameter('cudaDecon', false, @islogical);
ip.addParameter('cppDecon', ~false, @islogical);
ip.addParameter('cppDeconPath', '/global/home/groups/software/sl-7.x86_64/modules/RLDecon_CPU/20200718/build-cluster/cpuDeconv', @ischar);
ip.addParameter('loadModules', 'module load gcc/4.8.5; module load fftw/3.3.6-gcc; module load boost/1.65.1-gcc; module load libtiff/4.1.0; ', @ischar);
ip.addParameter('cudaDeconPath', '/global/home/groups/software/sl-7.x86_64/modules/cudaDecon/bin/cudaDeconv' , @ischar);
ip.addParameter('OTFGENPath', '/global/home/groups/software/sl-7.x86_64/modules/cudaDecon/bin/radialft' , @ischar); % point to radialft file
ip.addParameter('DS', true, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('dzPSF', 0.1, @isnumeric);
ip.addParameter('EdgeErosion', 8, @isnumeric);
ip.addParameter('ErodeByFTP', true, @islogical); % Edge erosion by the first time point (ranked the first in the inital file list for each dataset).
ip.addParameter('deconRotate', false, @islogical);
ip.addParameter('psfFullpaths', {'','',''}, @iscell);
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('rotatedPSF', false , @islogical); % psf is rotated (for dsr)
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('fixIter', false, @islogical); % CPU Memory in Gb
ip.addParameter('errThresh', [], @isnumeric); % error threshold for simplified code
ip.addParameter('debug', false, @islogical); % debug mode for simplified code
ip.addParameter('GPUJob', false, @islogical); % use gpu for chuck deconvolution. 
% job related parameters
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 1, @isnumeric);
ip.addParameter('minModifyTime', 1, @isnumeric); % the minimum during of last modify time of a file, in minute.
ip.addParameter('maxModifyTime', 10, @isnumeric); % the maximum during of last modify time of a file, in minute.
ip.addParameter('maxWaitLoopNum', 10, @isnumeric); % the max number of loops the loop waits with all existing files processed. 
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2022a; matlab -nodisplay -nosplash -nodesktop -nojvm -r', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root of XR_Repository or LLSM5DTools. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../../'];
cd(repo_rt);
if ~exist([repo_rt, 'setup.m'], 'file')
    repo_rt = [mpath, '/../'];
    cd(repo_rt);
end

pr = ip.Results;
Overwrite = pr.Overwrite;
Streaming = pr.Streaming;
% Resolution = pr.Resolution;
SkewAngle = pr.SkewAngle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
ObjectiveScan = pr.ObjectiveScan;
ZstageScan = pr.ZstageScan;
Reverse = pr.Reverse;
ChannelPatterns = pr.ChannelPatterns;
Save16bit = pr.Save16bit;
resampleType = pr.resampleType;
resample = pr.resample;
dzFromEncoder = pr.dzFromEncoder;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
save3DStack = pr.save3DStack; % only for DS and DSR for now
%deskew and rotate
Deskew = pr.Deskew;
Rotate = pr.Rotate;
parseSettingFile = pr.parseSettingFile;
flipZstack = pr.flipZstack;
DSRCombined = pr.DSRCombined;
LLFFCorrection = pr.LLFFCorrection;
BKRemoval = pr.BKRemoval;
LowerLimit = pr.LowerLimit;
constOffset = pr.constOffset;
LSImagePaths = pr.LSImagePaths;
BackgroundPaths = pr.BackgroundPaths;
InputBbox = pr.InputBbox;
% stitch parameters
Stitch = pr.Stitch;
stitchPipeline = pr.stitchPipeline;
stitchResultDir = pr.stitchResultDir;
imageListFullpaths = pr.imageListFullpaths;
axisOrder = pr.axisOrder;
BlendMethod = pr.BlendMethod;
xcorrShift = pr.xcorrShift;
xcorrMode = pr.xcorrMode;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
boundboxCrop = pr.boundboxCrop;
onlyFirstTP = pr.onlyFirstTP;
timepoints = pr.timepoints;
primaryCh = pr.primaryCh;
stitchMIP = pr.stitchMIP;
onlineStitch = pr.onlineStitch;
generateImageList = pr.generateImageList;
% decon parameters
Decon = pr.Decon;
cppDecon = pr.cppDecon;
cudaDecon = pr.cudaDecon;
cppDeconPath = pr.cppDeconPath;
loadModules = pr.loadModules;
cudaDeconPath = pr.cudaDeconPath;
OTFGENPath = pr.OTFGENPath;
EdgeErosion = pr.EdgeErosion;
ErodeByFTP = pr.ErodeByFTP;
DS = pr.DS;
DSR = pr.DSR;
Background = pr.Background;
dzPSF = pr.dzPSF;
psfFullpaths = pr.psfFullpaths;
deconRotate = pr.deconRotate;
RotateAfterDecon = pr.RotateAfterDecon;
DeconIter = pr.DeconIter;
rotatedPSF = pr.rotatedPSF;
RLMethod = pr.RLMethod;
GPUJob = pr.GPUJob;
% matlab decon simplified version related options
fixIter = pr.fixIter;
errThresh = pr.errThresh;
debug = pr.debug;
% job related
largeFile = pr.largeFile;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes = pr.cpuOnlyNodes;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
minModifyTime = pr.minModifyTime;
maxModifyTime = pr.maxModifyTime;
maxWaitLoopNum = pr.maxWaitLoopNum;
MatlabLaunchStr = pr.MatlabLaunchStr;
SlurmParam = pr.SlurmParam;

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
if Streaming
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
    psfFullpaths = cellfun(@(x) strrep(x, '\', '/'), psfFullpaths, 'unif', 0);    
end

if numel(Overwrite) == 1
    Overwrite = repmat(Overwrite, 1, 5);
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str, jobLogDir] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
end

% for stitching, enable DS and DSR
% For multiple datasets, if the image list is not provide for any dataset,
% set Stitch as false. 
if Stitch
    if ~Streaming
        if numel(imageListFullpaths) ~= nd
            warning("The number of image list files does not match that of the data paths, please make sure image list files are provided for each dataset!")
            Stitch = false;
        else
            for d = 1 : nd
                if ~exist(imageListFullpaths{d}, 'file')
                    warning('Image list filename %s does not exist, set Stitch option as False', imageListFullpaths{d});
                    Stitch = false;
                    break;
                end
            end
        end
    end
    
    Deskew = true;
    Rotate = true;
    stchPaths = cell(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};
        if isempty(stitchResultDir)
            if strcmp(BlendMethod, 'none')
                stitchResultDir = 'matlab_stitch';
            else
                stitchResultDir = sprintf('matlab_stitch_%s', BlendMethod);                
            end
        end
        stchPath = [dataPath, '/', stitchResultDir, '/'];
        if Overwrite(3) && exist(stchPath, 'dir')
            rmdir(stchPath, 's');
        end
        if ~exist(stchPath, 'dir')
            mkdir(stchPath);
            fileattrib(stchPath, '+w', 'g');
        end
        stchPaths{d} = stchPath;
    end
    % DSRDirstr = '/DSR/';
    if ~strcmp(xcorrMode, 'primaryFirst') && Decon && ErodeByFTP
        ErodeByFTP = false;
    end
    
    % check if axis order is valid
    axisOrder = strrep(axisOrder, ' ', '');
    pattern = '^(-?x,?-?y,?-?z|-?y,?-?x,?-?z|-?z,?-?y,?-?x|-?x,?-?z,?-?y|-?x,?-?z,?-?y|-?y,?-?z,?-?x)$';
    if ~regexpi(axisOrder, pattern)
        error("The axisOrder is not right, it must has the form like 'y,x,z' or '-x,y,z' (flipped in x-axis)!");
    end
end

% first check if DS and Deconvolution directories exist
if ~Rotate
    DSRCombined = false;
end

if DSRCombined
    Deskew = true;
    Rotate = true;
end

if Deskew
    if ~DSRCombined
        dsPaths = cell(nd, 1);
        for d = 1 : nd
            dataPath = dataPaths{d};
            dsPath = [dataPath, '/DS/'];
            if Overwrite(1) && exist(dsPath, 'dir')
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
    if LLFFCorrection
        if numel(ChannelPatterns) ~= numel(LSImagePaths) || numel(ChannelPatterns) ~= numel(BackgroundPaths) 
            error('The number of channels in ChannelPatterns does not match the number of LS files or Background Files!')
        end
        for c = 1 : numel(LSImagePaths)
            if ~exist(LSImagePaths{c}, 'file')
                error('LS Image file %s does not exist!', LSImagePaths{c});
            end
            if ~exist(BackgroundPaths{c}, 'file')
                error('LS Image file %s does not exist!', BackgroundPaths{c});
            end
        end
    end
else
    Rotate = false;
    DS = false;
    DSR = false;
end
    
if Rotate
    dsrPaths = cell(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};
        dsrPath = [dataPath, '/DSR/'];
        if Overwrite(2) && exist(dsrPath, 'dir')
            rmdir(dsrPath, 's');
        end
        if ~exist(dsrPath, 'dir')
            mkdir(dsrPath);
            fileattrib(dsrPath, '+w', 'g');
        end
        dsrPaths{d} = dsrPath;
    end
end

% define resample based on resample type, resample parameter and the microscope setting
if ~isempty(resample)
    resampleType = 'given';
end
[resample, zAniso] = XR_checkResampleSetting(resampleType, resample, ObjectiveScan, SkewAngle, xyPixelSize, dz);


% get actual dz from encoder positions
if dzFromEncoder && ~Streaming
    dz_all = zeros(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};        
        dz_actual = XR_estimate_actual_step_size_from_encoder(dataPath, 'dz', dz);
        dz_all(d) = dz_actual;
    end   
else
    dz_all = repmat(dz, nd, 1);
end


% for deconvolution, check whether there is a gpu in the node. if not, for
% cudaDecon, set parseCluster as true. 
deconPaths = cell(nd, 1);
if Decon
    % if both cudaDecon and cppDecon are true, use cppDecon
    if cudaDecon && cppDecon
        cudaDecon = false;
    end

    if cudaDecon && gpuDeviceCount() < 1 && ~parseCluster
        warning('There is no GPU in the node, and ther cluster is also not available. Set cudaDecon as false!');
        cudaDecon = false;
    end
        
    if cudaDecon
        deconName = 'GPUdecon';
    elseif cppDecon
        deconName = 'CPPdecon';
    else
        deconName = 'matlab_decon';
    end
        
    if RotateAfterDecon
        rdcPaths = cell(nd, 1);
    end

    for d = 1 : nd
        dataPath = dataPaths{d};

        deconPath = [dataPath, deconName, '/'];
        if DS
            if DSRCombined
                error('If using DS for deconvolution, "DSRCombined" must be set as false!')
            end
            dsPath = dsPaths{d};
            deconPath = [dsPath, deconName, '/'];
        end
        if DSR
            dsrPath = dsrPaths{d};
            deconPath = [dsrPath, deconName, '/'];
        end
        if Stitch
            stchPath = stchPaths{d};            
            deconPath = [stchPath, deconName, '/'];
            RotateAfterDecon = false;
        end

        if Overwrite(4) && exist(deconPath, 'dir')
            rmdir(deconPath, 's');
        end
        if ~exist(deconPath, 'dir')
            mkdir(deconPath);
            fileattrib(deconPath, '+w', 'g');
        end
        deconPaths{d} = deconPath;
        
        if RotateAfterDecon
            rdcPath = [deconPath '/' 'Rotated' '/'];
            if Overwrite(5) && exist(rdcPath, 'dir')
                rmdir(rdcPath, 's');
            end
            if ~exist(rdcPath, 'dir')
                mkdir(rdcPath);
                fileattrib(rdcPath, '+w', 'g');
            end
            rdcPaths{d} = rdcPath;

            if DSR 
                warning('The rotation is already performed before deconvolution! Please check the setting to make sure there is no duplicate rotation.');
            end
        end
    end
    
    % check whether a psf file exist, if not, show a warning
    for f = 1 : numel(psfFullpaths)
        if ~exist(psfFullpaths{f}, 'file')
            error('PSF file %s does not exist!', psfFullpaths{f});
        end
    end
    if rotatedPSF 
        if DS
            error('The PSF must be unrotated!');
        end
        rotPSFFullpaths = psfFullpaths;
    else
        if (DSR || Stitch)
            rotPSFFullpaths = cell(numel(psfFullpaths), 1);
        end

        for f = 1 : numel(psfFullpaths)
            [psfPath, fsname] = fileparts(psfFullpaths{f});        
            rotPSFFullpaths{f} = [psfPath, '/Rotated/', fsname, '.tif'];
            if (DSR || Stitch) && ~exist(rotPSFFullpaths{f}, 'file')
                XR_rotate_PSF(psfFullpaths{f}, 'Reverse', Reverse);
                % XR_rotate_PSF(psfFullpaths{f});
            end
        end
    end
else
    RotateAfterDecon = false;
end

%% check existing files and parse channels
[fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, latest_modify_times, FTP_inds, maskFullpaths] = ...
    XR_parseImageFilenames(dataPaths, ChannelPatterns, parseSettingFile, flipZstack, Decon, deconPaths, Streaming, minModifyTime, zarrFile);

nF = numel(fnames);

% flags: for thee: deskew w/o rotate, decon w/o rotate, rotate
is_done_flag = false(nF, 4);
% ensure only one stitch wrapper is running of a dataset
stitch_running = false(nd, 1); 
trial_counter = zeros(nF, 4);
waitLoopCounter = 0;

% For no-streaming computing, first skip the disabled options. 
if ~Streaming
    if ~Stitch
        is_done_flag(:, 2) = true;
    end
    if ~Decon
        is_done_flag(:, 3) = true;
    end
    if ~RotateAfterDecon
        is_done_flag(:, 4) = true;
    end
end

if parseCluster
    job_ids = -ones(nF, 4);
    imSize_mat = zeros(nF, 3, 2);
    dataSize_mat = zeros(nF, 2);
    dataSize_mat(:, 1) = dataSizes;
end

matlab_setup_str = 'setup([],true)';

% use while loop to perform computing for all images
ts = tic;
while ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') || ...
        (Streaming && (nF == 0 || any(latest_modify_times < maxModifyTime) || waitLoopCounter < maxWaitLoopNum))
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
        if Deskew
            if ~DSRCombined
                dsPath = dsPaths{fdind};
                dsFullpath = [dsPath, fsname, '.tif'];
                tmpFullpath = sprintf('%s.tmp', dsFullpath(1 : end - 4));
            end
            if Rotate
                dsrPath = dsrPaths{fdind};
                dsrFullpath = [dsrPath, fsname, '.tif'];
                tmpFullpath = sprintf('%s.tmp', dsrFullpath(1 : end - 4));                
            end

            if (DSRCombined || exist(dsFullpath, 'file')) && (~Rotate || exist(dsrFullpath, 'file'))
                is_done_flag(f, 1) = true;
                if exist(tmpFullpath, 'file')
                    delete(tmpFullpath);
                end
            end
        else
            is_done_flag(f, 1) = true;
        end
        
        if ~is_done_flag(f, 1) 
            if LLFFCorrection
                % LLFFMapping =  ~cellfun(@isempty, regexpi(fname, ChannelPatterns));
                % change to contains.m to unify the matching
                LLFFMapping =  cellfun(@(x) contains(frameFullpath, x), ChannelPatterns);
                LSImage = LSImagePaths{LLFFMapping};
                BackgroundImage = BackgroundPaths{LLFFMapping};
            else
                LSImage = '';
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
            
            func_str = sprintf(['XR_deskewRotateFrame(%s,%.20d,%.20d,''SkewAngle'',%.20d,', ...
                '''ObjectiveScan'',%s,''ZstageScan'',%s,''Reverse'',%s,''LLFFCorrection'',%s,', ...
                '''BKRemoval'',%s,''LowerLimit'',%.20d,''constOffset'',[%s],''LSImage'',''%s'',', ...
                '''BackgroundImage'',''%s'',''Rotate'',%s,''resample'',[%s],''DSRCombined'',%s,', ...
                '''InputBbox'',%s,''flipZstack'',%s,''Save16bit'',%s,''save3DStack'',%s)'], ...
                ds_input_str, xyPixelSize, dz_f, SkewAngle, string(ObjectiveScan), string(ZstageScan), ...
                string(Reverse),  string(LLFFCorrection), string(BKRemoval), LowerLimit, ...
                num2str(constOffset, '%0.10f'),  LSImage, BackgroundImage, string(Rotate), ...
                strrep(num2str(resample, '%d,'), ' ', ''),  string(DSRCombined), strrep(mat2str(InputBbox), ' ', ','), ...
                string(flipZstack), string(Save16bit(1)), string(save3DStack));
            
            if exist(tmpFullpath, 'file') || parseCluster
                if parseCluster
                    job_status = check_slurm_job_status(job_ids(f, 1), task_id);
                    
                    % if the job is still running, skip it. 
                    if job_status == 1 
                        continue;
                    end
                    
                    if job_status == -1
                        % first estimate file size and decide whether cpusPerTask
                        % is enough
                        if ~DSRCombined
                            estRequiredMemory = XR_estimateComputingMemory(frameFullpath, 'steps', {'deskew'});
                        else
                            estRequiredMemory = dataSize_mat(f, 1) / 2^30 * 2 * (7 + 3 / prod(resample) + trial_counter(f, 1) * 8);
                        end
                        cpusPerTask_ds = cpusPerTask;
                        if cpusPerTask_ds * 20 < estRequiredMemory
                            cpusPerTask_ds = min(24, ceil(estRequiredMemory / 20));
                        else
                            cpusPerTask_ds = min(cpusPerTask_ds, ceil(estRequiredMemory / 20));
                        end
                            
                        matlab_cmd = sprintf('%s;tic;%s;toc', matlab_setup_str, func_str);
                        process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                            '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], task_id, job_log_fname, ...
                            job_log_error_fname, cpusPerTask_ds, SlurmParam, slurm_constraint_str, ...
                            matlab_cmd, process_cmd);
                        [status, cmdout] = system(cmd, '-echo');

                        job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                        job_id = str2double(job_id{1}{1});
                        job_ids(f, 1) = job_id;
                        trial_counter(f, 1) = trial_counter(f, 1) + 1;
                    end
                else
                    temp_file_info = dir(tmpFullpath);
                    if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 < unitWaitTime
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
                trial_counter(f, 1) = trial_counter(f, 1) + 1;
            end

            % check if computing is done
            if (DSRCombined || exist(dsFullpath, 'file')) && (~Rotate || exist(dsrFullpath, 'file'))
                is_done_flag(f, 1) = true;
                if exist(tmpFullpath, 'file')
                    delete(tmpFullpath);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% then stitching
        if Stitch
            if is_done_flag(f, 2)
                continue;
            end
            stchPath = stchPaths{fdind};
            % dir_info = dir(sprintf('%s/%s*Abs.tif', stchPath, fname(regexp(fname, 'Scan.*') : end - 43)));
            switch stitchPipeline
                case 'zarr'
                    stch_dir_info = dir(sprintf('%s/%s*Abs.zarr', stchPath, fname(1 : end - 43)));                                                   
                case 'tiff'
                    stch_dir_info = dir(sprintf('%s/%s*Abs.tif', stchPath, fname(1 : end - 43)));                                
            end
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
            if Streaming && ~isempty(generateImageList) 
                stitch_flags_d = is_done_flag(:, 2) == 0 & fdinds == fdind;
                if any(stitch_flags_d) && (f == find(stitch_flags_d, 1, 'first'))
                    switch generateImageList
                        case 'from_encoder'
                            imageListFullpaths{fdind} = stitch_generate_imagelist_from_encoder(dataPath, dz_f, ChannelPatterns);
                        case 'from_sqlite'
                            imageListFullpaths{fdind} = stitch_generate_imagelist_from_sqlite(dataPath);                        
                    end
                end
            end
        else
            is_done_flag(f, 2) = true;
        end
        
        if ~is_done_flag(f, 2)
            % if the deskew result is not available, wait for the deskew data
            % to be finished        
            if ~parseCluster && (~exist(frameFullpath, 'file') || ~exist(dsrFullpath, 'file'))
                % if DS or DSR, it means the deskew and rotation not done.
                continue;
            end
            
            stitch_DS = false;
            stitch_DSR = true;
            useProcessedData = true;
            ProcessedDirStr = 'DSR';
            resampleType = 'isotropic';
            zNormalize = false;
            % onlyFirstTP = false;
            bbox = boundboxCrop;
            imageListFullpath = imageListFullpaths{fdind};
            ChannelPatterns_str = sprintf('{''%s''}', strjoin(ChannelPatterns, ''','''));

            func_str = sprintf(['XR_matlab_stitching_wrapper(''%s'',''%s'',''ResultDir'',''%s'',', ...
                '''Streaming'',%s,''DS'',%s,''DSR'',%s,''ChannelPatterns'',%s,''useProcessedData'',%s,', ...
                '''axisOrder'',''%s'',''resampleType'',''%s'',''resample'',[%s],''Reverse'',%s,', ...
                '''parseSettingFile'',%s,''xcorrShift'',%s,''xcorrMode'',''%s'',''xyMaxOffset'',%.10f,', ...
                '''zMaxOffset'',%.10f,''BlendMethod'',''%s'',''zNormalize'',%s,''onlyFirstTP'',%s,', ...
                '''timepoints'',[%s],''boundboxCrop'',[%s],''Save16bit'',%s,''primaryCh'',''%s'',', ...
                '''stitchMIP'',%s,''onlineStitch'',%s,''pipeline'',''%s'')'], dataPath, imageListFullpath, ...
                stitchResultDir, string(Streaming), string(stitch_DS), string(stitch_DSR), ChannelPatterns_str, ...
                string(useProcessedData), ProcessedDirStr, axisOrder, resampleType, strrep(num2str(resample, '%.10d,'), ' ', ''), ...
                string(Reverse), string(parseSettingFile), string(xcorrShift), xcorrMode, xyMaxOffset, ...
                zMaxOffset, BlendMethod, string(zNormalize), string(onlyFirstTP), strrep(num2str(timepoints, '%d,'), ' ', ''), ...
                strrep(num2str(bbox, '%d,'), ' ', ''), string(Save16bit(2)), primaryCh, ...
                strrep(mat2str(stitchMIP), ' ', ','), string(onlineStitch), stitchPipeline);

            if parseCluster
                dfirst_ind = find(fdinds == fdind, 1, 'first');
                job_status = check_slurm_job_status(job_ids(dfirst_ind, 2));

                if job_status == -1 % && ~stitch_running(fdind)
                    % first estimate file size and decide whether cpusPerTask
                    % is enough
                    memFactor = 20;
                    if any(stitchMIP)
                        memFactor = 2;
                    end
                    estRequiredMemory = dataSize_mat(f, 1) / 2^30 * 2 * memFactor;
                    cpusPerTask_stch = cpusPerTask;
                    if cpusPerTask_stch * 20 < estRequiredMemory
                        cpusPerTask_stch = min(24, ceil(estRequiredMemory / 20));
                    end

                    matlab_cmd = sprintf('%s;tic;%s;toc', matlab_setup_str, func_str);
                    process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                    cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                        '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], ...
                        task_id, job_log_fname, job_log_error_fname, cpusPerTask_stch, SlurmParam, ...
                        slurm_constraint_str, matlab_cmd, process_cmd);
                    [status, cmdout] = system(cmd, '-echo');

                    job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                    job_id = str2double(job_id{1}{1});
                    job_ids(dfirst_ind, 2) = job_id;
                    % stitch_running(fdind) = true;
                end
            else
                tic; feval(str2func(['@()', func_str])); toc;
            end
            
            if ~isempty(stch_dir_info) && exist(stchFullpath, 'file')
                is_done_flag(f, 2) = true;
                if exist(stchTmpFullpath, 'file')
                    delete(stchTmpFullpath);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% then deconvolution
        if Decon
            % if the deskew result is not available, wait for the deskew data
            % to be finished

            % input chosen order is dsr, ds, raw (in decreased priority order)
            dcframeFullpath = frameFullpath;
            
            if DS
                dcframeFullpath = dsFullpath;
                dc_dz = dz_f;
                dc_dzPSF = dzPSF;
                dc_psfFullpaths = psfFullpaths;
            end
            
            if DSR
                dcframeFullpath = dsrFullpath;
                dc_dz = xyPixelSize;
                if rotatedPSF
                    dc_dzPSF = dzPSF;
                else
                    dc_dzPSF = xyPixelSize;
                end
                dc_psfFullpaths = rotPSFFullpaths;
            end 

            if Stitch
                % in case fname not start with Scan_
                % dir_info = dir(sprintf('%s/%s*Abs.tif', stchPath, fname(regexp(fname, 'Scan.*') : end - 43)));
                if strcmp(stitchPipeline, 'zarr')
                    dir_info = dir(sprintf('%s/%s*Abs.zarr', stchPath, fname(1 : end - 43)));
                else
                    dir_info = dir(sprintf('%s/%s*Abs.tif', stchPath, fname(1 : end - 43)));
                end
                if isempty(dir_info) 
                    continue;
                end
                stch_fname = dir_info(1).name;
                [~, stch_fsname] = fileparts(stch_fname);
                dcframeFullpath = [stchPath, stch_fname];
                dc_dz = xyPixelSize;
                dc_dzPSF = xyPixelSize;
                dc_psfFullpaths = rotPSFFullpaths;
                
                % if stitch, only do computing when it is for the first
                % tile, to avoid replicate computing
                if ~contains(fname, stch_fsname)
                    is_done_flag(f, 3) = true;
                    continue;
                end
            end
            
            if ~exist(dcframeFullpath, 'file')
                % if DS or DSR, it means the deskew and rotation not done.
                if DS || DSR || Stitch
                    continue;
                end
            end
                            
            deconPath = deconPaths{fdind};
            deconFullpath = sprintf('%s/%s_decon.tif', deconPath, fsname);
            if Stitch
                deconFullpath = sprintf('%s/%s_decon.tif', deconPath, stch_fsname);
            end
            dctmpFullpath = sprintf('%s.tmp', deconFullpath(1 : end - 4));

            if exist(deconFullpath, 'file')
                is_done_flag(f, 3) = true;
                if exist(dctmpFullpath, 'file')
                    delete(dctmpFullpath);
                end
            end
            
            if parseCluster
                if dataSize_mat(f, 2) == 0 && (~all(is_done_flag(f, 3 : 4)) || f == FTP_ind)
                    dir_info = dir(dcframeFullpath);
                    dataSize_mat(f, 2) = dir_info.bytes;
                end
            end
            
            % for ErodeByFTP, check if the mask file exist
            maskFullpath = '';
            SaveMaskfile = false;
            if ErodeByFTP && ~cudaDecon
                if f == FTP_ind
                    SaveMaskfile = true;
                    % if decon result exist, but mask file not exist, rerun
                    % it to save the mask. 
                    if Stitch
                        maskFullpaths{fdind} = sprintf('%s/Masks/%s_eroded.tif', deconPath, stch_fsname);
                    end
                    if is_done_flag(f, 3) && ~exist(maskFullpaths{fdind}, 'file')
                        is_done_flag(f, 3) = false;
                        fprintf('Mask file %s does not exist, delete the deconvolved result\n', maskFullpaths{fdind});
                        delete(deconFullpath);
                    end
                else
                    % only check for the ones not finished. 
                    if ~is_done_flag(f, 3)
                        maskFullpath = maskFullpaths{fdind};
                        if ~exist(maskFullpath, 'file')
                            continue;
                        end

                        mask_sz = getImageSize(maskFullpath);
                        img_sz = getImageSize(dcframeFullpath);
                        if any(mask_sz ~= img_sz)
                            warning(['The image size [%s] does not match the defined ', ... 
                                'mask size [%s], use its own mask for edge erosion...'], ...
                                num2str(img_sz, '%d '), num2str(mask_sz, '%d '));
                            maskFullpath = '';
                        end
                    end
                end
            end            
        else
            is_done_flag(f, 3) = true;    
        end

        if ~is_done_flag(f, 3) 
            % psfMapping =  ~cellfun(@isempty, regexpi(fname, ChannelPatterns));
            % change to contains.m to unify the matching
            psfMapping =  ~cellfun(@isempty, regexpi(frameFullpath, ChannelPatterns));
            
            psfFullpath = dc_psfFullpaths{psfMapping};
            
            % do not use rotation in decon functions
            if cudaDecon
                func_str = sprintf(['XR_cudaDeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''cudaDeconPath'',''%s'',''OTFGENPath'',''%s'',''dzPSF'',%.10f,''Background'',[%d],', ...
                    '''SkewAngle'',%d,''Rotate'',%s,''DeconIter'',%d,''Save16bit'',%s,''largeFile'',%s)'], ...
                    dcframeFullpath, xyPixelSize, dc_dz, psfFullpath, cudaDeconPath, OTFGENPath, dc_dzPSF, ...
                    Background, SkewAngle, string(deconRotate), DeconIter, string(Save16bit(3)), string(largeFile));
            elseif cppDecon
                func_str = sprintf(['XR_cppDeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''cppDeconPath'',''%s'',''loadModules'',''%s'',''dzPSF'',%.10f,''Background'',[%d],', ...
                    '''SkewAngle'',%d,''EdgeErosion'',%d,''ErodeMaskfile'',''%s'',''SaveMaskfile'',%s,', ...
                    '''Rotate'',%s,''DeconIter'',%d,''Save16bit'',%s,''largeFile'',%s)'], dcframeFullpath, ...
                    xyPixelSize, dc_dz, psfFullpath, cppDeconPath, loadModules, dc_dzPSF, Background, ...
                    SkewAngle, EdgeErosion, maskFullpath, string(SaveMaskfile), string(deconRotate), ...
                    DeconIter, string(Save16bit(3)), string(largeFile));
            else
                func_str = sprintf(['XR_RLdeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''dzPSF'',%.10f,''Background'',[%d],''SkewAngle'',%d,''EdgeErosion'',%d,''ErodeMaskfile'',''%s'',', ...
                    '''SaveMaskfile'',%s,''Rotate'',%s,''DeconIter'',%d,''RLMethod'',''%s'',''fixIter'',%s,', ...
                    '''errThresh'',[%0.20f],''debug'',%s,''GPUJob'',%s,''Save16bit'',%s,''largeFile'',%s)'], ...
                    dcframeFullpath, xyPixelSize, dc_dz, psfFullpath,  dc_dzPSF, Background, SkewAngle, ...
                    EdgeErosion, maskFullpath, string(SaveMaskfile), string(deconRotate), DeconIter, RLMethod, ...
                    string(fixIter), errThresh, string(debug), string(GPUJob), string(Save16bit(3)), string(largeFile));
            end
           
            if exist(dctmpFullpath, 'file') || parseCluster
                if parseCluster
                    job_status = check_slurm_job_status(job_ids(f, 3), task_id);

                     % if the job is still running, skip it. 
                    if job_status == 1 
                        continue;
                    end

                    if job_status == -1
                        % for matlab decon,  decide how many cores. 
                        if ~cudaDecon
                            [estMem, estGPUMem] = XR_estimateComputingMemory('', {'deconvolution'}, ...
                                'dataSize', dataSize_mat(f, 2), 'cudaDecon', false);
                            cpusPerTask_dc = cpusPerTask;
                            if cpusPerTask_dc * 20 < estMem
                                cpusPerTask_dc = min(24, ceil(estMem / 20));
                            end
                        end

                        % do not use rotation in decon functions
                        matlab_cmd = sprintf('%s;tic;%s;toc', matlab_setup_str, func_str);
                        process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);

                        if cudaDecon
                            cmd = sprintf(['sbatch --array=%d -o %s -e %s -p abc --gres=gpu:1 --qos ', ...
                                'abc_normal -n1 --mem-per-cpu=33G --cpus-per-task=%d --wrap="%s"'], ...
                                task_id, job_log_fname, job_log_error_fname, 5, process_cmd);
                        else
                            cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                                '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], task_id, job_log_fname, ...
                                job_log_error_fname, cpusPerTask_dc, SlurmParam, slurm_constraint_str, ...
                                matlab_cmd, process_cmd);
                        end
                        [status, cmdout] = system(cmd, '-echo');

                        job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                        job_id = str2double(job_id{1}{1});
                        job_ids(f, 3) = job_id;
                        trial_counter(f, 3) = trial_counter(f, 3) + 1;
                    end
                else
                    temp_file_info = dir(dctmpFullpath);
                    if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 < unitWaitTime
                        continue; 
                    else
                        fclose(fopen(dctmpFullpath, 'w'));
                    end
                end
            else
                fclose(fopen(dctmpFullpath, 'w'));
            end

            if ~parseCluster
                % fileInfo = imfinfo(dsFullpath);
                tic; feval(str2func(['@()', func_str])); toc;
                trial_counter(f, 3) = trial_counter(f, 3) + 1;
            end
            
            % check if computing is done
            if exist(deconFullpath, 'file')
                is_done_flag(f, 3) = true;
                if exist(dctmpFullpath, 'file')
                    delete(dctmpFullpath);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% rotation after decon
        if RotateAfterDecon
            % input chosen order is dsr, ds, raw (in decreased priority order)
            if ~exist(deconFullpath, 'file')
                if Decon
                    continue;
                end
            end

            rdcPath = rdcPaths{fdind};
            rdcFullpath = sprintf('%s/%s_decon.tif', rdcPath, fsname);
            rdctmpFullpath = sprintf('%s.tmp', rdcFullpath(1 : end - 4));
            if exist(rdcFullpath, 'file')
                is_done_flag(f, 4) = true;
                if exist(rdctmpFullpath, 'file')
                    delete(rdctmpFullpath);
                end
            end
        else
            is_done_flag(f, 4) = true;    
        end

        if ~is_done_flag(f, 4)
            func_str = sprintf(['XR_RotateFrame3D(''%s'',%.20d,%.20d,''SkewAngle'',%.20d,', ...
                '''ObjectiveScan'',%s,''Reverse'',%s,''Save16bit'',%s)'], deconFullpath, ...
                xyPixelSize, dz_f, SkewAngle, string(ObjectiveScan), string(Reverse), string(Save16bit(4)));

            if exist(rdctmpFullpath, 'file') || parseCluster
                if parseCluster
                    job_status = check_slurm_job_status(job_ids(f, 4), task_id);
                    
                    % if the job is still running, skip it. 
                    if job_status == 1 
                        continue;
                    end
                    
                    if job_status == -1
                        [estMem, estGPUMem] = XR_estimateComputingMemory('', {'deconvolution'}, ...
                            'dataSize', dataSize_mat(f, 2), 'cudaDecon', false);
                        cpusPerTask_dcr = cpusPerTask;
                        if cpusPerTask_dcr * 20 < estMem
                            cpusPerTask_dcr = min(24, ceil(estMem / 20));
                        end

                        matlab_cmd = sprintf('%s;tic;%s;toc', matlab_setup_str, func_str);
                        process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                            '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask_dcr, SlurmParam, ...
                            slurm_constraint_str, matlab_cmd, process_cmd);

                        [status, cmdout] = system(cmd, '-echo');

                        job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                        job_id = str2double(job_id{1}{1});
                        job_ids(f, 4) = job_id;
                        trial_counter(f, 4) = trial_counter(f, 4) + 1;    
                    end
                else
                    temp_file_info = dir(rdctmpFullpath);
                    if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 < unitWaitTime
                    else
                        fclose(fopen(rdctmpFullpath, 'w'));
                    end
                end
            else
                fclose(fopen(rdctmpFullpath, 'w'));
            end
            
            if ~parseCluster
                tic; feval(str2func(['@()', func_str])); toc;
                trial_counter(f, 4) = trial_counter(f, 4) + 1;    
            end
            
            % check if computing is done
            if exist(rdcFullpath, 'file')
                is_done_flag(f, 4) = true;
                if exist(rdctmpFullpath, 'file')
                    delete(rdctmpFullpath);
                end
            end
        end
        toc
    end
    
    %% wait for running jobs finishing and checking for new coming images
    nF_done = sum(all(is_done_flag, 2));
    sprintf('Time %d s: %d / %d (%0.3f) are finished!\n', toc(ts), nF_done, nF, nF_done / nF);
    
    if ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
        waitLoopCounter = 0;
        pause(30);
        % if not streaming, do not check for new files. 
        if ~Streaming
            continue;
        end
    else
         % if not streaming, exit when all existing files are processed. 
        if Streaming 
            if waitLoopCounter < maxWaitLoopNum
                waitLoopCounter = waitLoopCounter + 1;
                pause(30);
            end
        else
            break;
        end
    end
            
    % check whether there are new coming images (only for streaming option)
    [fnames_new, fdinds_new, gfnames_new, partialvols_new, dataSizes_new, flipZstack_mat_new, latest_modify_times_new, FTP_inds_new, maskFullpaths_new] = ...
    XR_parseImageFilenames(dataPaths, ChannelPatterns, parseSettingFile, flipZstack, Decon, deconPaths, Streaming, minModifyTime, zarrFile);
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
    maskFullpaths = maskFullpaths_new;
    
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

