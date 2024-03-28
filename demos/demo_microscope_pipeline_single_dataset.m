% demo to run microscope automatic pipeline for a single dataset with stitch

% note: this demo is designed for quick processing of raw microscopy data 
% during acquistion, providing fast feedback and inspection. For final
% processing after acquistion, follow the functions outlined in the stitching, 
% deconvolution, geometric transformation, and large-scale processing demos.


clear, clc;

fprintf('Microscopy pipeline for a single dataset...\n\n');

% move to the LLSM5DTools root directory
curPath = pwd;
if ~endsWith(curPath, 'LLSM5DTools')
    mfilePath = mfilename('fullpath');
    if contains(mfilePath,'LiveEditorEvaluationHelper')
        mfilePath = matlab.desktop.editor.getActiveFilename;
    end
    
    mPath = fileparts(mfilePath);
    if endsWith(mPath, 'demos')
        cd(mPath);
        cd('..')
    end
end

setup();


%% Step 1: get our demo data from zenodo/Dropbox (skip this step if the data is already downloaded)
% download the example dataset from zenodo (https://doi.org/10.5281/zenodo.10471978) manually, 
% or use the code below to download the data from Dropbox
if ispc
    destPath = fullfile(getenv('USERPROFILE'), 'Downloads');
    destPath = strrep(destPath, '\', '/');
else
    destPath = '~/Downloads/';
end
demo_data_downloader(destPath);

dataPath = [destPath, '/LLSM5DTools_demo_cell_image_dataset/'];


%% Step 2: set parameters for each step. 
% add the software to the path
cd(fileparts(which(mfilename)));
addpath(genpath('../'));

% data path
% defined above

% general parameters.
% pixel size in xy-plane
xyPixelSize = 0.108;
% pixel size in z-axis
dz = 0.3;
% Inverse direction of z axis. 
Reverse = true;
% Channel identifiers and orders. 
ChannelPatterns = {'CamB_ch0', 'CamA_ch0'};
% Save 16bit for deskew/rotate, stitch, deconvolution, and rotation after decon
Save16bit = [true, true, true, true];
% not overwrite existing results
Overwrite = false;
% Processing for existing dataset
Streaming = false;
% use cluster for processing if true
parseCluster = true;
% number of cores for the pipeline (only for ABC cluster). For other
% cluster environment, the number of core is determined by the larger one of this
% setting and the estimated one based on the configFile. This setting is
% ignored for local environment. 
cpusPerTask = 24;
% configuration file for job submission
configFile = '';
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;

% define which step is included in the analysis
Deskew = true; % deskew
Rotate = true; % rotate after deskew
Stitch = true; % stitch
Decon = false; % deconvolution
RotateAfterDecon = false; % rotate after deconvolution

DSRCombined = true; % bypassing deskew and use combined processing

% flat field parameters (if LLFFCorrection is set to true) 
LLFFCorrection = true;
% lower bound cap for flat field
LowerLimit = 0.4;
% add a constant background value after correction to avoid zero background. 
% if empty, add the background image back after correction.
constOffset = 1;
% flat field image paths in the order of CamB_ch0 and CamA_ch0 as defined in ChannelPatterns
LSImagePaths = {[dataPath, '/FF/averaged/ff_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_000x_000y_000z_0017t.tif'],...
                [dataPath, '/FF/averaged/ff_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_000x_000y_000z_0017t.tif']};
% camera background image paths in the order of CamB and CamA as defined in ChannelPatterns
BackgroundPaths = {[dataPath, '/FF/KorraFusions/AVG_DF400_CamB_10ms.tif'],...
                   [dataPath, '/FF/KorraFusions/AVG_DF400_CamA_10ms.tif']};

% stitch parameters
% image list csv file path
ImageListFullpaths = [dataPath, 'ImageList_from_encoder.csv'];
% axis order for coordinates
axisOrder = '-x,y,z';

% parameters for real-time processing and feedback (when Streaming is true)
% The minimum time in minutes for the latest modified file to decide whether it is fully transferred.
minModifyTime = 1;
% The maximum time in minutes to check whether there are coming new files.
maxModifyTime = 10;
% Number of maximum loops without any computing.
maxWaitLoopNum = 10;


%% Step 3: run the analysis with given parameters. 

% result folders:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/DS/    (only when DSRCombined is false)
% {destPath}/LLSM5DTools_demo_cell_image_dataset/DSR/
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_stitch/

XR_microscopeAutomaticProcessing(dataPath, 'xyPixelSize', xyPixelSize, 'dz', dz,  ...
    'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'Save16bit', Save16bit, ...
    'Overwrite', Overwrite, 'Streaming', Streaming, 'Deskew', Deskew, 'Rotate', Rotate, ...
    'Stitch', Stitch, 'Decon', Decon, 'RotateAfterDecon', RotateAfterDecon, ...
    'ImageListFullpaths', ImageListFullpaths, 'axisOrder', axisOrder, 'LLFFCorrection', LLFFCorrection, ...
    'LowerLimit', LowerLimit, 'constOffset', constOffset, 'LSImagePaths', LSImagePaths, ...
    'BackgroundPaths', BackgroundPaths,'minModifyTime', minModifyTime, 'maxModifyTime', maxModifyTime, ...
    'maxWaitLoopNum', maxWaitLoopNum, 'parseCluster', parseCluster, 'cpusPerTask', cpusPerTask, ...
    'configFile', configFile, 'mccMode', mccMode);


