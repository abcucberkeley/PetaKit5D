% demo to run microscope automatic pipeline for multiple datasets with same 
% experiment settings with stitch

% note: this demo is designed for quick processing of raw microscopy data 
% during acquistion, providing fast feedback and inspection. For final
% processing after acquistion, follow the functions outlined in the stitching, 
% deconvolution, geometric transformation, and large-scale processing demos.


%% Step 1: set parameters for each step. 
% add the software to the path
cd(fileparts(which(mfilename)));
addpath(genpath('../'));

%  rt is the root directory of datasets (if all datasets are within this folder)
dataPaths = {'/clusterfs/fiona/Data/20200806_p35p4_Hex_Raptv_234-1_DLS/Exp01/', ...
            '/clusterfs/fiona/Data/20200806_p35p4_Hex_Raptv_234-1_DLS/Exp02/', ...
            '/clusterfs/fiona/Data/20200806_p35p4_Hex_Raptv_234-1_DLS/Exp03/'};            

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
LLFFCorrection = false;
% flat field image paths in the order of CamB_ch0 and CamA_ch0 as defined in ChannelPatterns
LSImagePaths = {'/clusterfs/fiona/Data/20200803_LLCPKI/20200803_Calibration/FF/AVG_488nm.tif',...
                '/clusterfs/fiona/Data/20200803_LLCPKI/20200803_Calibration/FF/AVG_560nm.tif'};

% camera background image paths in the order of CamB and CamA as defined in ChannelPatterns
BackgroundPaths = {'/clusterfs/fiona/OrcaDC/Aang/Orca_Aang_AVG_Cam_B.tif',...
                   '/clusterfs/fiona/OrcaDC/Aang/Orca_Aang_AVG_Cam_A.tif'};

% stitch parameters
% image list csv file path for each dataset
ImageListFullpaths = {[dataPaths{1}, 'ImageList_Exp01.csv'], ...
                      [dataPaths{2}, 'ImageList_Exp02.csv'], ...
                      [dataPaths{3}, 'ImageList_Exp03.csv']};
% axis order for coordinates
axisOrder = '-x,y,z';

% parameters for real-time processing and feedback (when Streaming is true)
% The minimum time in minutes for the latest modified file to decide whether it is fully transferred.
minModifyTime = 1;
% The maximum time in minutes to check whether there are coming new files.
maxModifyTime = 10;
% Number of maximum loops without any computing.
maxWaitLoopNum = 10;


%% Step 2: run the analysis with given parameters. 

% result folders:
% {dataPaths}/DS/    (only when DSRCombined is false)
% {dataPaths}/DSR/   
% {dataPaths}/matlab_stitch/    

XR_microscopeAutomaticProcessing(dataPaths, 'xyPixelSize', xyPixelSize, 'dz', dz,  ...
    'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'Save16bit', Save16bit, ...
    'Overwrite', Overwrite, 'Streaming', Streaming, 'Deskew', Deskew, 'Rotate', Rotate, ...
    'Stitch', Stitch, 'Decon', Decon, 'RotateAfterDecon', RotateAfterDecon, ...
    'ImageListFullpaths', ImageListFullpaths, 'axisOrder', axisOrder, 'DSRCombined', DSRCombined, ...
    'LLFFCorrection', LLFFCorrection, 'LSImagePaths', LSImagePaths, 'BackgroundPaths', BackgroundPaths, ...
    'minModifyTime', minModifyTime, 'maxModifyTime', maxModifyTime, 'maxWaitLoopNum', maxWaitLoopNum, ...
    'parseCluster', parseCluster, 'cpusPerTask', cpusPerTask, 'configFile', configFile, ...
    'mccMode', mccMode);


