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
dataPath = {'/clusterfs/fiona/Data/20200806_p35p4_Hex_Raptv_234-1_DLS/Exp01/', ...
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
% number of cores for the pipeline (only for ABC cluster). For other
% environment, this setting will be ignored. 
cpusPerTask = 24;

% define which step is included in the analysis
Deskew = true; % deskew
Rotate = true; % rotate after deskew
Stitch = true; % stitch
Decon = false; % deconvolution
RotateAfterDecon = false; % rotate after deconvolution


% flat field parameters (if LLFFCorrection is set to true) 
LLFFCorrection = true;
% flat field image paths in the order of CamB_ch0 and CamA_ch0
LSImagePaths = {'/clusterfs/fiona/Data/20200803_LLCPKI/20200803_Calibration/FF/AVG_488nm.tif',...
                '/clusterfs/fiona/Data/20200803_LLCPKI/20200803_Calibration/FF/AVG_560nm.tif'};

% camera background image paths in the order of CamB and CamA
BackgroundPaths = {'/clusterfs/fiona/OrcaDC/Aang/Orca_Aang_AVG_Cam_B.tif',...
                   '/clusterfs/fiona/OrcaDC/Aang/Orca_Aang_AVG_Cam_A.tif'};

% stitch parameters
% image list csv file path for each dataset
ImageListFullpaths = {[dataPath{1}, 'ImageList_Exp01.csv'], ...
                      [dataPath{2}, 'ImageList_Exp02.csv'], ...
                      [dataPath{3}, 'ImageList_Exp03.csv']};
% axis order for coordinates
axisOrder = '-x,y,z';


%% Step 2: run the analysis with given parameters. 
XR_microscopeAutomaticProcessing(dataPath, 'xyPixelSize', xyPixelSize, 'dz', dz,  ...
    'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'Save16bit', Save16bit, ...
    'Overwrite', Overwrite, 'Streaming', Streaming, 'Deskew', Deskew, 'Rotate', Rotate, ...
    'Stitch', Stitch, 'Decon', Decon, 'RotateAfterDecon', RotateAfterDecon, ...
    'ImageListFullpaths', ImageListFullpaths, 'axisOrder', axisOrder, 'cpusPerTask', cpusPerTask);

