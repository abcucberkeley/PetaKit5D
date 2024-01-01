% demo to run microscope automatic pipeline for a single dataset with
% stitch

%% Step 1: set parameters for each step. 
% add the software to the path
cd(fileparts(which(mfilename)));
addpath(genpath('../'));

% data path
dataPath = '/clusterfs/fiona/Data/20200806_p35p4_Hex_Raptv_234-1_DLS/Exp01/';

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
Decon = true; % deconvolution
RotateAfterDecon = true; % rotate after deconvolution


% flat field parameters (if LLFFCorrection is set to true) 
LLFFCorrection = true;
% flat field image paths in the order of CamB_ch0 and CamA_ch0
LSImagePaths = {'/clusterfs/fiona/Data/20200803_LLCPKI/20200803_Calibration/FF/AVG_488nm.tif',...
'/clusterfs/fiona/Data/20200803_LLCPKI/20200803_Calibration/FF/AVG_560nm.tif'};

% camera background image paths in the order of CamB and CamA
BackgroundPaths = {'/clusterfs/fiona/OrcaDC/Aang/Orca_Aang_AVG_Cam_B.tif',...
    '/clusterfs/fiona/OrcaDC/Aang/Orca_Aang_AVG_Cam_A.tif'};


% stitch parameters
% image list csv file path
ImageListFullpaths = [dataPath, 'ImageList_Exp01.csv'];
% axis order for coordinates
axisOrder = '-x,y,z';


% deconvolution parameters
% by default, cpp decon is used 
% cpp decon path
cppDeconPath = '/global/home/groups/software/sl-7.x86_64/modules/RLDecon_CPU/20200718/build-cluster/cpuDeconv';
% if the dependency libraries are not loaded, we may also need to load the
% libraries. 
% loadModules = '';
cppDecon = true;
% if cppDecon and cudaDecon (false by default) are false, it uses matlab decon
% uncomment the line to use matlab decon. 
% cppDecon = false;

% psf full paths in the order of CamB_ch0 and CamA_ch0
PSFFullpaths = {'/clusterfs/fiona/Data/20200806_p35p4_Hex_Raptv_234-1_DLS/20200806_Calibration/PSF/Hex/TotPSF_488_CamB_3.tif', ...
'/clusterfs/fiona/Data/20200806_p35p4_Hex_Raptv_234-1_DLS/20200806_Calibration/PSF/Hex/TotPSF_560_CamA_5.tif'};
% background for deconvolution
Background = 100;
% If Stitch is true, it will use stitched data as input. Otherwise, use DS,
% DSR to specify use DS or DSR, or raw (both DS and DSR false) as input


%% Step 2: run the analysis with given parameters. 
XR_microscopeAutomaticProcessing(dataPath, 'xyPixelSize', xyPixelSize, 'dz', dz,  ...
    'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'Save16bit', Save16bit, ...
    'Overwrite', Overwrite, 'Streaming', Streaming, 'Deskew', Deskew, 'Rotate', Rotate, ...
    'Stitch', Stitch, 'Decon', Decon, 'cppDecon', cppDecon, 'cppDeconPath', cppDeconPath, ...
    'RotateAfterDecon', RotateAfterDecon, 'ImageListFullpaths', ImageListFullpaths, ...
    'axisOrder', axisOrder, 'PSFFullpaths', PSFFullpaths, 'Background', Background, ...
    'cpusPerTask', cpusPerTask);

