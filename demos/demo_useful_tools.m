% demo to run some useful tools

% Note: useful tools for image processing, support tiff/zarr, and large-scale 
% image data. If you would like to test the slurm cluster running of the tools, 
% please update the parseCluster, configFile, mccMode.
%
% The parameters demonstrated here are usually a subset of those available 
% for the functions, with the rest using default values. For a comprehensive 
% list of parameters and their defaults, please see the function's parameter 
% list (or input parser) or refer to the parameter documentation (major_functions_documentation.txt).
%
%
% List of tools in the demo: 
%   resample
%   crop
%   max intensity projection(MIP)
%   FFT analysis
%   PSF analysis
%   Fourier shell correlation (FSC) analysis, 
%   Tiff to Zarr converter
%   Zarr to Tiff converter
%   Imaris converter


clear, clc;

fprintf('Useful tools demo...\n\n');

% move to the PetaKit5D root directory
curPath = pwd;
if ~endsWith(curPath, 'PetaKit5D')
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


%% get our demo data from zenodo/Dropbox (skip this step if the data is already downloaded)
% download the example dataset from zenodo (https://doi.org/10.5281/zenodo.1492027) manually, 
% or use the code below to download the data from Dropbox
if ispc
    destPath = fullfile(getenv('USERPROFILE'), 'Downloads');
    destPath = strrep(destPath, '\', '/');    
else
    destPath = '~/Downloads/';
end
demo_data_downloader(destPath);

dataPath = [destPath, '/PetaKit5D_demo_cell_image_dataset/'];


%% run deconvolution that are used for some tools (or running the deconvolution demo)
% use OMW deconvolution 

% add the software to the path
setup([]);

% root path
rt = dataPath;
% data path for data to be deconvolved, also support for multiple data folders
dataPaths = {[rt, '/']};

% xy pixel size in um
xyPixelSize = 0.108;
% z step size
dz = 0.3;
% psf z step size (we assume xyPixelSize also apply to psf)
dzPSF = 0.3;

% if true, check whether image is flipped in z using the setting files
parseSettingFile = false;

% channel patterns for the channels, the channel patterns should map the
% order of PSF filenames.
ChannelPatterns = {'CamB_ch0', ...
                   };  

% psf path
psf_rt = rt;            
PSFFullpaths = {
                [psf_rt, 'PSF/488_2_c.tif'], ...
                };            

% RL method
RLmethod = 'omw';
% wiener filter parameter
% alpha parameter should be adjusted based on SNR and data quality.
% typically 0.002 - 0.01 for SNR ~20; 0.02 - 0.1 or higher for SNR ~7
wienerAlpha = 0.005;
% OTF thresholding parameter
OTFCumThresh = 0.9;
% true if the PSF is in skew space
skewed = true;
% deconvolution result path string (within dataPath)
resultDirName = 'matlab_decon_omw';

% background to subtract
Background = 100;
% number of iterations
DeconIter = 2;
% decon to 80 iterations (not use the criteria for early stop)
fixIter = true;
% erode the edge after decon for number of pixels.
EdgeErosion = 0;
% save as 16bit; if false, save to single
Save16bit = true;
% use zarr file as input; if false, use tiff as input
zarrFile = false;
% save output as zarr file; if false,s ave as tiff
saveZarr = false;
% number of cpu cores
cpusPerTask = 24;
% use cluster computing for different images
parseCluster = false;
% set it to true for large files that cannot be fitted to RAM/GPU, it will
% split the data to chunks for deconvolution
largeFile = false;
% use GPU for deconvolution
GPUJob = false;
% if true, save intermediate results every 5 iterations.
debug = false;
% config file for the master jobs that runs on CPU node
ConfigFile = '';
% config file for the GPU job scheduling on GPU node
GPUConfigFile = '';
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_decon_omw/

XR_decon_data_wrapper(dataPaths, 'resultDirName', resultDirName, 'xyPixelSize', xyPixelSize, ...
    'dz', dz, 'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, ...
    'dzPSF', dzPSF, 'parseSettingFile', parseSettingFile, 'RLmethod', RLmethod, ...
    'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'skewed', skewed, ...
    'Background', Background, 'DeconIter', DeconIter, 'fixIter', fixIter, ...
    'EdgeErosion', EdgeErosion, 'Save16bit', Save16bit, 'zarrFile', zarrFile, ...
    'saveZarr', saveZarr, 'parseCluster', parseCluster, 'largeFile', largeFile, ...
    'GPUJob', GPUJob, 'debug', debug, 'cpusPerTask', cpusPerTask, 'ConfigFile', ConfigFile, ...
    'GPUConfigFile', GPUConfigFile, 'mccMode', mccMode);

% release GPU if using GPU computing
if GPUJob && gpuDeviceCount('available') > 0
    reset(gpuDevice);
end


%% run deskew/rotation to generate deskew/rotate results that are used for some tools

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_decon_omw/DS_for_useful_tools/
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_decon_omw/DSR_for_useful_tools/

dataPaths = {[dataPath, 'matlab_decon_omw/']};

% deskew result path string (within dataPath)
DSDirName = 'DS_for_useful_tools/';
% deskew/rotate result path string (within dataPath)
DSRDirName = 'DSR_for_useful_tools/';

% run deskew if true
Deskew = true;
% run rotation if true
Rotate = true;
% use combined processing if true. Note: if you only need deskew without
% rotation, set both Rotate and DSRCombined as false.
DSRCombined = false;
% xy pixel size
xyPixelSize = 0.108;
% z scan step size
dz = 0.3;
% Skew angle
skewAngle = 32.45;
% if true, objective scan; otherwise, sample scan
ObjectiveScan = false;
% scan direction
Reverse = true;
% channel patterns to map the files for processing
ChannelPatterns = {'CamA', 'CamB'};

% flat field related parameters
% disable flat field for this demo
FFCorrection = false;
% lower bound cap for flat field
LowerLimit = 0.4;
% constant offset after flat field correction
constOffset = 1;
% flat field image paths
FFImagePaths = {[dataPath, 'FF/averaged/ff_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_000x_000y_000z_0017t.tif'], ...
                [dataPath, 'FF/averaged/ff_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_000x_000y_000z_0017t.tif']};
% background image paths
BackgroundPaths = {[dataPath, 'FF/KorraFusions/AVG_DF400_CamA_10ms.tif'], ...
                   [dataPath, 'FF/KorraFusions/AVG_DF400_CamB_10ms.tif']};

% if true, use large scale processing pipeline (split, process, and then merge)
largeFile = false;
% true if input is in zarr format
zarrFile = false;
% save output as zarr if true
saveZarr = true;
% block size to save the result 
blockSize = [256, 256, 256];
% save output as uint16 if true
Save16bit = true;

% use slurm cluster if true, otherwise use the local machine (master job)
parseCluster = false;
% use master job for task computing or not. 
masterCompute = true;
% configuration file for job submission
configFile = '';
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;

XR_deskew_rotate_data_wrapper(dataPaths, DSDirName=DSDirName, DSRDirName=DSRDirName, ...
    Deskew=Deskew, Rotate=Rotate, DSRCombined=DSRCombined, xyPixelSize=xyPixelSize, ...
    dz=dz, SkewAngle=skewAngle, ObjectiveScan=ObjectiveScan, Reverse=Reverse, ...
    ChannelPatterns=ChannelPatterns, FFCorrection=FFCorrection, LowerLimit=LowerLimit, ...
    constOffset=constOffset, FFImagePaths=FFImagePaths, BackgroundPaths=BackgroundPaths, ...
    largeFile=largeFile, zarrFile=zarrFile, saveZarr=saveZarr, blockSize=blockSize, ...
    Save16bit=Save16bit, parseCluster=parseCluster, masterCompute=masterCompute, ...
    configFile=configFile, mccMode=mccMode);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% resampling tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resample data according to the resample factor in one or multiple folders. 

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/resampled/

fprintf('\nResampling tool... \n\n');

dataPaths = {[dataPath, '']};

% resampling factor: order yxz, downsample 2x if it is 2. 
resampleFactor = [2, 2, 2];

% result path string (within dataPath)
resultDirName = 'resampled';
% channel patterns to map the files for processing
channelPatterns = {'CamA_ch0', 'CamB_ch0'};
% interploation method: linear, cubic, nearest
interpMethod = 'linear';
% save as 16 bit if true
save16bit = true;
% true if input is in zarr format
zarrFile = false;
% true for large files that the processing cannot fit into memory
largeFile = false;
% save output as zarr if true
saveZarr = false;
% block size to define the zarr chunk size, only needed for zarr file
blockSize = [256, 256, 256];
% batch size to define the size for each task, only needed for large zarr file
batchSize = [512, 512, 512];
% border buffer size for the resampling, only needed for large zarr file
borderSize = [5, 5, 5];

% use cluster computing for different images
parseCluster = false;
% master compute
masterCompute = true;
% number of cpu cores
cpusPerTask = 4;
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;
% config file for the master jobs that runs on CPU node
configFile = '';

XR_resample_dataset(dataPaths, resampleFactor, resultDirName=resultDirName, ...
    channelPatterns=channelPatterns, interpMethod=interpMethod, save16bit=save16bit, ...
    zarrFile=zarrFile, largeFile=largeFile, saveZarr=saveZarr, blockSize=blockSize, ...
    batchSize=batchSize, borderSize=borderSize, parseCluster=parseCluster, ...
    masterCompute=masterCompute, cpusPerTask=cpusPerTask, mccMode=mccMode, ...
    configFile=configFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% crop tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crop images with given bounding box in one or multiple folders. 

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/Cropped/

fprintf('\nCropping tool... \n\n');

dataPaths = {[dataPath, '']};

% bounding box to crop: [ymin, xmin, zmin, ymax, xmax, zmax]
inputBbox = [201, 101, 101, 1200, 412, 400];

% result path string (within dataPath)
resultDirName = 'Cropped';
% channel patterns to map the files for processing
channelPatterns = {'CamA_ch0', 'CamB_ch0'};
% true if input is in zarr format
zarrFile = false;
% true for large files that the processing cannot fit into memory
largeFile = false;
% save output as zarr if true
saveZarr = false;
% block size to define the zarr chunk size, only needed for zarr file
blockSize = [256, 256, 256];
% save as 16 bit if true
save16bit = true;

% use cluster computing for different images
parseCluster = false;
% master compute
masterCompute = true;
% number of cpu cores
cpusPerTask = 4;
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;
% config file for the master jobs that runs on CPU node
configFile = '';

XR_crop_dataset(dataPaths, inputBbox, resultDirName=resultDirName, channelPatterns=channelPatterns, ...
    zarrFile=zarrFile, largeFile=largeFile, saveZarr=saveZarr, blockSize=blockSize, ...
    save16bit=save16bit, parseCluster=parseCluster, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, mccMode=mccMode, configFile=configFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MIP tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run max projection in one or multiple folders. 

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/MIPs/

fprintf('\nMax intensity projection tool... \n\n');

dataPaths = {[dataPath, '']};

% result path string (within dataPath)
resultDirName = 'MIPs';
% axes for the MIP: xyz
axis = [1, 1, 1];
% channel patterns to map the files for processing
channelPatterns = {'CamA_ch0', 'CamB_ch0'};
% true if input is in zarr format
zarrFile = false;
% true for large files that the processing cannot fit into memory
largeFile = false;
% batch size for a task, only needed for zarr file
batchSize = [1024, 1024, 1024];
% save as 16 bit if true
save16bit = true;

% use cluster computing for different images
parseCluster = false;
% master compute
masterCompute = true;
% number of cpu cores
cpusPerTask = 4;
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;
% config file for the master jobs that runs on CPU node
configFile = '';

XR_MIP_wrapper(dataPaths, resultDirName=resultDirName, axis=axis, channelPatterns=channelPatterns, ...
    zarrFile=zarrFile, largeFile=largeFile, batchSize=batchSize, save16bit=save16bit, ...
    parseCluster=parseCluster, masterCompute=masterCompute, cpusPerTask=cpusPerTask, ...
    mccMode=mccMode, configFile=configFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FFT analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute fourier power spectrum of the given images, by default only the central 
% slides are saved with gamma 0.5. If set save3DStack as true, also save 3D 
% stack with log10 applied.

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_decon_omw/DS_for_useful_tools/FFT/

fprintf('\nFourier space analysis tool... \n\n');

dataPaths = {[dataPath, 'matlab_decon_omw/DS_for_useful_tools/']};

% result path string (within dataPaths)
resultDirName = 'FFT';
% xy pixel size
xyPixelSize = 0.108;
% z scan step size (the input data is deskewed, so the z step size is
% original dz * sind(theta). 
dz = 0.3 * sind(32.45);
% true if input is in zarr format
zarrFile = true;
% output size for FFT spectrum
outSize = [1001, 1001, 1001];
% output voxel size, isotropic, here we use 0.1 medium wavelength
outPixelSize = 0.1 * 0.488 / 1.33;
% save 3D spectrum if true.
save3DStack = false;

% use cluster computing for different images
parseCluster = false;
% master compute
masterCompute = true;
% number of cpu cores
cpusPerTask = 24;
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;
% config file for the master jobs that runs on CPU node
configFile = '';

XR_fftSpectrumComputingWrapper(dataPaths, 'resultDirName', resultDirName, ...
    'xyPixelSize', xyPixelSize, 'dz', dz, 'zarrFile', zarrFile, 'outSize', outSize, ...
    'outPixelSize', outPixelSize, 'save3DStack', save3DStack, parseCluster=parseCluster, ...
    masterCompute=masterCompute, cpusPerTask=cpusPerTask, mccMode=mccMode, ...
    configFile=configFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSF analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform psf analysis and visualization for the real and fourier spaces, and
% compare the patterns with theoretical widefield PSF (Richards & Wolf)

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/PSF/DS/PSFAnalysis/

fprintf('\nPSF analysis tool... \n\n');

dataPaths = {[dataPath, 'PSF/']};

% xy pixel size
xyPixelSize = 0.108;
% z step size
dz = 0.3;
% Skew angle
skewAngle = 32.45;
% channel patterns to map the files for processing
ChannelPatterns = {'488'};
% channel wavelength
Channels = [488];
% deskew PSF if true. This is for PSF in the skewed space
Deskew = true;
% if true, objective scan; otherwise, sample scan
ObjectiveScan = false;

% Richards & Wolf theoretical PSF path for corresponding channels
RWFn = {[dataPath, 'PSF/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif'], ...
        };
% source description
sourceStr = 'Demo PSF';
% if true display the figures
visible = true;

% use cluster computing for different images
parseCluster = false;
% master compute
masterCompute = true;
% number of cpu cores
cpusPerTask = 24;
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;
% config file for the master jobs that runs on CPU node
configFile = '';

XR_psf_analysis_wrapper(dataPaths, 'xyPixelSize', xyPixelSize, 'dz', dz, ...
    'skewAngle', skewAngle, 'ChannelPatterns', ChannelPatterns, 'Channels', Channels, ...
    'Deskew', Deskew, 'ObjectiveScan', ObjectiveScan, 'sourceStr', sourceStr, ...
    'RWFn', RWFn, visible=visible, parseCluster=parseCluster, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, mccMode=mccMode, configFile=configFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FSC analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run deconvolution with different number of iterations as the demostration of FSC analysis

% add the software to the path
setup([]);

% root path
rt = dataPath;
% data path for data to be deconvolved, also support for multiple data folders
dataPaths = {[rt, '/']};


% xy pixel size in um
xyPixelSize = 0.108;
% z step size
dz = 0.3;
% psf z step size (we assume xyPixelSize also apply to psf)
dzPSF = 0.3;

% if true, check whether image is flipped in z using the setting files
parseSettingFile = false;

% channel patterns for the channels, the channel patterns should map the
% order of PSF filenames.
% we only deconvolve one tile for FSC computing and visualization
ChannelPatterns = {'Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t', ...
                   };  

% psf path
psf_rt = rt;            
PSFFullpaths = {
                [psf_rt, 'PSF/488_2_c.tif'], ...
                };            

% RL method
RLmethod = 'omw';
% wiener filter parameter
% alpha parameter should be adjusted based on SNR and data quality.
% typically 0.002 - 0.01 for SNR ~20; 0.02 - 0.1 or higher for SNR ~7
% here we use larger alpha for more stable results with more iterations in
% the debug mode
wienerAlpha = 0.05;
% OTF thresholding parameter
OTFCumThresh = 0.9;
% true if the PSF is in skew space
skewed = true;
% deconvolution result path string (within dataPath)
resultDirName = 'matlab_decon_omw_debug';

% background to subtract
Background = 100;
% number of iterations
DeconIter = 20;
% decon to 2 iterations
fixIter = true;
% erode the edge after decon for number of pixels.
EdgeErosion = 0;
% save as 16bit; if false, save to single
Save16bit = true;
% use zarr file as input; if false, use tiff as input
zarrFile = false;
% save output as zarr file; if false,s ave as tiff
saveZarr = false;
% number of cpu cores
cpusPerTask = 24;
% use cluster computing for different images
parseCluster = false;
% set it to true for large files that cannot be fitted to RAM/GPU, it will
% split the data to chunks for deconvolution
largeFile = false;
% use GPU for deconvolution
GPUJob = false;
% save Step 
saveStep = 2;
% if true, save intermediate results every #saveStep iterations.
debug = true;
% config file for the master jobs that runs on CPU node
ConfigFile = '';
% config file for the GPU job scheduling on GPU node
GPUConfigFile = '';
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_decon_omw/

XR_decon_data_wrapper(dataPaths, 'resultDirName', resultDirName, 'xyPixelSize', xyPixelSize, ...
    'dz', dz, 'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, ...
    'dzPSF', dzPSF, 'parseSettingFile', parseSettingFile, 'RLmethod', RLmethod, ...
    'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'skewed', skewed, ...
    'Background', Background, 'DeconIter', DeconIter, 'fixIter', fixIter, ...
    'EdgeErosion', EdgeErosion, 'Save16bit', Save16bit, 'zarrFile', zarrFile, ...
    'saveZarr', saveZarr, 'parseCluster', parseCluster, 'largeFile', largeFile, ...
    'GPUJob', GPUJob, 'saveStep', saveStep, 'debug', debug, 'cpusPerTask', cpusPerTask, ...
    'ConfigFile', ConfigFile, 'GPUConfigFile', GPUConfigFile, 'mccMode', mccMode);

% release GPU if using GPU computing
if GPUJob && gpuDeviceCount('available') > 0
    reset(gpuDevice);
end


%% FSC analysis
% perform fourier shell correlation (FSC) analysis for images in one or
% multiple folders. Here we run FSC analysis for the deconvolution with different 
% iterations to demonstrate how to use FSC to determine the optimal number
% of iterations.  

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_decon_omw_debug/Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t_debug/FSCs
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_decon_omw_debug/Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t_debug/figures

fprintf('\nFourier Shell Correlation (FSC) analysis tool... \n\n');

dataPaths = {[dataPath, 'matlab_decon_omw_debug/Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t_debug/']};

% result path string (within dataPaths)
resultDirName = 'FSCs';
% xy pixel size
xyPixelSize = 0.108;
% z scan step size (the input data is deskewed, so the z step size is
% original dz * sind(theta). 
dz = 0.3 * sind(32.45);
% radial interval for the calculation of correlation
dr = 5;
% angular interval for the calculation of correlation
dtheta = pi / 6;
% resolution thresholding method: one-bit, half-bit, or fixed 
resThreshMethod = 'one-bit';
% resolution threshold (only for fixed method)
resThresh = 0.2;
% half size of region for FSC (typically leave some extra space for border)
% must be the same across all three dimensions and better to be odd number
halfSize = [201, 201, 201];
% input bounding box to define the region for FSC computing. If it is
% empty, use the central region for the computing. 
inputBbox = [];
% channel patterns to map the files for processing
ChannelPatterns = {'Iter'};
% channel wavelength (only for figure naming purpose)
Channels = [488];
% suffix (only for figure naming purpose)
suffix = 'decon';
% iteration interval for figure plotting
iterInterval = 2;

% use cluster computing for different images
parseCluster = false;
% master compute
masterCompute = true;
% number of cpu cores
cpusPerTask = 24;
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;
% config file for the master jobs that runs on CPU node
configFile = '';

XR_FSC_analysis_wrapper(dataPaths, 'resultDirName', resultDirName, 'xyPixelSize', xyPixelSize, ...
    'dz', dz, 'dr', dr, 'dtheta', dtheta, 'resThreshMethod', resThreshMethod, ...
    'resThresh', resThresh, 'halfSize', halfSize, 'inputBbox', inputBbox, ...
    'ChannelPatterns', ChannelPatterns, 'Channels', Channels, 'suffix', suffix, ...
    'iterInterval', iterInterval, parseCluster=parseCluster, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, mccMode=mccMode, configFile=configFile)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tiff to zarr conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert tiff to zarr files in one or multiple folders

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/tiff_to_zarr/

fprintf('\nTiff to Zarr conversion tool... \n\n');

dataPaths = {[dataPath, '/']};

% result path string (within dataPaths)
resultDirName = 'tiff_to_zarr';
% channel patterns to map the files for processing
channelPatterns = {'CamA_ch0', 'CamB_ch0'};
% use cluster computing for different images
% bounding box to crop the files in the conversion: [ymin, xmin, zmin, ymax, xmax, zmax]
inputBbox = [];

parseCluster = false;
% master compute
masterCompute = true;
% number of cpu cores
cpusPerTask = 2;
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;
% config file for the master jobs that runs on CPU node
configFile = '';

XR_tiffToZarr_wrapper(dataPaths, resultDirName=resultDirName, channelPatterns=channelPatterns, ...
    inputBbox=inputBbox, parseCluster=parseCluster, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, mccMode=mccMode, configFile=configFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  zarr to tiff conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert zarr to tiff files in one or multiple folders

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_decon_omw/DSR_for_useful_tools/tiffs/

fprintf('\nZarr to Tiff conversion tool... \n\n');

dataPaths = {[dataPath, 'matlab_decon_omw/DSR_for_useful_tools/']};

% result path string (within dataPaths)
resultDirName = 'tiffs';
% channel patterns to map the files for processing
channelPatterns = {'CamA_ch0', 'CamB_ch0'};
% use cluster computing for different images

parseCluster = false;
% master compute
masterCompute = true;
% number of cpu cores
cpusPerTask = 2;
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;
% config file for the master jobs that runs on CPU node
configFile = '';

XR_zarrToTiff_wrapper(dataPaths, resultDirName=resultDirName, channelPatterns=channelPatterns, ...
    parseCluster=parseCluster, masterCompute=masterCompute, cpusPerTask=cpusPerTask, ...
    mccMode=mccMode, configFile=configFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Imaris converter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert zarr or tiff files to imaris file in one or multiple folders. It
% supports multiple channels and multiple time points. If there are
% multiple files in the same channel, it assumes they are time series data
% with the order of filenames. 

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_decon_omw/DSR_for_useful_tools/imaris/

fprintf('\nTiff/Zarr to Imaris conversion tool... \n\n');

dataPaths = {[dataPath, 'matlab_decon_omw/DSR_for_useful_tools/']};

% result path string (within dataPaths)
resultDirName = 'imaris';
% since the demo dataset has several tiles, we only convert the time series
% for one tile. 
ChannelPatterns = {'000x_003y_000z'};
% voxel sizes for the data: y, x, z
pixelSizes = [0.108, 0.108, 0.108];
% use zarr file as input if true
zarrFile = true;
% block size for the imaris file
blockSize = [64, 64, 64];

parseCluster = false;
% master compute
masterCompute = true;
% number of cpu cores
cpusPerTask = 24;
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;
% config file for the master jobs that runs on CPU node
configFile = '';

XR_imaris_conversion_data_wrapper(dataPaths, 'ChannelPatterns', ChannelPatterns, ...
    'pixelSizes', pixelSizes, 'zarrFile', zarrFile, 'blockSize', blockSize, ...
    parseCluster=parseCluster, masterCompute=masterCompute, cpusPerTask=cpusPerTask, ...
    mccMode=mccMode, configFile=configFile);

