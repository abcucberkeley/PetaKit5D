%% demo to demonstrate how to setup the deconvolution for widefield and confocal images
% This may also apply to other microscopy modalities

% Note: the OMW method may need more than 2 or 3 iterations to converge,
% like 10-100 iterations with relative larger wiener factor (0.01 - 1). 
% Confocal images may have similar parameter settings as light sheet images.
%
% The parameters demonstrated here are usually a subset of those available 
% for the functions, with the rest using default values. For a comprehensive 
% list of parameters and their defaults, please see the function's parameter 
% list (or input parser) or refer to the parameter documentation (major_functions_documentation.txt).


clear, clc;

fprintf('Stitching demo for Phase and 2-Photon images...\n\n');

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


%% Step 1: get our demo data from zenodo/Dropbox (skip this step if the data is already downloaded)
% download the example dataset from zenodo (https://doi.org/10.5281/zenodo.11500863) manually, 
% or use the code below to download the data from Dropbox

if ispc
    destPath = fullfile(getenv('USERPROFILE'), 'Downloads');
    destPath = strrep(destPath, '\', '/');
else
    destPath = '~/Downloads/';
end
demo_data_downloader(destPath, 'other_modalities');

dataPath = [destPath, '/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% deconvolution of the widefield image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step A2: test the parameters for OMW backward projector

psfFn = [dataPath, 'Widefield/WF_PSF.tif'];
% OTF thresholding parameter
OTFCumThresh = 0.65;
% true if the PSF is in skew space
skewed = false;
% minimum internsity threshold for the OTF segmentation
minIntThrsh = 1e-3;
XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed, minIntThrsh=minIntThrsh);


%% Step A3: OMW deconvolution 
%% Step A3.1: set parameters 
% add the software to the path
setup([]);

% root path
rt = dataPath;
% data path for data to be deconvolved, also support for multiple data folders
dataPaths = {[rt, 'Widefield/']};

% xy pixel size in um
xyPixelSize = 0.157;
% z step size
dz = 0.1;
% scan direction
Reverse = true;
% psf z step size (we assume xyPixelSize also apply to psf)
dzPSF = 0.1;

% if true, check whether image is flipped in z using the setting files
parseSettingFile = false;

% channel patterns for the channels, the channel patterns should map the
% order of PSF filenames.
ChannelPatterns = {'WF_raw', ...
                   };  

% psf path
psf_rt = rt;            
PSFFullpaths = {
                [psf_rt, 'Widefield/WF_PSF.tif'], ...
                };            

% RL method
RLmethod = 'omw';
% wiener filter parameter
% alpha parameter should be adjusted based on SNR and data quality.
% typically 0.002 - 0.01 for SNR ~20; 0.02 - 0.1 or higher for SNR ~7
wienerAlpha = 0.02;
% OTF thresholding parameter
OTFCumThresh = 0.65;
% true if the PSF is in skew space
skewed = true;
% hann window range applied to the distance transform, 0.0 means the center and 1.0 means border of OTF mask
hanWinBounds = [0.8, 1.2];
% deconvolution result path string (within dataPath)
resultDirName = 'matlab_decon_omw';

% background to subtract
Background = 100;
% number of iterations
DeconIter = 20;
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
GPUJob = true;
% if true, save intermediate results every 5 iterations.
debug = false;
% config file for the master jobs that runs on CPU node
ConfigFile = '';
% config file for the GPU job scheduling on GPU node
GPUConfigFile = '';
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;


%% Step A3.2: run the deconvolution with given parameters. 
% the results will be saved in matlab_decon under the dataPaths. 
% the next step is deskew/rotate (if in skewed space for x-stage scan) or 
% rotate (if objective scan) or other processings. 

% result folder:
% {destPath}/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Widefield/matlab_decon_omw/

XR_decon_data_wrapper(dataPaths, 'resultDirName', resultDirName, 'xyPixelSize', xyPixelSize, ...
    'dz', dz, 'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, ...
    'dzPSF', dzPSF, 'parseSettingFile', parseSettingFile, 'RLmethod', RLmethod, ...
    'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'skewed', skewed, ...
    'hanWinBounds', hanWinBounds, 'Background', Background, 'DeconIter', DeconIter, ...
    'fixIter', fixIter, 'EdgeErosion', EdgeErosion, 'Save16bit', Save16bit, ...
    'zarrFile', zarrFile, 'saveZarr', saveZarr, 'parseCluster', parseCluster, ...
    'largeFile', largeFile, 'GPUJob', GPUJob, 'debug', debug, 'cpusPerTask', cpusPerTask, ...
    'ConfigFile', ConfigFile, 'GPUConfigFile', GPUConfigFile, 'mccMode', mccMode);

% release GPU if using GPU computing
if GPUJob && gpuDeviceCount('available') > 0
    reset(gpuDevice);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% deconvolution of the Confocal image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step B2: test the parameters for OMW backward projector

psfFn = [dataPath, 'Confocal/Confocal_PSF.tif'];
% OTF thresholding parameter
OTFCumThresh = 0.875;
% true if the PSF is in skew space
skewed = false;
% minimum internsity threshold for the OTF segmentation
minIntThrsh = 1e-3;
XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed, minIntThrsh=minIntThrsh);


%% Step B3: OMW deconvolution 
%% Step B3.1: set parameters 
% add the software to the path
setup([]);

% root path
rt = dataPath;
% data path for data to be deconvolved, also support for multiple data folders
dataPaths = {[rt, 'Confocal/']};

% xy pixel size in um
xyPixelSize = 0.157;
% z step size
dz = 0.1;
% scan direction
Reverse = true;
% psf z step size (we assume xyPixelSize also apply to psf)
dzPSF = 0.1;

% if true, check whether image is flipped in z using the setting files
parseSettingFile = false;

% channel patterns for the channels, the channel patterns should map the
% order of PSF filenames.
ChannelPatterns = {'Confocal_raw', ...
                   };  

% psf path
psf_rt = rt;            
PSFFullpaths = {
                [psf_rt, 'Confocal/Confocal_PSF.tif'], ...
                };            

% RL method
RLmethod = 'omw';
% wiener filter parameter
% alpha parameter should be adjusted based on SNR and data quality.
% typically 0.002 - 0.01 for SNR ~20; 0.02 - 0.1 or higher for SNR ~7
wienerAlpha = 0.004;
% OTF thresholding parameter
OTFCumThresh = 0.875;
% true if the PSF is in skew space
skewed = true;
% deconvolution result path string (within dataPath)
resultDirName = 'matlab_decon_omw';

% background to subtract
Background = 100;
% number of iterations
DeconIter = 3;
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


%% Step B3.2: run the deconvolution with given parameters. 
% the results will be saved in matlab_decon under the dataPaths. 
% the next step is deskew/rotate (if in skewed space for x-stage scan) or 
% rotate (if objective scan) or other processings. 

% result folder:
% {destPath}/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Confocal/matlab_decon_omw/

XR_decon_data_wrapper(dataPaths, 'resultDirName', resultDirName, 'xyPixelSize', xyPixelSize, ...
    'dz', dz, 'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, ...
    'dzPSF', dzPSF, 'parseSettingFile', parseSettingFile, 'RLmethod', RLmethod, ...
    'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'skewed', skewed, ...
    'Background', Background, 'DeconIter', DeconIter, 'fixIter', fixIter, 'EdgeErosion', EdgeErosion, ...
    'Save16bit', Save16bit, 'zarrFile', zarrFile, 'saveZarr', saveZarr, 'parseCluster', parseCluster, ...
    'largeFile', largeFile, 'GPUJob', GPUJob, 'debug', debug, 'cpusPerTask', cpusPerTask, ...
    'ConfigFile', ConfigFile, 'GPUConfigFile', GPUConfigFile, 'mccMode', mccMode);

% release GPU if using GPU computing
if GPUJob && gpuDeviceCount('available') > 0
    reset(gpuDevice);
end

