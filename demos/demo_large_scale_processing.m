% demo to run large scale processing

% Note: Large-scale stitching with large tiles, large-scale deconvolution (inplace) 
% and deskew/rotation only works for zarr input. Large-scale stitching with small
% tiles, and large-scale deconvolution (inmemory option) works for tiff files 
% that can fit to memory. 

clear, clc;

fprintf('Large-scale processing demo...\n\n');

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


%% Step 2: stitching in the skewd space
% the code is in demo_skewed_space_stitching.m
%
% note: if you would like to test the slurm cluster running of the
% stitching, please update the parseCluster, configFile, mccMode in
% demo_skewed_space_stitching.m

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_stitch/

demo_skewed_space_stitching


%% Step 3: large-scale deconvolution

%% Step 3.1: test the parameter for OMW backward projector

% Cam A
psfFn = [dataPath, 'PSF/560_c.tif'];
% OTF thresholding parameter
OTFCumThresh = 0.9;
% true if the PSF is in skew space
skewed = true;
XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed);


%% 

% Cam B
psfFn = [dataPath, 'PSF/488_2_c.tif'];
% OTF thresholding parameter
OTFCumThresh = 0.9;
% true if the PSF is in skew space
skewed = true;
XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed);


%% Step 3.2: set parameters 
% add the software to the path
setup([]);

% root path
rt = dataPath;
% data path for data to be deconvolved, also support for multiple data folders
dataPaths = {[rt, 'matlab_stitch/']};

% xy pixel size in um
xyPixelSize = 0.108;
% z step size
dz = 0.3;
% scan direction
Reverse = true;
% psf z step size (we assume xyPixelSize also apply to psf)
dzPSF = 0.3;

% if true, check whether image is flipped in z using the setting files
parseSettingFile = false;

% channel patterns for the channels, the channel patterns should map the
% order of PSF filenames.
ChannelPatterns = {'CamA_ch0', ...
                   'CamB_ch0', ...
                   };  

% psf path
psf_rt = rt;            
PSFFullpaths = {
                [psf_rt, 'PSF/560_c.tif'], ...    
                [psf_rt, 'PSF/488_2_c.tif'], ...
                };            

% RL method
RLmethod = 'omw';
% wiener filter parameter for CamA and CamB, respectively
wienerAlpha = [0.008, 0.005];
% OTF thresholding parameter
OTFCumThresh = 0.9;
% true if the PSF is in skew space
skewed = true;
% deconvolution result path string (within dataPath)
deconPathstr = 'matlab_decon_omw';

% background to subtract
Background = 100;
% number of iterations
DeconIter = 2;
% decon to 80 iterations (not use the criteria for early stop)
fixIter = true;
% erode the edge after decon for number of pixels.
EdgeErosion = 20;
% save as 16bit, if false, save to single
Save16bit = true;
% use zarr file as input, if false, use tiff as input
zarrFile = true;
% number of cpu cores
cpusPerTask = 24;
% use cluster computing for different images
parseCluster = false;
% set it to true for large files that cannot be fitted to RAM/GPU, it will
% split the data to chunks for deconvolution
largeFile = true;
% large method: "inplace": in place decon (only load the region with border buffer for decon);
%               "inmemory": for data can load to memory (load whole the
%               data, and decon a small region with border buffer). 
largeMethod = 'inplace';
% batch size to define each basic region for deconvolution (without
% including border buffer), typically as multipler of blockSize, and can
% fit to GPU if use GPU. 
% This size fits to GPUs with 24 GB or more vRAM for the demo data.
batchSize = [512, 768, 768];
% block size to define the zarr chunk size
blockSize = [256, 256, 256];

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

%% Step 3.3: run the deconvolution with given parameters. 
% the results will be saved in matlab_decon under the dataPaths. 
% the next step is deskew/rotate (if in skewed space for x-stage scan) or 
% rotate (if objective scan) or other processings. 

% If using GPU, please make sure the batch can fit to GPU; otherwise, reduce batchSize

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_stitch/matlab_decon_omw/

XR_decon_data_wrapper(dataPaths, 'deconPathstr', deconPathstr, 'xyPixelSize', xyPixelSize, ...
    'dz', dz, 'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, ...
    'dzPSF', dzPSF, 'parseSettingFile', parseSettingFile, 'RLmethod', RLmethod, ...
    'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'skewed', skewed, ...
    'Background', Background, 'CPPdecon', false, 'CudaDecon', false, 'DeconIter', DeconIter, ...
    'fixIter', fixIter, 'EdgeErosion', EdgeErosion, 'Save16bit', Save16bit, ...
    'zarrFile', zarrFile, 'batchSize', batchSize, 'blockSize', blockSize, ...
    'parseCluster', parseCluster, 'largeFile', largeFile, 'largeMethod', largeMethod, ...
    'GPUJob', GPUJob, 'debug', debug, 'cpusPerTask', cpusPerTask, 'ConfigFile', ConfigFile, ...
    'GPUConfigFile', GPUConfigFile, 'mccMode', mccMode);

% release GPU if using GPU computing
if GPUJob && gpuDeviceCount('available') > 0
    reset(gpuDevice);
end


%% Step 4: large-scale deskew/rotation

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_stitch/matlab_decon_omw/DSR/

% data path for data to be deconvolved, also support for multiple data folders
dataPath_exps = {[dataPath, 'matlab_stitch/matlab_decon_omw/']};

% xy pixel size
xyPixelSize = 0.108;
% z scan step size
dz = 0.3;
% scan direction
Reverse = true;
% channel patterns to map the files for processing
ChannelPatterns = {'CamA', 'CamB'};

% if true, use large scale processing pipeline (split, process, and then merge)
largeFile = true;
% true if input is in zarr format
zarrFile = true;
% save output as zarr if true
saveZarr = true;
% batch size for individual task, only the size in y is used, and it should
% be the multiplier of blocksize in y. Also need to adjust accordingly
% based on the available memory 
BatchSize = [512, 512, 512];
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

XR_deskew_rotate_data_wrapper(dataPath_exps, xyPixelSize=xyPixelSize, dz=dz, ...
    Reverse=Reverse, largeFile=largeFile, zarrFile=zarrFile, saveZarr=saveZarr, ...
    BatchSize=BatchSize, blockSize=blockSize, Save16bit=Save16bit, parseCluster=parseCluster, ...
    masterCompute=masterCompute, configFile=configFile, mccMode=mccMode);


