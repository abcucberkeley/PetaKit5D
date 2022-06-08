% demo to run RL deconvolution

%% Step 1: set parameters 
% add the software to the path
setup([],true);

% root path
rt = '/clusterfs/nvme/Data/rawData/';
% data path for data to be deconvolved, also support for multiple data folders
dataPaths = {[rt, '/']};

% xy pixel size in um
xyPixelSize = 0.108;
% z step size
dz = 0.5;
% scan direction
Reverse = true;
% psf z step size (we assume xyPixelSize also apply to psf)
dzPSF = 0.5;

% if true, check whether image is flipped in z using the setting files
parseSettingFile = false;

% channel patterns for the channels, the channel patterns should map the
% order of PSF filenames.
ChannelPatterns = {'CamA_ch0', ...
                   };  

% psf path
psf_rt = '/clusterfs/nvme/Data/PSFs/';            
PSFFullpaths = {
                [psf_rt, '/488nm.tif'], ...
                };            

% background to subtract
Background = 100;
% number of iterations
DeconIter = 80;
% decon to 80 iterations (not use the criteria for early stop)
fixIter = true;
% erode the edge after decon for number of pixels.
EdgeErosion = 8;
% save as 16bit, if false, save to single
Save16bit = false;
% use zarr file as input, if false, use tiff as input
zarrFile = false;
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

%% Step 2: run the deconvolution with given parameters. 
% the results will be saved in matlab_decon under the dataPaths. 
% the next step is deskew/rotate (if in skewed space for x-stage scan) or 
% rotate (if objective scan) or other processings. 

XR_decon_data_wrapper(dataPaths, 'xyPixelSize', xyPixelSize, 'dz', dz, 'Reverse', Reverse, ...
    'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, 'dzPSF', dzPSF, ...
    'parseSettingFile', parseSettingFile, 'Background', Background, 'CPPdecon', false, ...
    'CudaDecon', false, 'DeconIter', DeconIter, 'fixIter', fixIter, 'EdgeErosion', EdgeErosion, ...
    'Save16bit', Save16bit, 'zarrFile', zarrFile, 'Save16bit', ~false, 'parseCluster', parseCluster, ...
    'largeFile', largeFile,  'GPUJob', GPUJob, 'debug', debug, 'cpusPerTask', cpusPerTask);

