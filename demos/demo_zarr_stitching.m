% demo to run ZarrStitcher on both skewed data and deskew/rotated data.

clear, clc;
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


%% skewed space stitching
% the skewed space stitching is in another demo demo_skewed_space_stitching.m
% please refer to that demo for the details, here we directly run that demo

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_stitch/

demo_skewed_space_stitching


%% deskew/rotate stitched data

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_stitch/DSR/

dataPath_exps = {[dataPath, 'matlab_stitch/']};

% xy pixel size
xyPixelSize = 0.108;

% z scan step size
dz = 0.3;

% scan direction
Reverse = true;

% channel patterns to map the files for processing
ChannelPatterns = {'CamA', 'CamB'};

% if true, use large scale processing pipeline (split, process, and then merge)
largeFile = false;

% true if input is in zarr format
zarrFile = true;

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

XR_deskew_rotate_data_wrapper(dataPath_exps, xyPixelSize=xyPixelSize, dz=dz, ...
    Reverse=Reverse, ChannelPatterns=ChannelPatterns, largeFile=largeFile, ...
    zarrFile=zarrFile, saveZarr=saveZarr, blockSize=blockSize, Save16bit=Save16bit, ...
    parseCluster=parseCluster, masterCompute=masterCompute, configFile=configFile, ...
    mccMode=mccMode);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stitch in the DSR space (proper geometric space)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deskew/rotate data if DSR for tiles do not exist
% save results in Zarr format

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/DSR/

dataPath_exps = {dataPath};

% xy pixel size
xyPixelSize = 0.108;

% z scan step size
dz = 0.3;

% scan direction
Reverse = true;

% channel patterns to map the files for processing
ChannelPatterns = {'CamA', 'CamB'};

% if true, use large scale processing pipeline (split, process, and then merge)
largeFile = false;

% true if input is in zarr format
zarrFile = false;

% save output as zarr if true
saveZarr = true;

% block size to save the result 
BlockSize = [256, 256, 256];

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
    Reverse=Reverse, ChannelPatterns=ChannelPatterns, largeFile=largeFile, ...
    zarrFile=zarrFile, saveZarr=saveZarr, BlockSize=BlockSize, Save16bit=Save16bit, ...
    parseCluster=parseCluster, masterCompute=masterCompute, configFile=configFile, ...
    mccMode=mccMode);


%% stitching in DSR space

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_stitch_dsr/

% Step 1: set parameters 
% add the software to the path
setup([]);

% data path
dataPath = [destPath, '/LLSM5DTools_demo_cell_image_dataset/'];

% image list path: csv file
% if not available, run stitch_generate_imagelist_from_encoder(dataPath, dz)
ImageListFullpath = [dataPath, 'ImageList_from_encoder.csv'];

% stitch in DS space
DS = false;
% stitch in Deskew/rotated space, if both DS and DSR false, stitch in skewed space
DSR = true;

% directory string for the processed data within dataPath
ProcessedDirStr = 'DSR';

% true if input is in Zarr format
zarrFile = true;

% skew angle
SkewAngle = 32.45;

% scan direction. 
Reverse = true;

% resolution [xyPixelsize, dz]
% the pixel size for the deskew/rotate data in the demo is isotropic at 0.108 um
Resolution = [0.108, 0.108];

% check the setting files to see if a tile is flipped (bidirectional scan)
parseSettingFile = false;

% axis order
axisOrder = '-x,y,z';

% stitch blending method, 'none': no blending, 'feather': feather blending
BlendMethod = 'feather';

% stitch dir string, inside dataPath
stitchResultDir = 'matlab_stitch_dsr';

% cross correlation registration, if false, directly stitch by the coordinates.
xcorrShift = true;

%if true, only stitch first time point
onlyFirstTP = false;

% stitch pipeline, no need to change, zarr pipeline is the mostly used one.
stitchPipeline = 'zarr';

% xcorr registration on the primary channel, and other channels uses the
% registration informaion from the primary channel
PrimaryCh = 'CamB_ch0';

% channels to stitch
ChannelPatterns = {'CamA_ch0', 'CamB_ch0'};

% if true, save result as 16bit, if false, save as single
Save16bit = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters below are used in special cases or for advanced fine tuning

% resample data, if empty, stitch in original resolution. If a 1X3 array,
% resample data and then stitch, mostly used for initial inspection
resample = [];

% resample type: isotropic, xy_isotropic; effective when 'resample' is
% not set. Keep it as 'isotropic'. 
resampleType = 'isotropic';

% max allowed shift (in pixel) in xy axes between neighboring tiles in xcorr registration. 
xyMaxOffset = 150;

% max allowed shift in z (in pixel) axis between neighboring tiles in xcorr registration. 
zMaxOffset = 150;

% downsampling factors for overlap regions to calculate xcorr in xcorr registration. 
% with larger downsamplinf factors, the faster of the xcorr computing, yet
% lower accuracy will be. 
xcorrDownsample = [2, 2, 2];

% crop input tile before any processing, empty (no crop) or a 1X6 array [ymin, xmin, zmin, ymax, xmax, zmax]
InputBbox = [];

% crop input tile after processing, empty (no crop) or a 1X6 array [ymin, xmin, zmin, ymax, xmax, zmax]
tileOutBbox = [];

% chunk size in zarr
blockSize = [256, 256, 256];

% user defined processing function in tiff to zarr conversion, i.e., flat field correction
processFunPath = '';

% tile offset: counts add to the image, used when the image background is
% 0, to differentiate between image background and empty space. 
TileOffset = 0;

% if true, use slurm cluster for the computing; otherwise, use local machine
parseCluster = false;


% Step 2: run the stitching with given parameters. 
% the stitched results will be the zarr files in 'matlab_stitch_xcorr_feather_skewed_zarr'
% They are in skewed space, the next step is deconvolution/deskew rotation,
% or just deskew rotation. 

tic
XR_matlab_stitching_wrapper(dataPath, ImageListFullpath, 'DS', DS, 'DSR', DSR, ...
    'ProcessedDirStr', ProcessedDirStr, 'zarrFile', zarrFile, 'SkewAngle', SkewAngle, ...
    'Reverse', Reverse, 'axisOrder', axisOrder, 'xcorrShift', xcorrShift,  ...
    'ChannelPatterns', ChannelPatterns, 'PrimaryCh', PrimaryCh, 'blockSize', blockSize, ...
    'TileOffset', TileOffset, 'resampleType', resampleType, 'resample', resample, ...
    'Resolution', Resolution,'BlendMethod', BlendMethod, 'resultDir', stitchResultDir, ...
    'onlyFirstTP', onlyFirstTP, 'xyMaxOffset', xyMaxOffset, 'zMaxOffset', zMaxOffset, ...
    'xcorrDownsample', xcorrDownsample,  'InputBbox', InputBbox, 'tileOutBbox', tileOutBbox, ...
    'pipeline', stitchPipeline,  'parseSettingFile', parseSettingFile, ...
    'processFunPath', processFunPath, 'Save16bit', Save16bit, 'parseCluster', false);
toc


%% 





