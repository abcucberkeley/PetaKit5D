% demo to run skewed space stitching
%
% Note: The parameters demonstrated here are usually a subset of those available 
% for the functions, with the rest using default values. For a comprehensive 
% list of parameters and their defaults, please see the function's parameter 
% list (or input parser) or refer to the parameter documentation (major_functions_documentation.txt).

%% Step 1: set parameters 
% add the software to the path
setup([]);

% data path
if ispc
    destPath = fullfile(getenv('USERPROFILE'), 'Downloads');   
    destPath = strrep(destPath, '\', '/');    
else
    destPath = '~/Downloads/';
end

dataPath = [destPath, '/PetaKit5D_demo_cell_image_dataset/'];

% image list path: csv file
% if not available, run stitch_generate_imagelist_from_encoder(dataPath, dz)
ImageListFullpath = [dataPath, 'ImageList_from_encoder.csv'];

% stitch in DS space
DS = false;
% stitch in Deskew/rotated space, if both DS and DSR false, stitch in skewed space
DSR = false;

% skew angle
SkewAngle = 32.45;

% scan direction. 
Reverse = true;

% xy pixel size and z scan step size
xyPixelSize = .108;
dz = 0.3;

% check the setting files to see if a tile is flipped (bidirectional scan)
parseSettingFile = false;

% axis order for the coordinates in the image list file
axisOrder = '-x,y,z';

% axis order for the image data
dataOrder = 'y,x,z';

% stitch blending method, 'none': no blending, 'feather': feather blending
BlendMethod = 'feather';

% stitch dir string, inside dataPath
stitchResultDir = 'matlab_stitch';

% cross correlation registration, if false, directly stitch by the coordinates.
xcorrShift = true;

% xcorr registration on the primary channel, and other channels uses the
% registration informaion from the primary channel
PrimaryCh = 'CamB_ch0';

% channels to stitch
ChannelPatterns = {'CamA_ch0', 'CamB_ch0'};

% use zarr file as input, if false, use tiff as input
zarrFile = false;

% if true, save result as 16bit, if false, save as single
Save16bit = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters below are used in special cases or for advanced fine tuning

% resample data, if empty, stitch in original resolution. If a 1X3 array,
% resample data and then stitch, mostly used for initial inspection
resampleFactor = [];

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

% batch size for stitching, should be multipler for block size
batchSize = [512, 512, 512];

% chunk size in zarr
blockSize = [256, 256, 256];

% user defined processing function in tiff to zarr conversion, i.e., flat-field correction
% here is the example code to enable flat field correction
% define path to put the user defined function string (for the flat-field correction here).
tmpPath = [dataPath, '/tmp/processFunc/'];
mkdir(tmpPath);
processFunPath = cell(numel(ChannelPatterns), 1);
% lower bound cap for flat-field
LowerLimit = 0.4;
% if true, subtract background image from flat-field image first;
% otherwise, use flat-field image as it is.
LSBackground = false;
% flat field image paths
LSImagePaths = {[dataPath, 'FF/averaged/ff_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_000x_000y_000z_0017t.tif'], ...
                [dataPath, 'FF/averaged/ff_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_000x_000y_000z_0017t.tif']};
% background image paths
BackgroundPaths = {[dataPath, 'FF/KorraFusions/AVG_DF400_CamA_10ms.tif'], ...
                   [dataPath, 'FF/KorraFusions/AVG_DF400_CamB_10ms.tif']};
% offset to add after flat-field correction
TileOffset = 1;
% date time to create a subfolder with the time stamp for the function string mat file.
dt = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));
for i = 1 : numel(ChannelPatterns)
    usrFun = sprintf("@(x)processFFCorrectionFrame(x,'%s','%s',%d,%d,%s)", ...
        LSImagePaths{i}, BackgroundPaths{i}, TileOffset, LowerLimit, string(LSBackground));

    fn = sprintf('%s/processFunction_c%04d_%s.mat', tmpPath, i, dt);
    save('-v7.3', fn, 'usrFun');
    processFunPath{i} = fn;
end

% Currently, flat field correction is disabled in the demo. Comment this line to enable it.
processFunPath = '';

% tile offset: counts add to the image, used when the image background is
% 0, to differentiate between image background and empty space. 
TileOffset = 0;

% use slurm cluster if true, otherwise use the local machine (master job)
parseCluster = false;
% use master job for task computing or not. 
masterCompute = true;
% configuration file for job submission
configFile = '';
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;


%% Step 2: run the stitching with given parameters. 
% the stitched results will be the zarr files in 'matlab_stitch_xcorr_feather_skewed_zarr'
% They are in skewed space, the next step is deconvolution/deskew rotation,
% or just deskew rotation. 

t0 = tic;
XR_matlab_stitching_wrapper(dataPath, ImageListFullpath, 'DS', DS, 'DSR', DSR, ...
    'SkewAngle', SkewAngle, 'xyPixelSize', xyPixelSize, 'dz', dz, 'Reverse', Reverse, ...
    'axisOrder', axisOrder, 'dataOrder', dataOrder, 'xcorrShift', xcorrShift, ...
    'ChannelPatterns', ChannelPatterns, 'PrimaryCh', PrimaryCh, 'zarrFile', zarrFile, ...
    'batchSize', batchSize, 'blockSize', blockSize, 'TileOffset', TileOffset, ...
    'resampleType', resampleType, 'resampleFactor', resampleFactor, 'BlendMethod', BlendMethod, ...
    'resultDir', stitchResultDir, 'xyMaxOffset', xyMaxOffset, 'zMaxOffset', zMaxOffset, ...
    'xcorrDownsample', xcorrDownsample, 'InputBbox', InputBbox, 'tileOutBbox', tileOutBbox, ...
    'parseSettingFile', parseSettingFile, 'processFunPath', processFunPath, ...
    'Save16bit', Save16bit, 'parseCluster', parseCluster, 'masterCompute', masterCompute, ...
    'configFile', configFile, 'mccMode', mccMode);
toc(t0);




