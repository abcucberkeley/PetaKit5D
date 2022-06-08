% demo to run skewed space stitching

%% Step 1: set parameters 
% add the software to the path
setup([],true);

% data path
dataPath = '/clusterfs/fiona/Data/20220115_Korra_LLCPK_LFOV_0p1PSAmpKan/run1/';

% image list path: csv file
% if not available, run stitch_generate_imagelist_from_encoder(dataPath, dz)
ImageListFullpath = '/clusterfs/fiona/Data/20220115_Korra_LLCPK_LFOV_0p1PSAmpKan/run1/ImageList_from_encoder.csv';

% stitch in DS space
DS = false;
% stitch in Deskew/rotated space, if both DS and DSR false, stitch in skewed space
DSR = false;

% skew angle
SkewAngle = 32.45;

% scan direction. 
Reverse = true;

% resolution [xyPixelsize, dz]
Resolution = [0.108, 0.30];

% check the setting files to see if a tile is flipped (bidirectional scan)
parseSettingFile = false;

% axis order
axisOrder = '-x,y,z';

% stitch blending method, 'none': no blending, 'feather': feather blending
BlendMethod = 'feather';

% stitch dir string, inside dataPath
stitchResultDir = 'matlab_stitch_xcorr_feather_skewed_zarr_test';

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

% crop input tile after processing, empty (no crop) or a 1X6 array [ysize, xsize, zsize]
CropToSize = [];

% chunk size in zarr
blockSize = [256, 256, 256];

% user defined processing function in tiff to zarr conversion, i.e., flat field correction
processFunPath = '';

% tile offset: counts add to the image, used when the image background is
% 0, to differentiate between image background and empty space. 
TileOffset = 0;


%% Step 2: run the stitching with given parameters. 
% the stitched results will be the zarr files in 'matlab_stitch_xcorr_feather_skewed_zarr'
% They are in skewed space, the next step is deconvolution/deskew rotation,
% or just deskew rotation. 

tic
XR_matlab_stitching_wrapper(dataPath, ImageListFullpath, 'useExistDecon', ~true, ...
    'DS', DS, 'DSR', DSR, 'SkewAngle', SkewAngle, 'Reverse', Reverse, 'axisOrder', axisOrder, ...
    'xcorrShift', xcorrShift, 'ChannelPatterns', ChannelPatterns, 'PrimaryCh', PrimaryCh, ...
    'blockSize', blockSize, 'TileOffset', TileOffset, 'resampleType', resampleType, ...
    'resample', resample, 'Resolution', Resolution,'BlendMethod', BlendMethod, ...
    'resultDir', stitchResultDir, 'onlyFirstTP', onlyFirstTP, 'xyMaxOffset', xyMaxOffset, ...
    'zMaxOffset', zMaxOffset, 'xcorrDownsample', xcorrDownsample,  'InputBbox', InputBbox, ...
    'CropToSize', CropToSize,  'pipeline', stitchPipeline,  'parseSettingFile', parseSettingFile, ...
    'processFunPath', processFunPath, 'Save16bit', Save16bit, 'parseCluster', false);
toc




