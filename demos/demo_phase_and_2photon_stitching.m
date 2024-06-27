% demo to run stitching for non-light-sheet modalities (phase and 2-Photon) images
% 
% This may also apply to other microscopy modalities for both 2D and 3D images
%
% Note: The parameters demonstrated here are usually a subset of those available 
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


%% flat field correction for phase images (optional)
% here we use BaSiC (https://github.com/marrlab/BaSiC) to estimate and correct flat field 
% we also included flat field corrected images in 
% {destPath}/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Phase/ff_corrected/


% result folder:
% {destPath}/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Phase/ff_corrected_1/


currDataPath = [dataPath, '/Phase/'];
dir_info = dir([currDataPath, '*.tif']);
fns = cellfun(@(x) [currDataPath, x], {dir_info.name}', 'unif', 0);
outPath = sprintf('%s/ff_corrected_1/', currDataPath);
mkdir(outPath);
% background offset of the phase image
const = 6000;
% number of time points
nT = 2;
inpuFullpaths = cell(nT, 1);
outputFullpaths = cell(nT, 1);
func_strs = cell(nT, 1);
for i = 1 : nT
    fns_i = fns(contains(fns, sprintf('Iter_%04d', i - 1)));
    inpuFullpaths{i} = fns_i{1};
    [~, fsns_i] = fileparts(fns_i);
    if iscell(fsns_i)
        fsns_i = fsns_i{end};
    end
    outputFullpaths{i} = sprintf('%s/%s.tif', outPath, fsns_i);
    fns_strs = strjoin(fns_i, ''',''');
    
    XR_phase_image_flat_field_correction(fns_i, outPath, const);
end


%% generate image list for Phase images from tile positions (optional)
%
% We also include the image list csv file in
% {destPath}/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Phase/ImageList_all_timepoints.csv

% result file:
% {destPath}/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Phase/ImageList_from_tile_positions.csv

currDataPath = [dataPath, '/Phase/'];

% generation method
generationMethod = 'tile_position';

% channel pattern to include for image list generation
channelPatterns = {'Scan'};

% tile patterns for t, c, x, y, z to determine the time point, channel and
% tile indices. It must be a unique pattern for them. Typically, it should
% be a number combined with some string with or without underscore/dash.
% The program will use this pattern to determine the tile information for
% all the tiles with the given patterns. 
% time or channel patterns can be left as emtpy string if there is only one
% time point or one channel. For channel, it can also include cam, like CamA_ch0
tilePatterns = {'0000t', 'ch0', '000x', '000y', '000z'};

% xy pixel size
% the pixel size for the deskew/rotate data in the demo is isotropic at 0.108 um
xyPixelSize = 0.108;

% objective scan: scanning in the DS space
ObjectiveScan = false;

% inverted objective scan: scanning in the DSR space. 
IOScan = true;

% overlap size between tiles: axis order: yxz in pixels, xyz in um
overlapSize = [100, 100, 1];

% overlap size unit: pixel or um
overlapSizeType = 'pixel';

XR_generate_image_list_wrapper(currDataPath, generationMethod, xyPixelSize=xyPixelSize, ...
    channelPatterns=channelPatterns, tilePatterns=tilePatterns, overlapSize=overlapSize, ...
    overlapSizeType=overlapSizeType, ObjectiveScan=ObjectiveScan, IOScan=IOScan);


%% phase stitching
% The data is 2d and the coordinates are in the stage coordinates without any conversion
% so we use the 2D stitching scheme

% result folder:
% {destPath}/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Phase/ff_corrected/matlab_stitch_2d_phase/


currDataPath = [dataPath, 'Phase/ff_corrected/'];

% image list path: csv file
ImageListFullpath = [dataPath, 'Phase/ImageList_all_timepoints.csv'];

% stitch dir string, inside currDataPath
resultDirName = 'matlab_stitch_2d_phase';

% stitch in DS space
DS = false;
% stitch in Deskew/rotated space, if both DS and DSR false, stitch in skewed space
DSR = false;

% objective scan: scanning in the DS space
ObjectiveScan = false;

% inverted objective scan: scanning in the DSR space. 
IOScan = true;

% xy pixel size
% the pixel size for the deskew/rotate data in the demo is isotropic at 0.108 um
xyPixelSize = 0.108;

% scan direction. 
Reverse = true;

% check the setting files to see if a tile is flipped (bidirectional scan)
parseSettingFile = false;

% axis order
axisOrder = 'x,y,z';

% channels to stitch
ChannelPatterns = {'CamB_ch0'};

% stitch blending method, 'none': no blending, 'feather': feather blending
BlendMethod = 'feather';

% cross correlation registration, if false, directly stitch by the coordinates.
xcorrShift = true;
xcorrMode = 'primary';

% max allowed shift (in pixel) in xy axes between neighboring tiles in xcorr registration. 
xyMaxOffset = 100;

% max allowed shift in z (in pixel) axis between neighboring tiles in xcorr registration. 
zMaxOffset = 10;

% downsampling factors for overlap regions to calculate xcorr in xcorr registration. 
% with larger downsamplinf factors, the faster of the xcorr computing, yet
% lower accuracy will be. 
xcorrDownSample = [1, 1, 1];

% if true, save result as 16bit, if false, save as single
Save16bit = true;

% chunk size in zarr
blockSize = [512, 512, 1];

% batch size for stitching, should be multipler for block size
batchSize = [2048, 2048, 1];

% tile offset: counts add to the image, used when the image background is
% 0, to differentiate between image background and empty space. 
TileOffset = 0;

% erode the edge for given number of pixels
EdgeArtifacts = 0;

% here we use processed data that has already been flatfield correcetd, so we leave the processing function as empty.
processFunPath = '';

% if true, use slurm cluster for the computing; otherwise, use local machine
parseCluster = false;

XR_matlab_stitching_wrapper(currDataPath, ImageListFullpath, 'DS', DS, 'DSR', DSR, ...
    'ProcessedDirStr', '', 'Reverse', Reverse, 'ObjectiveScan', ObjectiveScan, ...
    'IOScan', IOScan, 'axisOrder', axisOrder, 'xcorrShift', xcorrShift, ...
    'xcorrMode', xcorrMode, 'xyMaxOffset', xyMaxOffset, 'zMaxOffset', zMaxOffset, ...
    'xcorrDownSample',xcorrDownSample, 'ChannelPatterns', ChannelPatterns, ...
    'resampleType', 'isotropic', 'xyPixelSize', xyPixelSize, 'BlendMethod', BlendMethod, ...
    'blockSize', blockSize, 'batchSize', batchSize, 'resultDirName', resultDirName, ...
    'parseSettingFile', parseSettingFile, 'processFunPath', processFunPath, ...
    'TileOffset', TileOffset, 'EdgeArtifacts', EdgeArtifacts, 'Save16bit', Save16bit, ...
    'parseCluster', parseCluster);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate image list for 2P images from tile positions

% result file:
% {destPath}/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/2P/ImageList_from_tile_positions.csv

currDataPath = [dataPath, '/2P/'];

% generation method
generationMethod = 'tile_position';

% channels to include for image list generation
channelPatterns = {'Tile'};

% tile patterns for t, c, x, y, z to determine the time point, channel and
% tile indices. It must be a unique pattern for them. Typically, it should
% be a number combined with some string with or without underscore/dash.
% The program will use this pattern to determine the tile information for
% all the tiles with the given patterns. 
% time or channel patterns can be left as emtpy string if there is only one
% time point or one channel. For channel, it can also include cam, like CamA_ch0
tilePatterns = {'0000t', 'ch0', '000x', '000y', '000z'};

% xy pixel size
% the pixel size for the deskew/rotate data in the demo is isotropic at 0.108 um
xyPixelSize = 0.108;

% scan step size
dz = 0.5;

% objective scan: scanning in the DS space
ObjectiveScan = false;

% inverted objective scan: scanning in the DSR space. 
IOScan = true;

% overlap size between tiles: axis order: yxz in pixels, xyz in um
overlapSize = [10, 10, 5];

% overlap size unit: pixel or um
overlapSizeType = 'um';

XR_generate_image_list_wrapper(currDataPath, generationMethod, xyPixelSize=xyPixelSize, ...
    dz=dz, channelPatterns=channelPatterns, tilePatterns=tilePatterns, overlapSize=overlapSize, ...
    overlapSizeType=overlapSizeType, ObjectiveScan=ObjectiveScan, IOScan=IOScan);


%% stitching of 2P image

% result folder:
% {destPath}/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/2P/matlab_stitch/


currDataPath = [dataPath, '/2P/'];

ImageListFullpath = [currDataPath, 'ImageList_from_tile_positions.csv'];

% stitch dir string, inside dataPath
resultDirName = 'matlab_stitch';

% stitch in DS space
DS = false;
% stitch in Deskew/rotated space, if both DS and DSR false, stitch in skewed space
DSR = false;

% objective scan: scanning in the DS space
ObjectiveScan = false;

% inverted objective scan: scanning in the DSR space. 
IOScan = true;

% scan direction. 
Reverse = true;

% intensally change xy pixelsize from 0.108 to 0.15 to count for the large
% overlaps between tiles. 
xyPixelSize = 0.15;

% scan step size
dz = 0.5;

% axis order: the axis order for x and z are actually inversed for this dataset
axisOrder = '-x,y,-z';

% primary channel for cross registration
PrimaryCh = 'ch0';

% channels to stitch
ChannelPatterns = {'ch0'};

% stitch blending method, 'none': no blending, 'feather': feather blending
BlendMethod = 'feather';

% cross correlation registration, if false, directly stitch by the coordinates.
xcorrShift = true;

% max allowed shift (in pixel) in xy axes between neighboring tiles in xcorr registration. 
xyMaxOffset = 200;

% max allowed shift in z (in pixel) axis between neighboring tiles in xcorr registration. 
zMaxOffset = 10;

% downsampling factors for overlap regions to calculate xcorr in xcorr registration. 
% with larger downsamplinf factors, the faster of the xcorr computing, yet
% lower accuracy will be. 
xcorrDownsample = [1, 1, 1];

% grid: only consider registration with direct neighbors
shiftMethod = 'grid';

% if true, save result as 16bit, if false, save as single
Save16bit = false;

% chunk size in zarr
blockSize = [256, 256, 256];

% batch size for stitching, should be multipler for block size
batchSize = [512, 512, 512];

% user defined processing function in tiff to zarr conversion, i.e., flat-field correction
% here is the example code to enable flat field correction
% define path to put the user defined function string (for the flat-field correction here).
tmpPath = [currDataPath, '/tmp/processFunc/'];
mkdir(tmpPath);
processFunPath = cell(numel(ChannelPatterns), 1);
% lower bound cap for flat-field
LowerLimit = 0.4;
% if true, subtract background image from flat-field image first;
% otherwise, use flat-field image as it is.
LSBackground = false;
% flat field image path
LSImagePaths = {
                    [currDataPath '/flatfield/flatfield.tif'], ...
                   };
% background image path
backgroundPaths = {
                [currDataPath, '/flatfield/background.tif']
                };
% offset to add after flat-field correction
TileOffset = 1;
% erode given number of pixels in the edge
EdgeArtifacts = 1;

dt = datetime('now', 'format', 'yyyymmddHHMMSS');
% number of channels
c = 1;
% define user-defined function
usrFun = sprintf("@(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'%s','%s',%d,%d,%s),%0.10d)", ...
    LSImagePaths{c}, backgroundPaths{c}, TileOffset, LowerLimit, string(LSBackground), EdgeArtifacts);

% save user defined function to path
fn = sprintf('%s/processFunction_c%04d_%s.mat', tmpPath, 1, dt);
save('-v7.3', fn, 'usrFun');
processFunPath{1} = fn;

% use slurm cluster if true, otherwise use the local machine (master job)
parseCluster = false;
% use master job for task computing or not. 
masterCompute = true;
% configuration file for job submission
configFile = '';
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;

XR_matlab_stitching_wrapper(currDataPath, ImageListFullpath, 'DS', DS, 'DSR', DSR, 'ProcessedDirStr', '', ...
    'Reverse', Reverse, 'ObjectiveScan', ObjectiveScan, 'IOScan', IOScan, ...
    'axisOrder', axisOrder, 'ChannelPatterns', ChannelPatterns, 'PrimaryCh', PrimaryCh, ...
    'xcorrShift', xcorrShift, 'xcorrDownsample', xcorrDownsample, 'shiftMethod', shiftMethod, ...
    'xyMaxOffset', xyMaxOffset, 'zMaxOffset', zMaxOffset,'xyPixelSize', xyPixelSize, ...
    'dz', dz, 'BlendMethod', BlendMethod, 'resultDirName', resultDirName, 'processFunPath', processFunPath, ...
    Save16bit=Save16bit, blockSize=blockSize, batchSize=batchSize, parseCluster=parseCluster, ...
    masterCompute=masterCompute, configFile=configFile, mccMode=mccMode);


