% demo to run geometric transformation

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


%% use the data wrapper to run deskew/rotation for all tiles in the data path
% if no resampling factor is provided, the output will be in isotropic
% voxel size as xyPixelSize.

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
    blockSize=blockSize, Save16bit=Save16bit, parseCluster=parseCluster, ...
    masterCompute=masterCompute, configFile=configFile, mccMode=mccMode);


%% compare separate deskew/rotation vs combined deskew/rotation

% read a tiff file
fsn = 'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t';
fn = [dataPath, fsn, '.tif'];

% fast Tiff reader
tic
im = readtiff(fn);
toc

%% deskew rotation
% size 1800 x 512 x 500

% conventional method (separate)
xyPixelSize = 0.108;
dz = 0.3;

fprintf('Conventional deskew/rotation: \n')
tic
XR_deskewRotateFrame_old(fn, xyPixelSize, dz, 'rotate', true, 'Reverse', true, 'Save16bit', true);
toc

% rename result folders 
movefile([dataPath, 'DS'], [dataPath, 'DS_separate']);
movefile([dataPath, 'DSR'], [dataPath, 'DSR_separate']);

% current method (combined)

fprintf('Combined deskew/rotation: \n')
tic
XR_deskewRotateFrame(fn, xyPixelSize, dz, 'rotate', true, 'Reverse', true, 'Save16bit', true);
toc

movefile([dataPath, 'DSR'], [dataPath, 'DSR_combined']);


%% deskew/rotate for larger data with more frames

% replicate the data to 2000 frames
im_rep = repmat(im, 1, 1, 4);
size(im_rep)

% save the data to disk
outPath = [dataPath, 'replicated/'];
mkdir(outPath);

% write to disk as a Tiff file
fnrep = [outPath, fsn, '_2k.tif'];

tic
writetiff(im_rep, fnrep);
toc
clear im_rep;


%% deskew rotation
% size 1800 x 512 x 2000
% note: please make sure there is enough memory for old method (>125 GB)

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/replicated/DSR/

% conventional method (separate)
xyPixelSize = 0.108;
dz = 0.3;

fprintf('Conventional deskew/rotation: \n')
tic
XR_deskewRotateFrame_old(fnrep, xyPixelSize, dz, 'rotate', true, 'Reverse', true, 'Save16bit', true);
toc

% rename result folders 
movefile([outPath, 'DS'], [outPath, 'DS_separate']);
movefile([outPath, 'DSR'], [outPath, 'DSR_separate']);

% current method (combined)
fprintf('Combined deskew/rotation: \n')
tic
XR_deskewRotateFrame(fnrep, xyPixelSize, dz, 'rotate', true, 'Reverse', true, 'Save16bit', true);
toc

movefile([outPath, 'DSR'], [outPath, 'DSR_combined']);


%% deskew/rotatin with even larger data for combined method (failed with old method in a desktop with 512 GB RAM)
% Please make sure your system has at least 60 GB RAM available

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/replicated/DSR/

% replicate the data to 5000 frames
im_rep = repmat(im, 1, 1, 10);
size(im_rep)

% save the data to disk
outPath = [dataPath, 'replicated/'];
mkdir(outPath);

% write to disk as a Tiff file
fnrep = [outPath, fsn, '_5k.tif'];

writetiff(im_rep, fnrep);
clear im_rep;

xyPixelSize = 0.108;
dz = 0.3;

tic
XR_deskewRotateFrame(fnrep, xyPixelSize, dz, 'rotate', true, 'Reverse', true, 'Save16bit', true);
toc

movefile([outPath, 'DSR'], [outPath, 'DSR_combined']);



