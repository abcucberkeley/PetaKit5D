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
% Note: please make sure there is enough memory for conventional method (>60 GB)

xyPixelSize = 0.108;
dz = 0.3;

% remove old results
if exist([dataPath, 'DS/', fsn, '.tif'], 'file') 
    delete([dataPath, 'DS/', fsn, '.tif']);
end
if exist([dataPath, 'DSR/', fsn, '.tif'], 'file') 
    delete([dataPath, 'DSR/', fsn, '.tif']);
end

fprintf('\nConventional deskew/rotation for 500 frames: \n')
tic
XR_deskewRotateFrame_old(fn, xyPixelSize, dz, 'rotate', true, 'Reverse', true, 'Save16bit', true);
toc

% move results to method specific DS/DSR folders 
if exist([dataPath, 'DS_separated'], 'dir') 
    rmdir([dataPath, 'DS_separated'], 's');
end
if exist([dataPath, 'DSR_separated'], 'dir') 
    rmdir([dataPath, 'DSR_separated'], 's');
end
mkdir([dataPath, 'DS_separated']);
mkdir([dataPath, 'DSR_separated']);
movefile([dataPath, 'DS/', fsn, '.tif'], [dataPath, 'DS_separated']);
movefile([dataPath, 'DSR/', fsn, '.tif'], [dataPath, 'DSR_separated']);

% current method (combined)

fprintf('Combined deskew/rotation for 500 frames: \n')
tic
XR_deskewRotateFrame(fn, xyPixelSize, dz, 'rotate', true, 'Reverse', true, 'Save16bit', true);
toc

% move results to method specific DSR folder
if exist([dataPath, 'DSR_combined'], 'dir')
    rmdir([dataPath, 'DSR_combined'], 's');
end
mkdir([dataPath, 'DSR_combined'])
movefile([dataPath, 'DSR/', fsn, '.tif'], [dataPath, 'DSR_combined']);


%% deskew/rotate for larger data with more frames

% replicate the data to 2000 frames
nframe = 2000;
fprintf('\nReplicate data to %d frames... ', nframe);
im_rep = repmat(im, 1, 1, ceil(nframe / 500));
if size(im_rep, 3) > nframe
    im_rep = crop3d_mex(im_rep, [1, 1, 1, size(im_rep, 1), size(im_rep, 2), nframe]);
end
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
% Note: please make sure there is enough memory for conventional method (>125 GB)

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/replicated/DSR/

% conventional method (separate)
xyPixelSize = 0.108;
dz = 0.3;

fprintf('\nConventional deskew/rotation for %d frames: \n', nframe)
tic
XR_deskewRotateFrame_old(fnrep, xyPixelSize, dz, 'rotate', true, 'Reverse', true, 'Save16bit', true);
toc

% move the files to method specific DS/DSR folders
if ~exist([outPath, 'DS_separate/'], 'dir') ||  ~exist([outPath, 'DSR_separate/'], 'dir')
    mkdir([outPath, 'DS_separate/']);
    mkdir([outPath, 'DSR_separate/']);
end
movefile([outPath, 'DS/', fsn, '_2k.tif'], [outPath, 'DS_separate/']);
movefile([outPath, 'DSR/', fsn, '_2k.tif'], [outPath, 'DSR_separate/']);

% current method (combined)
fprintf('Combined deskew/rotation for %d frames: \n', nframe)
tic
XR_deskewRotateFrame(fnrep, xyPixelSize, dz, 'rotate', true, 'Reverse', true, 'Save16bit', true);
toc

% move the file to method specific DSR folder
if ~exist([outPath, 'DSR_combined/'], 'dir')
    mkdir([outPath, 'DSR_combined/']);
end
movefile([outPath, 'DSR/', fsn, '_2k.tif'], [outPath, 'DSR_combined/']);


%% deskew/rotatin with even larger data for combined method (failed with old method in a desktop with 512 GB RAM)
% Note: please make sure your system has at least 60 GB RAM available

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/replicated/DSR/

% replicate the data to 5000 frames
nframe = 5000;
fprintf('\nReplicate data to %d frames... ', nframe);
im_rep = repmat(im, 1, 1, ceil(nframe / 500));
if size(im_rep, 3) > nframe
    im_rep = crop3d_mex(im_rep, [1, 1, 1, size(im_rep, 1), size(im_rep, 2), nframe]);
end
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

fprintf('\nCombined deskew/rotation for %d frames: \n\n', nframe)
tic
XR_deskewRotateFrame(fnrep, xyPixelSize, dz, 'rotate', true, 'Reverse', true, 'Save16bit', true);
toc

% move the file to method specific DSR folder
movefile([outPath, 'DSR/', fsn, '_5k.tif'], [outPath, 'DSR_combined']);



