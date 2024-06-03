% demo to run ZarrStitcher on both skewed data and deskew/rotated data.

clear, clc;

fprintf('Stitching demo...\n\n');

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
% download the example dataset from zenodo (https://doi.org/10.5281/zenodo.10471978) manually, 
% or use the code below to download the data from Dropbox
if ispc
    destPath = fullfile(getenv('USERPROFILE'), 'Downloads'); 
    destPath = strrep(destPath, '\', '/');
else
    destPath = '~/Downloads/';
end
demo_data_downloader(destPath);

dataPath = [destPath, '/PetaKit5D_demo_cell_image_dataset/'];


%% skewed space stitching
% the skewed space stitching is in another demo demo_skewed_space_stitching.m
% please refer to that demo for the details, here we directly run that demo
%
% note: if you would like to test the slurm cluster running of the
% stitching, please update the parseCluster, configFile, mccMode in
% demo_skewed_space_stitching.m

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_stitch/

demo_skewed_space_stitching


%% deskew/rotate stitched data

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_stitch/DSR/

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
% {destPath}/PetaKit5D_demo_cell_image_dataset/DSR/

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
% {destPath}/PetaKit5D_demo_cell_image_dataset/matlab_stitch_dsr/

% Step 1: set parameters 
% add the software to the path
setup([]);

% data path
dataPath = [destPath, '/PetaKit5D_demo_cell_image_dataset/'];

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

% xy pixel size
% the pixel size for the deskew/rotate data in the demo is isotropic at 0.108 um
xyPixelSize = 0.108;

% z scan step size
dz = 0.3;

% scan direction. 
Reverse = true;

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
    'xyPixelSize', xyPixelSize, 'dz', dz, 'Reverse', Reverse, 'axisOrder', axisOrder, ...
    'xcorrShift', xcorrShift, 'ChannelPatterns', ChannelPatterns, 'PrimaryCh', PrimaryCh, ...
    'blockSize', blockSize, 'TileOffset', TileOffset, 'resampleType', resampleType, ...
    'resampleFactor', resampleFactor, 'BlendMethod', BlendMethod, 'resultDir', stitchResultDir, ...
    'xyMaxOffset', xyMaxOffset, 'zMaxOffset', zMaxOffset, ...
    'xcorrDownsample', xcorrDownsample, 'InputBbox', InputBbox, 'tileOutBbox', tileOutBbox, ...
    'parseSettingFile', parseSettingFile, ...
    'processFunPath', processFunPath, 'Save16bit', Save16bit, 'parseCluster', false);
toc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stitching benchmarks with BigStitcher and Sitching-spark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, we demonstrate the benchmark setup using the demo dataset with BigStitcher
% and Stitching-spark.

% To do son, we need to perform deskewing and rotation on the input data. 
% The deskewed and rotated data will be used for stitching, along with 
% registration information required by BigStitcher. Assuming the previous 
% two sessions have already completed, the deskewed/rotated data is located 
% at {destPath}/DSR/, and the stitched result is at {destPath}/matlab_stitch_dsr/


fprintf(['The following steps will not be automatically run, and requires the ', ...
         'setup of BigStitcher-Spark and Stitching-Spark.\n']);  
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup of stitching with Spark version of BigStitcher
% Step A1. Install Apache Spark on all systems you wish to use: https://spark.apache.org/downloads.html
% In our benchmarks, we installed Spark 3.4.1 with java 11.0.0.1 
% Step A2. Install BigSitcher-Spark: https://github.com/JaneliaSciComp/BigStitcher-Spark
%% Step A3. Convert the deskew/rotated data from Zarr to N5 using the code below:
% Step A3.1. Add python enviroment to Matlab with zarr installed, same way as in
% demo_fast_tiff_zarr_readers_writers.m
pythonPath = '/path/to/your/python';
setup([], true, pythonPath);

%% Step A3.2 Convert Zarr files to N5 format
dataPath = [destPath, 'PetaKit5D_demo_cell_image_dataset/'];
DSRPath = [dataPath, '/DSR/'];
%   only include the first time point of CamB
channel = 'Scan_Iter_0000_0000_CamB_ch0';
blockSize = [256, 256, 256];
masterCompute = true;

dir_info = dir([DSRPath, '*', channel, '*zarr']);
fsns = {dir_info.name}';
nF = numel(fsns);

for f = 1 : nF
    fsn = fsns{f};
    fn = [DSRPath, fsn]
    
    zarrToN5(fn, blockSize=blockSize, masterCompute=masterCompute);
end

%% Step A3.3. generate xml file for BigStitcher
% Set the follwoing variables as needed
fileName = 'test_cell_data_DSR_CamB.n5';
outFile = 'test_cell_data_DSR_CamB_with_xcorr.xml';
dir_info = dir(DSRPath);
currDir = [dir_info(1).folder, '/'];
fn = [currDir, fileName];
outFn = [currDir, outFile];

csvFile = [dataPath, 'ImageList_from_encoder.csv'];
voxelSize = [0.108, 0.108, 0.108];
ChannelPatterns = {'Scan_Iter_0000_0000_CamB_ch0'};
zarrFile = true;
axisOrder = '-xyz';
dataOrder = 'yxz';

xcorrInfoFn = [dataPath, 'matlab_stitch_dsr/xcorr/Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106126432msecAbs.mat'];

[s, t, voxelSize_n5] = constructXmlFile(currDir, fn,csvFile,voxelSize, ChannelPatterns, zarrFile, axisOrder, dataOrder, xcorrInfoFn);
writestruct(s,outFn,"StructNodeName","SpimData");
removeTextTags(outFn);
cleanXmlFile(outFn);

%% Step A3.4. reorganize data via symbolic links as the input for BigStitcher
% Set the follwoing variables as needed
mainFolder = [currDir, 'test_cell_data_DSR_CamB.n5'];
n5Path = [currDir, 'n5'];

nFiles = size(t,1);
for i = 1:nFiles    
    inFile = [n5Path '/' t.Filename{i}(1:end-4) 'n5/'];
    outFile = char(sprintf("%s/setup%d/timepoint0/",mainFolder,i-1));
    mkdir(outFile);

    % symlink
    system(['ln -s "' inFile '" "' outFile 's0"']);

    timePointAttrFile = sprintf("%s/setup%d/attributes.json",mainFolder,i-1);
    fileID = fopen(timePointAttrFile,'w');
    fprintf(fileID,"{""downsamplingFactors"":[[1,1,1]],""dataType"":""uint16""}");

    resAttrFile = sprintf("%s/setup%d/timepoint0/attributes.json",mainFolder,i-1);
    fileID = fopen(resAttrFile,'w');
    fprintf(fileID,"{""resolution"":[%f,%f,%f],""saved_completely"":true,""multiScale"":true}",voxelSize_n5);
end

%% Step A4: run BigStitcher using the Slurm cluster
% Step A4.1: Start your master spark node and initialize your worker nodes
%{
module load blosc/1.21.5
module load spark/3.4.1
module load java/11.0.0.1
start-master.sh
%}
% example code in slurm submit script initWorker.sh with one node
%{
#!/bin/sh
#SBATCH --qos=abc_high
#SBATCH --partition=abc
#SBATCH --account=co_abc
#SBATCH --nodes=1
#SBATCH --time=36:00:00
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=21000
#SBATCH --array=1

export SPARK_NO_DAEMONIZE=1

module purge
module load spark/3.4.1
module load java/11.0.0.1

start-worker.sh spark://10.17.209.11:7077
%}
% run the code below in terminal
%{
sbatch initWorker.sh
%}
% Step A4.2: Submit jobs to your spark cluster using net.preibisch.bigstitcher.spark.AffineFusion
% and BigStitcher-Spark/target/BigStitcher-Spark-0.0.2-SNAPSHOT.jar in the terminal
% Main Spark CMD
%{
destPath='/string/for/destPath/'
spark-submit --master spark://10.17.209.11:7077 \
            --driver-cores 2 \
            --driver-memory 2G \
            --executor-memory 500G \
            --executor-cores 24 \
            --class net.preibisch.bigstitcher.spark.AffineFusion \
             /your/BigStitcher/Installation/path/BigStitcher-Spark/target/BigStitcher-Spark-0.0.2-SNAPSHOT.jar \
            -x  ${destPath}/PetaKit5D_demo_cell_image_dataset/DSR/test_cell_data_DSR_CamB_with_xcorr.xml \
            -o ${destPath}/PetaKit5D_demo_cell_image_dataset/DSR/test_cell_data_DSR_CamB.n5 -d '/ch488/s0' \
            --UINT16 --minIntensity 1 --maxIntensity 65535 --channelId 0
%}


%% Step A5: convert n5 to zarr (optional)
% if you would like to convert the stitched result to zarr, run the code below

% add   ,"n5":"2.0.0"   in the attributes.json before conversion
fn = [currDir, 'test_cell_data_DSR_CamB.n5/ch488/s0'];
flipEmptyValue = true;
masterCompute = true;
N5ToZarr(fn, flipEmptyValue=flipEmptyValue, masterCompute=masterCompute);

% output zarr is in 
% {destPath}/PetaKit5D_demo_cell_image_dataset/DSR/test_cell_data_DSR_CamB.n5/ch488/zarr/s0.zarr


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup of stitching with Stitching-Spark
% Step B1. Install Apache Spark on all systems you wish to use: https://spark.apache.org/downloads.html
% In our benchmarks, we installed Spark 3.4.1 with java 11.0.0.1 and Blosc 1.21.0 (compatible with Zarr)
% Step B2. Install stitching-spark: https://github.com/saalfeldlab/stitching-spark/tree/master
%% Step B3. Run DSR on the demo data with demo_geometric_transformation.m and save the result as tiff format. 
% below is the code to run the deskew/rotation for CamB time point 0

dataPath = [destPath, 'PetaKit5D_demo_cell_image_dataset/'];
dataPath_exps = {dataPath};
% xy pixel size
xyPixelSize = 0.108;
% z scan step size
dz = 0.3;
% scan direction
Reverse = true;
% channel patterns to map the files for processing
ChannelPatterns = {'Scan_Iter_0000_0000_CamB_ch0'};

% if true, use large scale processing pipeline (split, process, and then merge)
largeFile = false;
% true if input is in zarr format
zarrFile = false;
% save output as zarr if true
saveZarr = false;
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
    zarrFile=zarrFile, saveZarr=saveZarr, BlockSize=blockSize, Save16bit=Save16bit, ...
    parseCluster=parseCluster, masterCompute=masterCompute, configFile=configFile, ...
    mccMode=mccMode);


%% Step B4. Generate the image list for the DSR data

csvFile = [dataPath, 'ImageList_from_encoder.csv'];
tmp = readtable(csvFile, 'Delimiter', 'comma');
% only keep CamB in time point 0
tmp = tmp(contains(tmp.Filename, '0000_0000_CamB'), :);

outFn = [dataPath, 'DSR/ImageList_from_encoder_t0_CamB.csv'];
writetable(tmp, outFn);


%% Step B5. Run the stitching-spark conversion script to generate a json file from the image list
% example (all the following is a one line cmd): python ./startup-scripts/spark-local/parse-imagelist-metadata.py 
% -b ./PetaKit5D_demo_cell_image_dataset/DSR/ -a x,y,z -i ./PetaKit5D_demo_cell_image_dataset/DSR/ImageList_from_encoder_t0_CamB.csv 
% -r .108,.108,.108
% Step B6. Start your master spark node and initialize your worker nodes
%{
module load spark/3.4.1
module load java/11.0.0.1
start-master.sh
%}
% example code in slurm submit script initWorker.sh with one node
%{
#!/bin/sh 
#SBATCH --qos=abc_high
#SBATCH --partition=abc
#SBATCH --account=co_abc
#SBATCH --nodes=1
#SBATCH --time=36:00:00
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=21000
#SBATCH --array=1

export SPARK_NO_DAEMONIZE=1

module purge
module load spark/3.4.1
module load java/11.0.0.1

start-worker.sh spark://10.17.209.11:7077
%}
% run the code below in terminal
%{
sbatch initWorker.sh
%}
% Step B7. Follow the instructions on the stitching spark wiki: https://github.com/saalfeldlab/stitching-spark/wiki
% Step B8. Flatfield can be skipped. Using the wiki instructions you can
% submit your Stitching, N5 export, and tiff conversion jobs to your spark cluster

% Flatfield Spark CMD
%{
stitching_install_path='/installation/path/for/stitching-spark'
destPath='/string/for/destPath/'
spark-submit --master spark://10.17.209.11:7077 \
            --driver-cores 2 \
            --driver-memory 2G \
            --executor-memory 500G \
            --executor-cores 24 \
            --class org.janelia.flatfield.FlatfieldCorrection \
            ${stitching_install_path}/target/stitching-spark-1.9.1-SNAPSHOT.jar -i ${destPath}/PetaKit5D_demo_cell_image_dataset/DSR/488nm.json
%}

% Stitching Spark CMD
%{
spark-submit --master spark://10.17.209.11:7077 \
            --driver-cores 2 \
            --driver-memory 2G \
            --executor-memory 500G \
            --executor-cores 24 \
            --class org.janelia.stitching.StitchingSpark \
            ${stitching_install_path}/target/stitching-spark-1.9.1-SNAPSHOT.jar --stitch -i ${destPath}/PetaKit5D_demo_cell_image_dataset/DSR/488nm.json
%}

% Export Spark CMD
%{
spark-submit --master spark://10.17.209.11:7077 \
            --driver-cores 2 \
            --driver-memory 2G \
            --executor-memory 500G \
            --executor-cores 24 \
            --class org.janelia.stitching.StitchingSpark \
            ${stitching_install_path}/target/stitching-spark-1.9.1-SNAPSHOT.jar --fuse -i ${destPath}/PetaKit5D_demo_cell_image_dataset/DSR/488nm-final.json
%}

% N5 to Tiff Spark CMD
%{
spark-submit --master spark://10.17.209.11:7077 \
            --driver-cores 2 \
            --driver-memory 2G \
            --executor-memory 500G \
            --executor-cores 24 \
            --class org.janelia.stitching.N5ToSliceTiffSpark \
            ${stitching_install_path}/target/stitching-spark-1.9.1-SNAPSHOT.jar -i ${destPath}/PetaKit5D_demo_cell_image_dataset/DSR/export.n5
%}


%% Step B9: convert n5 to zarr (optional)
% if you would like to convert the stitched result to zarr, run the code below

% add   ,"n5":"2.0.0"   in the attributes.json before conversion
fn = [dataPath, 'DSR/export.n5/c0/s0'];
flipEmptyValue = true;
masterCompute = true;
N5ToZarr(fn, flipEmptyValue=flipEmptyValue, masterCompute=masterCompute);

% output zarr is in 
% {destPath}/PetaKit5D_demo_cell_image_dataset/DSR/export.n5/c0/zarr/s0.zarr


