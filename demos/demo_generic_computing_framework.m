% demo to run a user defined function using the generic computing framework
% in a distributed way.
% 
% Note: Distributed computing with multiple computing nodes only works on a 
% Linux cluster with slurm. It works on local machines on Windows, Linux and MacOS. 
%
% The parameters demonstrated here are usually a subset of those available 
% for the functions, with the rest using default values. For a comprehensive 
% list of parameters and their defaults, please see the function's parameter 
% list (or input parser) or refer to the parameter documentation (major_functions_documentation.txt).


clear, clc;

fprintf('Generic computing framework demo...\n\n');

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
% download the example dataset from zenodo (https://doi.org/10.5281/zenodo.1492027) manually, 
% or use the code below to download the data from Dropbox
if ispc
    destPath = fullfile(getenv('USERPROFILE'), 'Downloads');
    destPath = strrep(destPath, '\', '/');    
else
    destPath = '~/Downloads/';
end
demo_data_downloader(destPath);

dataPath = [destPath, '/PetaKit5D_demo_cell_image_dataset/'];


%% example define configuration files for the slurm cluster (skip this step if running locally)
% The framework needs a configuration file for the slurm cluster job submission
% for the whole list of parameters that can be configured, please refer to
% the parameters in generic_computing_frameworks_wrapper.m
% Below is an example for the configuration file, and you may include a
% subset of parameters based on your slurm cluster. 

% if using MATLAB runtime, setup the runtime path in your cluster. 
MCRParam = '/global/home/users/user/bin/MATLAB_Runtime/R2023a';

% the mcc master script path in your cluster. It is within PetaKit5D/mcc/linux/run_mccMaster.sh
MCCMasterStr = '/clusterfs/fiona/user/PetaKit5D/mcc/linux/run_mccMaster.sh';

% bash commands to be executed before running for worker jobs, i.e., setup environmental variables. 
% here we show an example for loading the GNU-parallel function for parallel running within a job. 
BashLaunchStr = 'module load parallel/20190822';

% slurm parameter for node allocation, just include partition, qos, and mem-per-cpus
% also each worker job will only run within a single node (or several cores within a node)
SlurmParam = '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M';

% path of the directory to store job logs to avoid jobs logs in the code repo
jobLogDir = '/clusterfs/fiona/user/Projects/job_logs';

% max number of cpu cores, equal or less than the number of cores in a node
maxCPUNum = 24;

% minimum number of cpu cores, minimum 1. 
minCPUNum = 1;

% allocate a whole node or just several cores for job running. If
% allocating a whole node, usually use GNU parallel to run several tasks in
% parallel within a node. This is mostly used in clusters that disallow
% allocating several cores (partial node). 
wholeNodeJob = true;

% memory mapping per cpu core in GB. If there are 502 GB in 24 core node,
% then it is 20.9 GB per core. 
MemPerCPU = 20.9;

% use GNU parallel for computing, only available in mccMode 
GNUparallel = true;

% job time limit in hr, i.e. 12 hrs
jobTimeLimit = 12;

% max number of jobs submitted at once. 
maxJobNum = 5000;

% format as a struct
a = struct();
a.MCRParam = MCRParam;
a.MCCMasterStr = MCCMasterStr;
a.BashLaunchStr = BashLaunchStr;
a.SlurmParam = SlurmParam;
a.jobLogDir = jobLogDir;
a.maxCPUNum = maxCPUNum;
a.minCPUNum = minCPUNum;
a.wholeNodeJob = wholeNodeJob;
a.MemPerCPU = MemPerCPU;
a.GNUparallel = GNUparallel;
a.jobTimeLimit = jobTimeLimit;
a.maxJobNum = maxJobNum;

% save as json file
s = jsonencode(a, PrettyPrint=true);
configFile = [pwd, '/demos/test_slurm_cpu_config_file.json'];
fid = fopen(configFile, 'w');
fprintf(fid, s);
fclose(fid);


%% example to define a user-defined function and use the generic computing framework
% 
% Here we define the crop and deskew/rotate in this user defined function
% Here we put the function in this session as an illustration, it should be saved as
% an independent function (crop_deskew_rotate_crop_demo_function.m), so worker jobs can call it. 
%

fprintf(['The code below is just for the illustration of how to set up a function ', ...
        'for the generic computing framework. \nThe actual function is crop_deskew_rotate_crop_demo_function.m\n'])

% disable the code 
if false

% function [] = crop_deskew_rotate_demo_function(inputFullpath, outputFullpath, skewAngle, xyPixelSize, dz, inBbox)
% the function needs the path for the input and output, along with
% parameters. The function needs to first check if input and output exist,
% it should skip the output file if it exists.

if ~exist(inputFullpath, 'file')
    error('The input file %s does not exist!', inputFullpath);
end

if exist(outputFullpath, 'file')
    fprintf('The output file %s already exist, skip it!\n', outputFullpath);
end

% read input file
im = readtiff(inputFullpath);

% crop the image
im = crop3d_mex(im, inBbox);

% deskew/rotate for the cropped image
% these are default parameters for the demo images
Reverse = true;
ObjectiveScan = false;
resample = [];
Interp = 'linear';

dsr = deskewRotateFrame3D(single(im), Angle, dz, xyPixelSize, 'reverse', Reverse, ...
    'Crop', true, 'ObjectiveScan', ObjectiveScan, 'resample', resample, 'Interp', Interp);

% convert image to uint16
dsr = uint16(dsr);

% write the image to disk with the given file name
% also to ensure the result file is complete, we first write to an
% temporary file, and then rename the temporary file to the final output file

% get uuid
uuid = get_uuid();
tmpFn = sprintf('%s_%s.tif', outputFullpath(1 : end - 4));

% write to temporary file
writetiff(dsr, tmpFn);

% rename to the final output
movefile(tmpFn, outputFullpath);

% end

end

%% setup running for multiple images

% result folder:
% {destPath}/PetaKit5D_demo_cell_image_dataset/Cropped_DSR/

% get filenames
dir_info = dir([dataPath, '*.tif']);
fsns = {dir_info.name};
disp(fsns);

% concatenate folder and file names
inputFullnames = cellfun(@(x) [dataPath, x], fsns, 'UniformOutput', false);

% define output filename
outPath = [dataPath, 'Cropped_DSR/'];
mkdir(outPath);

outputFullnames = cellfun(@(x) [outPath, x], fsns, 'UniformOutput', false);

% define parameters

% skew angle
skewAngle = 32.45;
% xy pixel size
xyPixelSize = 0.108;
% z scan step size
dz = 0.3;
% input bounding box for the crop: [ymin, xmin, zmin, ymax, xmax, zmax]
inBbox = [101, 1, 1, 1700, 512, 500];

% define function strings for the processing
nF = numel(inputFullnames);
funcStrs = cell(nF, 1);

for f = 1 : nF
    funcStrs{f} = sprintf(['crop_deskew_rotate_demo_function(''%s'',''%s'',%d,', ...
                    '%d,%d,%s)'], inputFullnames{f}, outputFullnames{f}, skewAngle, ...
                    xyPixelSize, dz, strrep(mat2str(inBbox), ' ', ','));
end

% setup parameters and run the function using the generic computing framework

% use slurm cluster if true, otherwise use the local machine (master job)
parseCluster = false;

% use master job for task computing or not. 
masterCompute = true;

% estimated memory requirement in GB
memAllocate = 20;

% configuration file for job submission
configFile = '';

% if true, use Matlab runtime (for the situation without matlab license),
% in this case, you can create a parser for the function, and add it to
% mccMaster.m and compile the code with compile_mccMaster.m
% the parser function for crop_deskew_rotate_demo_function.m is included in
% the demo directory: crop_deskew_rotate_demo_function_parser.m
mccMode = false;

generic_computing_frameworks_wrapper(inputFullnames, outputFullnames, funcStrs, ...
    parseCluster=parseCluster, masterCompute=masterCompute, memAllocate=memAllocate, ...
    configFile=configFile, mccMode=mccMode);

% note: if there is no slurm cluster, it will run the task sequentially in
% the current matlab session. In this situation, please delete any .tmp
% file in the result folder from old incomplete runs if there are any.
% Otherwise, the function will assume other workers are working on these tasks. 





