function [] = XR_resave_zarr_wrapper(dataPaths, varargin)
% Dataset level wrapper for resaving zarr files with different parameters.
%
%
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%
% Parameters (as 'specifier'-value pairs):
%       tiffFullpaths : empty or char or cell array (default: ''). A single Tiff file path in char, or a list of Tiff file paths in a cell array. If empty, get the Tiff files from dataPaths.
%       resultDirName : char (default: 'resaved_zarr'). Result directory under data paths.
%           batchSize : 1x3 vector (default: [512, 512, 512]). Batch size per task.
%           blockSize : 1x3 vector (default: [256, 256, 256]). Block/chunk size for zarr output.
%           shardSize : empty or 1x3 vector (default: []). If not empty, set shard size in zarr output. (reserved for future use)
%     channelPatterns : a cell array (default: {'.zarr'}).  Channel identifiers for included channels. 
%           inputBbox : empty or 1x6 vector (default: []). Input bounding box for crop. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%           largeFile : true|false (default: true). Large file processing.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%           jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%         cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%         maxTrialNum : a number (default: 3). The max number of retries for a task.
%        unitWaitTime : a number (default: 1). The wait time per file in minutes to check whether the computing is done.
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%
%
% Author: Xiongtao Ruan (07/31/2025)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) iscell(x) || ischar(x));
ip.addParameter('resultDirName', 'resaved_zarr', @ischar);
ip.addParameter('batchSize', [512, 512, 512], @isnumeric); % size to process in one batch
ip.addParameter('blockSize', [256, 256, 256], @isnumeric);
ip.addParameter('shardSize', [], @isnumeric); % reserved for future use
ip.addParameter('channelPatterns', {'.zarr'}, @iscell);
ip.addParameter('inputBbox', [], @isnumeric);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('largeFile', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 3, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
shardSize= pr.shardSize;
channelPatterns = pr.channelPatterns;
inputBbox = pr.inputBbox;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
largeFile = pr.largeFile;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;


% suppress directory exists warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), filesep)
        dataPaths{d} = [dataPath, filesep];
    end

    dataPaths{d} = [simplifyPath(dataPaths{d}), '/'];
    resultPath = [dataPaths{d}, '/' resultDirName, '/'];
    resultPaths{d} = resultPath;
    mkdir(resultPath);
    fileattrib(resultPath, '+w', 'g');

    save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
    writeJsonFile(pr, [resultPath, '/parameters.json']);
end

if isempty(uuid)
    uuid = get_uuid();
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath, jobLogDir);
end

% parse image filenames
zarrFile = true;
[fnames, fsns, fd_inds, filepaths] = parseImageFilenames(dataPaths, zarrFile, channelPatterns);
nF = numel(fnames);
sz = getImageSize(filepaths{1});
dtype = getImageDataType(filepaths{1});
byteNum = dataTypeToByteNumber(dtype);

if ~largeFile
    % if there are too many batches, use large scale strategy, otherwise, set batchSize as image size
    if isempty(inputBbox)
        outSize = sz;
    else
        outSize = inputBbox(4 : 6) - inputBbox(1 : 3) + 1;
    end
    if prod(ceil(outSize ./ max(512, max(batchSize, blockSize)))) > 10
        largeFile = true;
    else
        batchSize = outSize;
    end
end

%% use generic framework for the computing

frameFullpaths = filepaths;

ext = '.zarr';
resultFullpaths = arrayfun(@(x) sprintf('%s/%s%s', resultPaths{fd_inds(x)}, fsns{x}, ext), ...
    1 : nF, 'UniformOutput', false);

func_strs = arrayfun(@(x) sprintf(['XR_resaveSingleZarr(''%s'',''%s'',''inputBbox'',%s,', ...
    '''blockSize'',%s,''batchSize'',%s,''shardSize'',%s,''parseCluster'',%s,', ...
    '''masterCompute'',%s,''cpusPerTask'',%d,''uuid'',''%s'',''mccMode'',%s,''configFile'',''%s'')'], ...
    frameFullpaths{x}, resultFullpaths{x}, mat2str_comma(inputBbox), mat2str_comma(blockSize), ...
    mat2str_comma(batchSize), mat2str_comma(shardSize), string(parseCluster), ...
    string(masterCompute), cpusPerTask, uuid, string(mccMode), string(configFile)), ...
    1 : nF, 'unif', false);

if largeFile
    memAllocate = prod(batchSize) * byteNum / 2^30 * 1.5;
else
    memAllocate = prod(sz) * byteNum / 2^30 * 1.5;    
end

generic_computing_frameworks_wrapper(frameFullpaths, resultFullpaths, func_strs, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, memAllocate=memAllocate, mccMode=mccMode, configFile=configFile);

end

