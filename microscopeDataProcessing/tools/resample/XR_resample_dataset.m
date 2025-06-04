function [] = XR_resample_dataset(dataPaths, resampleFactor, varargin)
% Dataset level wrapper for resampling tools. 
% ResampleFactor [5,5,5] (yxz) means downsample data by 5x5x5.
% 
%
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%      resampleFactor : 1x1, 1x2 or 1x3 vector. Resampling factor in yxz order. Greater than 1 means downsampling. The first number for y, and the last number for z, the middle (if exist) or the first number for x.
%
% Parameters (as 'specifier'-value pairs):
%       resultDirName : char (default: 'matlab_stitch'). Resampling result directory under data path.
%     channelPatterns : a cell array (default: {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}).  Channel identifiers for included channels. 
%           inputBbox : empty or 1x6 vector (default: []). Input bounding box for crop. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%        interpMethod : 'linear'|'cubic'|'nearest'|'max'|'mean' (default: 'linear'). Interpolation method. 'max' and 'mean' only support integer resample factors.
%           save16bit : true|false (default: true). Save 16bit result for the result. 
%            zarrFile : true|false (default: false). Use Zarr file as input.
%           largeFile : true|false (default: false). Use large scale resampling strategy with batch processing. Only for Zarr files.
%            saveZarr : true|false (default: false). Save results as Zarr files.
%           blockSize : 1x3 vector (default: [256, 256, 256]). Block/chunk size for zarr output.
%           batchSize : 1x3 vector (default: [512, 512, 512]). Batch size per stitching task.
%          borderSize : 1x3 vector (default: [5, 5, 5]. Padded border for each batch.
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%           jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%         cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
% 
%
% Author: Xiongtao Ruan (01/15/2021)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('resampleFactor', @isnumeric);
ip.addParameter('resultDirName', 'resampled', @ischar);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @iscell);
ip.addParameter('inputBbox', [], @isnumeric); % bbox for input
ip.addParameter('interpMethod', 'linear', @(x) ischar(x) && any(strcmpi(x, {'cubic', 'linear', 'nearest', 'max', 'mean'})));
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('saveZarr', false, @islogical); % use zarr file as output
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @isnumeric); % size to process in one batch
ip.addParameter('borderSize', [5, 5, 5], @isnumeric); % padded boarder for each batch
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 2, @isscalar);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, resampleFactor, varargin{:});

warning('off', 'MATLAB:MKDIR:DirectoryExists');

pr = ip.Results;
resultDirName = pr.resultDirName;
channelPatterns = pr.channelPatterns;
inputBbox = pr.inputBbox;
interpMethod = pr.interpMethod;
save16bit = pr.save16bit;
zarrFile = pr.zarrFile;
largeFile = pr.largeFile;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
borderSize = pr.borderSize;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if isempty(uuid)
    uuid = get_uuid();
end

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPaths{d} = simplifyPath(dataPaths{d});
    resultPath = [dataPaths{d}, '/' resultDirName, '/'];
    resultPaths{d} = resultPath;
    mkdir(resultPath);
    fileattrib(resultPath, '+w', 'g');

    save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
    writeJsonFile(pr, [resultPath, '/parameters.json']);
end

% parse image filenames
[fnames, fsns, fd_inds, filepaths] = parseImageFilenames(dataPaths, zarrFile, channelPatterns);
nF = numel(fnames);
sz = getImageSize(filepaths{1});
dtype = getImageDataType(filepaths{1});
byteNum = dataTypeToByteNumber(dtype);

%% use generic framework for the cropping do computing

frameFullpaths = filepaths;

ext = '.tif';
if saveZarr
    ext = '.zarr';
end
resultFullpaths = arrayfun(@(x) sprintf('%s/%s%s', resultPaths{fd_inds(x)}, fsns{x}, ext), ...
    1 : nF, 'UniformOutput', false);

if largeFile
    func_strs = arrayfun(@(x) sprintf(['XR_resampleSingleZarr(''%s'',''%s'',%s,', ...
        '''inputBbox'',%s,''blockSize'',%s,''batchSize'',%s,''borderSize'',%s,', ...
        '''interpMethod'',''%s'',''parseCluster'',%s,''cpusPerTask'',%d,''uuid'',''%s'',', ...
        '''mccMode'',%s,''configFile'',''%s'')'], frameFullpaths{x}, ...
        resultFullpaths{x}, mat2str_comma(resampleFactor, 10), mat2str_comma(inputBbox), ...
        mat2str_comma(blockSize), mat2str_comma(batchSize), mat2str_comma(borderSize), ...
        interpMethod, string(parseCluster),cpusPerTask, uuid, string(mccMode), string(configFile)), ...
        1 : nF, 'unif', false);

    memAllocate = prod(batchSize) * byteNum / 2^30 * (prod(resampleFactor) + 1) * 3;
else
    func_strs = arrayfun(@(x) sprintf(['XR_resampleFrame(''%s'',''%s'',[%s],', ...
        '''inputBbox'',%s,''zarrFile'',%s,''saveZarr'',%s,''interpMethod'',''%s'',''save16bit'',%s,', ...
        '''uuid'',''%s'')'], frameFullpaths{x}, resultFullpaths{x}, mat2str_comma(resampleFactor, 10), ...
        mat2str_comma(inputBbox), string(zarrFile), string(saveZarr), interpMethod, string(save16bit), ...
        uuid), 1 : nF, 'unif', false);

    memAllocate = prod(sz) * byteNum / 1024^3 * (2 + 2 / prod(resampleFactor));
end

generic_computing_frameworks_wrapper(frameFullpaths, resultFullpaths, func_strs, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, memAllocate=memAllocate, mccMode=mccMode, configFile=configFile);

end
