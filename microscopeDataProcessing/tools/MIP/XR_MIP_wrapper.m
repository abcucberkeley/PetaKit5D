function [] = XR_MIP_wrapper(dataPaths, varargin)   
% dataset level wrapper for MIP generation
%
%
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%
% Parameters (as 'specifier'-value pairs):
%       resultDirName : char (default: 'matlab_decon'). Result directory under data paths.
%                axis : 1x3 vector (default: [0, 0, 1]). MIP generation for given axis in yxz order. It generate the MIP when the number is nonzero. 
%     channelPatterns : a cell array (default: {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}).  Channel identifiers for included channels. 
%            zarrFile : true|false (default: false). Use Zarr file as input.
%           largeFile : true|false (default: false). Use large scale deskew/rotation strategy.
%           batchSize : 1x3 vector (default: [2048, 2048, 2048]). Batch size per stitching task.
%           save16bit : true|false (default: true). Save 16bit result for deskew/rotate. 
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%         parseParfor : true|false (default: false). Use matlab parfor for paralle processing.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%         cpusPerTask : a number (default: 3). The number of cpu cores per task for slurm job submission.
%           jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%               debug : true|false (default: false). Debug mode. Not actually used in this function. Reserved for future use.
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%
%
% Author: Xiongtao Ruan


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'MIPs', @ischar);
ip.addParameter('axis', [0, 0, 1], @isnumeric);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @iscell);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('batchSize', [2048, 2048, 2048] , @isvector);
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); 
ip.addParameter('cpusPerTask', 3, @isscalar);
ip.addParameter('jobLogDir', '../job_logs/', @ischar);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
axis = pr.axis;
channelPatterns =  pr.channelPatterns;
zarrFile = pr.zarrFile;
largeFile = pr.largeFile;
batchSize = pr.batchSize;
save16bit = pr.save16bit;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
configFile = pr.configFile;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

if save16bit
    dtype = 'uint16';
else
    dtype = 'single';
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    resultPath = sprintf('%s/%s/', dataPath, resultDirName);
    mkdir(resultPath);
    
    dataPaths{d} = simplifyPath(dataPath);
    resultPaths{d} = simplifyPath(resultPath);
    try 
        fileattrib(resultPath, '+w', 'g');
    catch ME
        disp(ME);
    end

    save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
    writeJsonFile(pr, [resultPath, '/parameters.json']);
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPaths{1}, jobLogDir);
end

% parse image filenames
[~, fsns, fd_inds, filepaths] = parseImageFilenames(dataPaths, zarrFile, channelPatterns);
nF = numel(fsns);
sz = getImageSize(filepaths{1});

%% use generic framework for the MIP computing

frameFullpaths = filepaths;
MIPFullpaths = cell(nF, 1);
func_strs = cell(nF, 1);

% use the last select axis as flag if more than one axis for MIP
axis_strs = {'y', 'x', 'z'};
fidx = find(axis, 1, 'last');
axis_str = axis_strs{fidx};

for f = 1 : nF
    frameFullpath = frameFullpaths{f};

    d = fd_inds(f);
    resultPath = resultPaths{d};
    MIPFullpath = sprintf('%s/%s_MIP_%s.tif', resultPath, fsns{f}, axis_str);
    MIPFullpaths{f} = MIPFullpath;
    
    if zarrFile
        if largeFile
            func_strs{f} = sprintf(['XR_MIP_zarr(''%s'',''resultDirName'',''%s'',', ...
                '''axis'',%s,''batchSize'',%s,''parseCluster'',%s,''parseParfor'',%s,', ...
                '''jobLogDir'',''%s'',''mccMode'',%s,''configFile'',''%s'')'], ...
                frameFullpath, resultDirName, mat2str_comma(axis), mat2str_comma(batchSize), ...
                string(parseCluster), string(parseParfor), jobLogDir, string(mccMode), configFile);
        else
            func_strs{f} = sprintf(['saveMIP_zarr(''%s'',''%s'',''%s'',%s)'], frameFullpath, ...
                MIPFullpath, dtype, mat2str_comma(axis));
        end
    else
        func_strs{f} = sprintf(['saveMIP_tiff(''%s'',''%s'',''dtype'',''%s'',''axis'',%s)'], ...
            frameFullpath, MIPFullpath, dtype, mat2str_comma(axis));
    end
end

if parseParfor && ~largeFile
    matlab_parfor_generic_computing_wrapper(frameFullpaths, MIPFullpaths, func_strs, 'GPUJob', false, 'nworker', 12)
end
taskBatchNum = max(1, round(nF /1000));
memAllocate = prod(sz) * 4 / 1024^3 * 4;
generic_computing_frameworks_wrapper(frameFullpaths, MIPFullpaths, func_strs, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    taskBatchNum=taskBatchNum, cpusPerTask=cpusPerTask, memAllocate=memAllocate, ...
    mccMode=mccMode, configFile=configFile);

end

