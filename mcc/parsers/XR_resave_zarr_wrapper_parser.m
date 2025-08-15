function [] = XR_resave_zarr_wrapper_parser(dataPaths, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) iscell(x) || ischar(x));
ip.addParameter('resultDirName', 'resaved_zarr', @ischar);
ip.addParameter('batchSize', [512, 512, 512], @(x) isnumeric(x) || ischar(x)); % size to process in one batch
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x));
ip.addParameter('shardSize', [], @(x) isnumeric(x) || ischar(x)); % reserved for future use
ip.addParameter('channelPatterns', {'.zarr'}, @(x) iscell(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('largeFile', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
shardSize = pr.shardSize;
channelPatterns = pr.channelPatterns;
inputBbox = pr.inputBbox;
parseCluster = pr.parseCluster;
largeFile = pr.largeFile;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(shardSize)
    shardSize = str2num(shardSize);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(largeFile)
    largeFile = str2num(largeFile);
end
if ischar(masterCompute)
    masterCompute = str2num(masterCompute);
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_resave_zarr_wrapper(dataPaths, resultDirName=resultDirName, batchSize=batchSize, ...
    blockSize=blockSize, shardSize=shardSize, channelPatterns=channelPatterns, ...
    inputBbox=inputBbox, parseCluster=parseCluster, largeFile=largeFile, masterCompute=masterCompute, ...
    jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, uuid=uuid, mccMode=mccMode, ...
    configFile=configFile);

end

