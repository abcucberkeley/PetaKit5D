function [] = XR_tiffToZarr_wrapper_parser(dataPaths, varargin)


%#function tiffToZarr

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) iscell(x) || ischar(x));
ip.addParameter('tiffFullpaths', '', @(x) iscell(x) || ischar(x));
ip.addParameter('resultDirName', 'zarr', @ischar);
ip.addParameter('locIds', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x));
ip.addParameter('shardSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('flippedTile', [], @(x) isempty(x) || islogical(x) || ischar(x));
ip.addParameter('resampleFactor', [], @(x) isempty(x) || isnumeric(x) || ischar(x));
ip.addParameter('partialFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('channelPatterns', {'tif'}, @(x) iscell(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('tileOutBbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('processFunPath', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x) || isstring(x) || iscell(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('bigData', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
tiffFullpaths = pr.tiffFullpaths;
resultDirName = pr.resultDirName;
locIds = pr.locIds;
blockSize = pr.blockSize;
shardSize = pr.shardSize;
flippedTile = pr.flippedTile;
resampleFactor = pr.resampleFactor;
partialFile = pr.partialFile;
channelPatterns = pr.channelPatterns;
inputBbox = pr.inputBbox;
tileOutBbox = pr.tileOutBbox;
processFunPath = pr.processFunPath;
parseCluster = pr.parseCluster;
bigData = pr.bigData;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(tiffFullpaths) && ~isempty(tiffFullpaths) && strcmp(tiffFullpaths(1), '{')
    tiffFullpaths = eval(tiffFullpaths);
end
if ischar(locIds)
    locIds = str2num(locIds);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(shardSize)
    shardSize = str2num(shardSize);
end
if ischar(flippedTile)
    flippedTile = str2num(flippedTile);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(partialFile)
    partialFile = str2num(partialFile);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(tileOutBbox)
    tileOutBbox = str2num(tileOutBbox);
end
if ischar(processFunPath) && ~isempty(processFunPath) && (strcmp(processFunPath(1), '{') || strcmp(processFunPath(1), '[') || strcmp(processFunPath(1), '@'))
    processFunPath = eval(processFunPath);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(bigData)
    bigData = str2num(bigData);
end
if ischar(masterCompute)
    masterCompute = str2num(masterCompute);
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(maxTrialNum)
    maxTrialNum = str2num(maxTrialNum);
end
if ischar(unitWaitTime)
    unitWaitTime = str2num(unitWaitTime);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_tiffToZarr_wrapper(dataPaths, tiffFullpaths=tiffFullpaths, resultDirName=resultDirName, ...
    locIds=locIds, blockSize=blockSize, shardSize=shardSize, flippedTile=flippedTile, ...
    resampleFactor=resampleFactor, partialFile=partialFile, channelPatterns=channelPatterns, ...
    inputBbox=inputBbox, tileOutBbox=tileOutBbox, processFunPath=processFunPath, ...
    parseCluster=parseCluster, bigData=bigData, masterCompute=masterCompute, ...
    jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, uuid=uuid, maxTrialNum=maxTrialNum, ...
    unitWaitTime=unitWaitTime, mccMode=mccMode, configFile=configFile);

end

