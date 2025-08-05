function [] = XR_resaveSingleZarr_parser(zarrFullpath, resultFullpath, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath', @ischar); 
ip.addRequired('resultFullpath', @ischar); 
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x)); % inputBbox for input
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x)); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @(x) isnumeric(x) || ischar(x)); % size to process in one batch 
ip.addParameter('shardSize', [512, 512, 512], @(x) isnumeric(x) || ischar(x)); % size in one shard 
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x));
ip.addParameter('cpusPerTask', 1, @(x) isscalar(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(zarrFullpath, resultFullpath, varargin{:});

pr = ip.Results;
inputBbox = pr.inputBbox;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
shardSize = pr.shardSize;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(shardSize)
    shardSize = str2num(shardSize);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
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

XR_resaveSingleZarr(zarrFullpath, resultFullpath, inputBbox=inputBbox, blockSize=blockSize, ...
    batchSize=batchSize, shardSize=shardSize, parseCluster=parseCluster, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, uuid=uuid, mccMode=mccMode, configFile=configFile);

end

