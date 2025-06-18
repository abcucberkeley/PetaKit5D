function [] = XR_resampleSingleZarr_parser(zarrFullpath, rsFullpath, resampleFactor, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath', @ischar); 
ip.addRequired('rsFullpath', @ischar); 
ip.addRequired('resampleFactor', @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x)); % inputBbox for input
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x)); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @(x) isnumeric(x) || ischar(x)); % size to process in one batch 
ip.addParameter('borderSize', [5, 5, 5], @(x) isnumeric(x) || ischar(x)); % padded boarder for each batch
ip.addParameter('interpMethod', 'linear', @(x) ischar(x) && any(strcmpi(x, {'cubic', 'linear', 'nearest', 'max', 'mean'})));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('cpusPerTask', 1, @(x) isscalar(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(zarrFullpath, rsFullpath, resampleFactor, varargin{:});

pr = ip.Results;
inputBbox = pr.inputBbox;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
borderSize = pr.borderSize;
interpMethod = pr.interpMethod;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(borderSize)
    borderSize = str2num(borderSize);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_resampleSingleZarr(zarrFullpath, rsFullpath, resampleFactor, inputBbox=inputBbox, ...
    blockSize=blockSize, batchSize=batchSize, borderSize=borderSize, interpMethod=interpMethod, ...
    parseCluster=parseCluster, cpusPerTask=cpusPerTask, uuid=uuid, mccMode=mccMode, ...
    configFile=configFile);

end

