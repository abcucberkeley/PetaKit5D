function XR_unmix_channels_zarr_parser(zarrFullpaths, unmixFactors, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('unmixFactors', @(x) isnumeric(x) || ischar(x));
ip.addParameter('mode', 'linear', @ischar); % linear vs gaussian
ip.addParameter('unmixSigmas', [], @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('resultDirName', 'Unmixed', @ischar); 
ip.addParameter('channelInd', 1, @(x) isnumeric(x) || ischar(x)); % unmix for which channel
ip.addParameter('batchSize', [1024, 1024, 1024] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('borderSize', [0, 0, 0] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(zarrFullpaths, unmixFactors, varargin{:});

pr = ip.Results;
mode = pr.mode;
unmixSigmas = pr.unmixSigmas;
resultDirName = pr.resultDirName;
channelInd = pr.channelInd;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
borderSize = pr.borderSize;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
configFile = pr.configFile;
mccMode = pr.mccMode;
uuid = pr.uuid;
debug = pr.debug;

if ischar(zarrFullpaths) && ~isempty(zarrFullpaths) && strcmp(zarrFullpaths(1), '{')
    zarrFullpaths = eval(zarrFullpaths);
end
if ischar(unmixFactors)
    unmixFactors = str2num(unmixFactors);
end
if ischar(unmixSigmas)
    unmixSigmas = str2num(unmixSigmas);
end
if ischar(channelInd)
    channelInd = str2num(channelInd);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(borderSize)
    borderSize = str2num(borderSize);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(parseParfor)
    parseParfor = str2num(parseParfor);
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
if ischar(debug)
    debug = str2num(debug);
end

XR_unmix_channels_zarr(zarrFullpaths, unmixFactors, mode=mode, unmixSigmas=unmixSigmas, ...
    resultDirName=resultDirName, channelInd=channelInd, batchSize=batchSize, ...
    blockSize=blockSize, borderSize=borderSize, parseCluster=parseCluster, ...
    parseParfor=parseParfor, masterCompute=masterCompute, jobLogDir=jobLogDir, ...
    cpusPerTask=cpusPerTask, configFile=configFile, mccMode=mccMode, uuid=uuid, ...
    debug=debug);

end

