function [] = XR_unmix_channels_data_wrapper_parser(dataPaths, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('unmixFactors', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('mode', 'linear', @ischar); % linear vs gaussian
ip.addParameter('unmixSigmas', [], @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('resultDirName', 'Unmixed', @ischar); 
ip.addParameter('channelPatterns', {'CamA', 'CamB'}, @(x) iscell(x) || ischar(x));
ip.addParameter('channelInd', 1, @(x) isnumeric(x) || ischar(x)); % unmix for which channel
ip.addParameter('FFCorrection', false, @(x) islogical(x) || ischar(x));
ip.addParameter('lowerLimit', 0.4, @(x) isnumeric(x) || ischar(x));
ip.addParameter('FFImagePaths', {'',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('backgroundPaths', {'',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('constBackground', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('constOffset', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x));
ip.addParameter('save16bit', true, @(x) islogical(x) || ischar(x));
ip.addParameter('batchSize', [1024, 1024, 1024] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('borderSize', [0, 0, 0] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
unmixFactors = pr.unmixFactors;
mode = pr.mode;
unmixSigmas = pr.unmixSigmas;
resultDirName = pr.resultDirName;
channelPatterns = pr.channelPatterns;
channelInd = pr.channelInd;
FFCorrection = pr.FFCorrection;
lowerLimit = pr.lowerLimit;
FFImagePaths = pr.FFImagePaths;
backgroundPaths = pr.backgroundPaths;
constBackground = pr.constBackground;
constOffset = pr.constOffset;
zarrFile = pr.zarrFile;
largeFile = pr.largeFile;
saveZarr = pr.saveZarr;
save16bit = pr.save16bit;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
borderSize = pr.borderSize;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
configFile = pr.configFile;
mccMode = pr.mccMode;
uuid = pr.uuid;
debug = pr.debug;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(unmixFactors)
    unmixFactors = str2num(unmixFactors);
end
if ischar(unmixSigmas)
    unmixSigmas = str2num(unmixSigmas);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(channelInd)
    channelInd = str2num(channelInd);
end
if ischar(FFCorrection)
    FFCorrection = str2num(FFCorrection);
end
if ischar(lowerLimit)
    lowerLimit = str2num(lowerLimit);
end
if ischar(FFImagePaths) && ~isempty(FFImagePaths) && strcmp(FFImagePaths(1), '{')
    FFImagePaths = eval(FFImagePaths);
end
if ischar(backgroundPaths) && ~isempty(backgroundPaths) && strcmp(backgroundPaths(1), '{')
    backgroundPaths = eval(backgroundPaths);
end
if ischar(constBackground)
    constBackground = str2num(constBackground);
end
if ischar(constOffset)
    constOffset = str2num(constOffset);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(largeFile)
    largeFile = str2num(largeFile);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
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

XR_unmix_channels_data_wrapper(dataPaths, unmixFactors=unmixFactors, mode=mode, ...
    unmixSigmas=unmixSigmas, resultDirName=resultDirName, channelPatterns=channelPatterns, ...
    channelInd=channelInd, FFCorrection=FFCorrection, lowerLimit=lowerLimit, ...
    FFImagePaths=FFImagePaths, backgroundPaths=backgroundPaths, constBackground=constBackground, ...
    constOffset=constOffset, zarrFile=zarrFile, largeFile=largeFile, saveZarr=saveZarr, ...
    save16bit=save16bit, batchSize=batchSize, blockSize=blockSize, borderSize=borderSize, ...
    parseCluster=parseCluster, masterCompute=masterCompute, jobLogDir=jobLogDir, ...
    cpusPerTask=cpusPerTask, configFile=configFile, mccMode=mccMode, uuid=uuid, ...
    debug=debug);

end

