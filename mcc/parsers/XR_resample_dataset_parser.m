function [] = XR_resample_dataset_parser(dataPaths, resampleFactor, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('resampleFactor', @(x) isnumeric(x) || ischar(x));
ip.addParameter('resultDirName', 'resampled', @ischar);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @(x) iscell(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x)); % bbox for input
ip.addParameter('interpMethod', 'linear', @(x) ischar(x) && any(strcmpi(x, {'cubic', 'linear', 'nearest', 'max', 'mean'})));
ip.addParameter('save16bit', true, @(x) islogical(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); % use zarr file as output
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x)); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @(x) isnumeric(x) || ischar(x)); % size to process in one batch
ip.addParameter('borderSize', [5, 5, 5], @(x) isnumeric(x) || ischar(x)); % padded boarder for each batch
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 2, @(x) isscalar(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, resampleFactor, varargin{:});

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

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
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
if ischar(masterCompute)
    masterCompute = str2num(masterCompute);
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_resample_dataset(dataPaths, resampleFactor, resultDirName=resultDirName, ...
    channelPatterns=channelPatterns, inputBbox=inputBbox, interpMethod=interpMethod, ...
    save16bit=save16bit, zarrFile=zarrFile, largeFile=largeFile, saveZarr=saveZarr, ...
    blockSize=blockSize, batchSize=batchSize, borderSize=borderSize, parseCluster=parseCluster, ...
    jobLogDir=jobLogDir, masterCompute=masterCompute, cpusPerTask=cpusPerTask, ...
    uuid=uuid, mccMode=mccMode, configFile=configFile);

end

