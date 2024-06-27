function [] = XR_crop_dataset_parser(dataPaths, inputBbox, varargin)


%#function XR_crop_frame

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('inputBbox', @(x) isnumeric(x) || ischar(x));
ip.addParameter('resultDirName', 'Cropped', @ischar);
ip.addParameter('cropType', 'fixed', @ischar);
ip.addParameter('pad', false, @(x) islogical(x) || ischar(x));
ip.addParameter('lastStartCoords', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('zarrFile', false , @(x) islogical(x) || ischar(x));
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveZarr', false , @(x) islogical(x) || ischar(x));
ip.addParameter('batchSize', [1024, 1024, 1024] , @(x) isnumeric(x) || ischar(x));
ip.addParameter('blockSize', [256, 256, 256] , @(x) isnumeric(x) || ischar(x));
ip.addParameter('save16bit', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, inputBbox, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
cropType = pr.cropType;
pad = pr.pad;
lastStartCoords = pr.lastStartCoords;
channelPatterns = pr.channelPatterns;
zarrFile = pr.zarrFile;
largeFile = pr.largeFile;
saveZarr = pr.saveZarr;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
save16bit = pr.save16bit;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(pad)
    pad = str2num(pad);
end
if ischar(lastStartCoords)
    lastStartCoords = str2num(lastStartCoords);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
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
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
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

XR_crop_dataset(dataPaths, inputBbox, resultDirName=resultDirName, cropType=cropType, ...
    pad=pad, lastStartCoords=lastStartCoords, channelPatterns=channelPatterns, ...
    zarrFile=zarrFile, largeFile=largeFile, saveZarr=saveZarr, batchSize=batchSize, ...
    blockSize=blockSize, save16bit=save16bit, parseCluster=parseCluster, masterCompute=masterCompute, ...
    jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, uuid=uuid, mccMode=mccMode, ...
    configFile=configFile);

end

