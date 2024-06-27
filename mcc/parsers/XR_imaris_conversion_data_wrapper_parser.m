function XR_imaris_conversion_data_wrapper_parser(dataPaths, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'imaris',  @(x) ischar(x));
ip.addParameter('overwrite', false,  @(x) islogical(x) || ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('pixelSizes', [0.108, 0.108, 0.108], @(x) isnumeric(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('blockSize', [64, 64, 64], @(x) isnumeric(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('converterPath', '', @ischar);
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 24, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
overwrite = pr.overwrite;
channelPatterns = pr.channelPatterns;
pixelSizes = pr.pixelSizes;
zarrFile = pr.zarrFile;
blockSize = pr.blockSize;
inputBbox = pr.inputBbox;
converterPath = pr.converterPath;
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
if ischar(overwrite)
    overwrite = str2num(overwrite);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(pixelSizes)
    pixelSizes = str2num(pixelSizes);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
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

XR_imaris_conversion_data_wrapper(dataPaths, resultDirName=resultDirName, overwrite=overwrite, ...
    channelPatterns=channelPatterns, pixelSizes=pixelSizes, zarrFile=zarrFile, ...
    blockSize=blockSize, inputBbox=inputBbox, converterPath=converterPath, ...
    parseCluster=parseCluster, masterCompute=masterCompute, jobLogDir=jobLogDir, ...
    cpusPerTask=cpusPerTask, uuid=uuid, mccMode=mccMode, configFile=configFile);

end

