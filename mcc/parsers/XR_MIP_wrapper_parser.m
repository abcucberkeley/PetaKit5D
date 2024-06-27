function [] = XR_MIP_wrapper_parser(dataPaths, varargin)   


%#function XR_MIP_zarr
%#function saveMIP_zarr
%#function saveMIP_tiff

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'MIPs', @ischar);
ip.addParameter('axis', [0, 0, 1], @(x) isnumeric(x) || ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @(x) iscell(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('batchSize', [2048, 2048, 2048] , @(x) isvector(x) || ischar(x));
ip.addParameter('save16bit', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); 
ip.addParameter('cpusPerTask', 3, @(x) isscalar(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs/', @ischar);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
axis = pr.axis;
channelPatterns = pr.channelPatterns;
zarrFile = pr.zarrFile;
largeFile = pr.largeFile;
batchSize = pr.batchSize;
save16bit = pr.save16bit;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
jobLogDir = pr.jobLogDir;
uuid = pr.uuid;
debug = pr.debug;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(axis)
    axis = str2num(axis);
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
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
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
if ischar(debug)
    debug = str2num(debug);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_MIP_wrapper(dataPaths, resultDirName=resultDirName, axis=axis, channelPatterns=channelPatterns, ...
    zarrFile=zarrFile, largeFile=largeFile, batchSize=batchSize, save16bit=save16bit, ...
    parseCluster=parseCluster, parseParfor=parseParfor, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, jobLogDir=jobLogDir, uuid=uuid, debug=debug, mccMode=mccMode, ...
    configFile=configFile);

end

