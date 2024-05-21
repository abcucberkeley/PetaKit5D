function XR_FSC_analysis_wrapper_parser(dataPaths, varargin)


%#function XR_one_image_FSC_analysis_frame

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'FSCs', @ischar);
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dz', 0.1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dr', 10, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dtheta', pi / 12 , @(x) isnumeric(x) || ischar(x));
ip.addParameter('resThreshMethod', 'fixed', @ischar);
ip.addParameter('resThresh', 0.2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('halfSize', [251, 251, 251], @(x) isnumeric(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('resAxis', 'xz', @ischar);
ip.addParameter('skipConeRegion', true, @(x) islogical(x) || ischar(x));
ip.addParameter('channelPatterns', {'tif'}, @(x) iscell(x) || ischar(x));
ip.addParameter('Channels', [488, 560], @(x) isnumeric(x) || ischar(x));
ip.addParameter('multiRegions', false, @(x) islogical(x) || ischar(x));
ip.addParameter('regionInterval', [50, 50, -1], @(x) isnumeric(x) || ischar(x)); % yxz, -1 means only center
ip.addParameter('regionGrid', [], @(x) isnumeric(x) || ischar(x)); % user provided grid for region centers, N x 3
ip.addParameter('clipPer', [], @(x) isnumeric(x) || ischar(x)); % clip intensity higher than the given percentile
ip.addParameter('suffix', 'decon', @ischar);
ip.addParameter('iterInterval', 5, @(x) isnumeric(x) || ischar(x)); % iteration interval for FSC resolution plot
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x));
ip.addParameter('cpusPerTask', 4, @(x) isscalar(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
dr = pr.dr;
dtheta = pr.dtheta;
resThreshMethod = pr.resThreshMethod;
resThresh = pr.resThresh;
halfSize = pr.halfSize;
inputBbox = pr.inputBbox;
resAxis = pr.resAxis;
skipConeRegion = pr.skipConeRegion;
channelPatterns = pr.channelPatterns;
Channels = pr.Channels;
multiRegions = pr.multiRegions;
regionInterval = pr.regionInterval;
regionGrid = pr.regionGrid;
clipPer = pr.clipPer;
suffix = pr.suffix;
iterInterval = pr.iterInterval;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(dr)
    dr = str2num(dr);
end
if ischar(dtheta)
    dtheta = str2num(dtheta);
end
if ischar(resThresh)
    resThresh = str2num(resThresh);
end
if ischar(halfSize)
    halfSize = str2num(halfSize);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(skipConeRegion)
    skipConeRegion = str2num(skipConeRegion);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(Channels)
    Channels = str2num(Channels);
end
if ischar(multiRegions)
    multiRegions = str2num(multiRegions);
end
if ischar(regionInterval)
    regionInterval = str2num(regionInterval);
end
if ischar(regionGrid)
    regionGrid = str2num(regionGrid);
end
if ischar(clipPer)
    clipPer = str2num(clipPer);
end
if ischar(iterInterval)
    iterInterval = str2num(iterInterval);
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

XR_FSC_analysis_wrapper(dataPaths, resultDirName=resultDirName, xyPixelSize=xyPixelSize, ...
    dz=dz, dr=dr, dtheta=dtheta, resThreshMethod=resThreshMethod, resThresh=resThresh, ...
    halfSize=halfSize, inputBbox=inputBbox, resAxis=resAxis, skipConeRegion=skipConeRegion, ...
    channelPatterns=channelPatterns, Channels=Channels, multiRegions=multiRegions, ...
    regionInterval=regionInterval, regionGrid=regionGrid, clipPer=clipPer, ...
    suffix=suffix, iterInterval=iterInterval, parseCluster=parseCluster, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, mccMode=mccMode, configFile=configFile);

end

