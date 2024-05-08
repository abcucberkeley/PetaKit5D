function [] = XR_psf_analysis_wrapper_parser(dataPaths, varargin)


%#function XR_psf_analysis_plot

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dz', 0.1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('skewAngle', 32.45, @(x) isnumeric(x) || ischar(x));
ip.addParameter('deskew', true, @(x) islogical(x) || ischar(x));
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('objectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zStageScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('channels', [488, 560], @(x) isnumeric(x) || ischar(x));
ip.addParameter('save16bit', true, @(x) islogical(x) || ischar(x));
ip.addParameter('bgFactor', 1.5, @(x) isnumeric(x) || ischar(x));
ip.addParameter('RWFn', {'/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif', '/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_605em_128_128_101_100nmSteps.tif'}, @(x) iscell(x) || ischar(x));
ip.addParameter('sourceStr', 'test', @ischar);
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', false, @(x) islogical(x) || ischar(x));
ip.addParameter('cpusPerTask', 8, @(x) isscalar(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
skewAngle = pr.skewAngle;
deskew = pr.deskew;
flipZstack = pr.flipZstack;
objectiveScan = pr.objectiveScan;
zStageScan = pr.zStageScan;
channelPatterns = pr.channelPatterns;
channels = pr.channels;
save16bit = pr.save16bit;
bgFactor = pr.bgFactor;
RWFn = pr.RWFn;
sourceStr = pr.sourceStr;
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
if ischar(skewAngle)
    skewAngle = str2num(skewAngle);
end
if ischar(deskew)
    deskew = str2num(deskew);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(objectiveScan)
    objectiveScan = str2num(objectiveScan);
end
if ischar(zStageScan)
    zStageScan = str2num(zStageScan);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(channels)
    channels = str2num(channels);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(bgFactor)
    bgFactor = str2num(bgFactor);
end
if ischar(RWFn) && ~isempty(RWFn) && strcmp(RWFn(1), '{')
    RWFn = eval(RWFn);
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

XR_psf_analysis_wrapper(dataPaths, xyPixelSize=xyPixelSize, dz=dz, skewAngle=skewAngle, ...
    deskew=deskew, flipZstack=flipZstack, objectiveScan=objectiveScan, zStageScan=zStageScan, ...
    channelPatterns=channelPatterns, channels=channels, save16bit=save16bit, ...
    bgFactor=bgFactor, RWFn=RWFn, sourceStr=sourceStr, parseCluster=parseCluster, ...
    masterCompute=masterCompute, cpusPerTask=cpusPerTask, mccMode=mccMode, ...
    configFile=configFile);

end

