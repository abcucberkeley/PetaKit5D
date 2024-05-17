function [] = XR_psf_detection_and_analysis_wrapper_parser(dataPaths, varargin)


%#function XR_psf_detection_and_cropping
%#function XR_psf_analysis_plot

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dz', 0.1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('skewAngle', 32.45, @(x) isnumeric(x) || ischar(x));
ip.addParameter('cropSize', [256, 128, 201], @(x) isnumeric(x) || ischar(x));
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('distThresh', [256, 128, 201], @(x) isnumeric(x) || ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('channels', [488, 560], @(x) isnumeric(x) || ischar(x));
ip.addParameter('RWFn', {'/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif', '/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_605em_128_128_101_100nmSteps.tif'}, @(x) iscell(x) || ischar(x));
ip.addParameter('sourceStr', 'test', @ischar);
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
skewAngle = pr.skewAngle;
cropSize = pr.cropSize;
flipZstack = pr.flipZstack;
distThresh = pr.distThresh;
channelPatterns = pr.channelPatterns;
channels = pr.channels;
RWFn = pr.RWFn;
sourceStr = pr.sourceStr;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
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
if ischar(cropSize)
    cropSize = str2num(cropSize);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(distThresh)
    distThresh = str2num(distThresh);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(channels)
    channels = str2num(channels);
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
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_psf_detection_and_analysis_wrapper(dataPaths, xyPixelSize=xyPixelSize, dz=dz, ...
    skewAngle=skewAngle, cropSize=cropSize, flipZstack=flipZstack, distThresh=distThresh, ...
    channelPatterns=channelPatterns, channels=channels, RWFn=RWFn, sourceStr=sourceStr, ...
    parseCluster=parseCluster, masterCompute=masterCompute, mccMode=mccMode, ...
    configFile=configFile);

end

