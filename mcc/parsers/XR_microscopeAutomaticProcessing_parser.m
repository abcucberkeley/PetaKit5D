function [] = XR_microscopeAutomaticProcessing_parser(dataPaths, varargin)


%#function XR_deskewRotateFrame
%#function XR_matlab_stitching_wrapper
%#function XR_stitching_frame_zarr_dev_v1
%#function XR_cudaDeconFrame3D
%#function XR_cppDeconFrame3D
%#function XR_RLdeconFrame3D
%#function XR_RotateFrame3D

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 3) && islogical(x) || ischar(x));
ip.addParameter('streaming', true,  @(x) islogical(x) || ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('skewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('dz', 0.5, @(x) isscalar(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isscalar(x) || ischar(x));
ip.addParameter('reverse', true, @(x) islogical(x) || ischar(x));
ip.addParameter('objectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zStageScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('save16bit', [true, true], @(x) (numel(x) == 1 || numel(x) == 2) && islogical(x) || ischar(x));
ip.addParameter('dzFromEncoder', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x));
ip.addParameter('save3DStack', true , @(x) islogical(x) || ischar(x));
ip.addParameter('deskew', true, @(x) islogical(x) || ischar(x));
ip.addParameter('rotate', true, @(x) islogical(x) || ischar(x));
ip.addParameter('stitch', false, @(x) islogical(x) || ischar(x));
ip.addParameter('parseSettingFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSRCombined', true, @(x) islogical(x) || ischar(x)); 
ip.addParameter('FFCorrection', false, @(x) islogical(x) || ischar(x));
ip.addParameter('BKRemoval', false, @(x) islogical(x) || ischar(x));
ip.addParameter('lowerLimit', 0.4, @(x) isnumeric(x) || ischar(x));
ip.addParameter('constOffset', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('FFImagePaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('backgroundPaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('resampleType', 'isotropic', @ischar);
ip.addParameter('resampleFactor', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('stitchResultDirName', '', @ischar);
ip.addParameter('imageListFullpaths', '', @(x) ischar(x) || iscell(x));
ip.addParameter('axisOrder', 'xyz', @(x) ischar(x));
ip.addParameter('blendMethod', 'none', @ischar);
ip.addParameter('xcorrShift', false, @(x) islogical(x) || ischar(x));
ip.addParameter('xcorrMode', 'primaryFirst', @(x) ismember(lower(x), {'primary', 'primaryfirst', 'all'}) || ischar(x));
ip.addParameter('xyMaxOffset', 300, @(x) isnumeric(x) || ischar(x));
ip.addParameter('zMaxOffset', 50, @(x) isnumeric(x) || ischar(x));
ip.addParameter('edgeArtifacts', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('primaryCh', '', @ischar);
ip.addParameter('stitchMIP', [], @(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3)) || ischar(x));
ip.addParameter('onlineStitch', false, @(x) islogical(x) || ischar(x));
ip.addParameter('generateImageList', '', @(x) ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('minModifyTime', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('maxModifyTime', 10, @(x) isnumeric(x) || ischar(x));
ip.addParameter('maxWaitLoopNum', 10, @(x) isnumeric(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
overwrite = pr.overwrite;
streaming = pr.streaming;
channelPatterns = pr.channelPatterns;
skewAngle = pr.skewAngle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
reverse = pr.reverse;
objectiveScan = pr.objectiveScan;
zStageScan = pr.zStageScan;
save16bit = pr.save16bit;
dzFromEncoder = pr.dzFromEncoder;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
save3DStack = pr.save3DStack;
deskew = pr.deskew;
rotate = pr.rotate;
stitch = pr.stitch;
parseSettingFile = pr.parseSettingFile;
flipZstack = pr.flipZstack;
DSRCombined = pr.DSRCombined;
FFCorrection = pr.FFCorrection;
BKRemoval = pr.BKRemoval;
lowerLimit = pr.lowerLimit;
constOffset = pr.constOffset;
FFImagePaths = pr.FFImagePaths;
backgroundPaths = pr.backgroundPaths;
resampleType = pr.resampleType;
resampleFactor = pr.resampleFactor;
inputBbox = pr.inputBbox;
stitchResultDirName = pr.stitchResultDirName;
imageListFullpaths = pr.imageListFullpaths;
axisOrder = pr.axisOrder;
blendMethod = pr.blendMethod;
xcorrShift = pr.xcorrShift;
xcorrMode = pr.xcorrMode;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
edgeArtifacts = pr.edgeArtifacts;
primaryCh = pr.primaryCh;
stitchMIP = pr.stitchMIP;
onlineStitch = pr.onlineStitch;
generateImageList = pr.generateImageList;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
minModifyTime = pr.minModifyTime;
maxModifyTime = pr.maxModifyTime;
maxWaitLoopNum = pr.maxWaitLoopNum;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(overwrite)
    overwrite = str2num(overwrite);
end
if ischar(streaming)
    streaming = str2num(streaming);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(skewAngle)
    skewAngle = str2num(skewAngle);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(reverse)
    reverse = str2num(reverse);
end
if ischar(objectiveScan)
    objectiveScan = str2num(objectiveScan);
end
if ischar(zStageScan)
    zStageScan = str2num(zStageScan);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(dzFromEncoder)
    dzFromEncoder = str2num(dzFromEncoder);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(save3DStack)
    save3DStack = str2num(save3DStack);
end
if ischar(deskew)
    deskew = str2num(deskew);
end
if ischar(rotate)
    rotate = str2num(rotate);
end
if ischar(stitch)
    stitch = str2num(stitch);
end
if ischar(parseSettingFile)
    parseSettingFile = str2num(parseSettingFile);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(DSRCombined)
    DSRCombined = str2num(DSRCombined);
end
if ischar(FFCorrection)
    FFCorrection = str2num(FFCorrection);
end
if ischar(BKRemoval)
    BKRemoval = str2num(BKRemoval);
end
if ischar(lowerLimit)
    lowerLimit = str2num(lowerLimit);
end
if ischar(constOffset)
    constOffset = str2num(constOffset);
end
if ischar(FFImagePaths) && ~isempty(FFImagePaths) && strcmp(FFImagePaths(1), '{')
    FFImagePaths = eval(FFImagePaths);
end
if ischar(backgroundPaths) && ~isempty(backgroundPaths) && strcmp(backgroundPaths(1), '{')
    backgroundPaths = eval(backgroundPaths);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(imageListFullpaths) && ~isempty(imageListFullpaths) && strcmp(imageListFullpaths(1), '{')
    imageListFullpaths = eval(imageListFullpaths);
end
if ischar(xcorrShift)
    xcorrShift = str2num(xcorrShift);
end
if ischar(xyMaxOffset)
    xyMaxOffset = str2num(xyMaxOffset);
end
if ischar(zMaxOffset)
    zMaxOffset = str2num(zMaxOffset);
end
if ischar(edgeArtifacts)
    edgeArtifacts = str2num(edgeArtifacts);
end
if ischar(stitchMIP)
    stitchMIP = str2num(stitchMIP);
end
if ischar(onlineStitch)
    onlineStitch = str2num(onlineStitch);
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
if ischar(maxTrialNum)
    maxTrialNum = str2num(maxTrialNum);
end
if ischar(unitWaitTime)
    unitWaitTime = str2num(unitWaitTime);
end
if ischar(minModifyTime)
    minModifyTime = str2num(minModifyTime);
end
if ischar(maxModifyTime)
    maxModifyTime = str2num(maxModifyTime);
end
if ischar(maxWaitLoopNum)
    maxWaitLoopNum = str2num(maxWaitLoopNum);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_microscopeAutomaticProcessing(dataPaths, overwrite=overwrite, streaming=streaming, ...
    channelPatterns=channelPatterns, skewAngle=skewAngle, dz=dz, xyPixelSize=xyPixelSize, ...
    reverse=reverse, objectiveScan=objectiveScan, zStageScan=zStageScan, save16bit=save16bit, ...
    dzFromEncoder=dzFromEncoder, zarrFile=zarrFile, saveZarr=saveZarr, save3DStack=save3DStack, ...
    deskew=deskew, rotate=rotate, stitch=stitch, parseSettingFile=parseSettingFile, ...
    flipZstack=flipZstack, DSRCombined=DSRCombined, FFCorrection=FFCorrection, ...
    BKRemoval=BKRemoval, lowerLimit=lowerLimit, constOffset=constOffset, FFImagePaths=FFImagePaths, ...
    backgroundPaths=backgroundPaths, resampleType=resampleType, resampleFactor=resampleFactor, ...
    inputBbox=inputBbox, stitchResultDirName=stitchResultDirName, imageListFullpaths=imageListFullpaths, ...
    axisOrder=axisOrder, blendMethod=blendMethod, xcorrShift=xcorrShift, xcorrMode=xcorrMode, ...
    xyMaxOffset=xyMaxOffset, zMaxOffset=zMaxOffset, edgeArtifacts=edgeArtifacts, ...
    primaryCh=primaryCh, stitchMIP=stitchMIP, onlineStitch=onlineStitch, generateImageList=generateImageList, ...
    parseCluster=parseCluster, masterCompute=masterCompute, jobLogDir=jobLogDir, ...
    cpusPerTask=cpusPerTask, uuid=uuid, maxTrialNum=maxTrialNum, unitWaitTime=unitWaitTime, ...
    minModifyTime=minModifyTime, maxModifyTime=maxModifyTime, maxWaitLoopNum=maxWaitLoopNum, ...
    mccMode=mccMode, configFile=configFile);

end

