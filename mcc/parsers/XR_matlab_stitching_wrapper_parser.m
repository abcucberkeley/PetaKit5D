function [] = XR_matlab_stitching_wrapper_parser(dataPath, imageListFileName, varargin)


%#function XR_stitching_frame_zarr_dev_v1

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @(x) ischar(x) || iscell(x));
ip.addRequired('imageListFileName', @(x) ischar(x) || iscell(x));
ip.addParameter('streaming', false, @(x) islogical(x) || ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('multiLoc', false, @(x) islogical(x) || ischar(x)); % use subregions from different folders
ip.addParameter('processedDirStr', '', @ischar); % path for using existing exist processed data (i.e., DSR, decon)
ip.addParameter('stitchInfoFullpath', '', @ischar); % use exist stitch info for stitching
ip.addParameter('DS', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSR', false, @(x) islogical(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dz', 0.5, @(x) isnumeric(x) || ischar(x));
ip.addParameter('skewAngle', 32.45, @(x) isnumeric(x) || ischar(x));
ip.addParameter('reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('parseSettingFile', false, @(x) islogical(x) || ischar(x)); % use setting file to decide whether filp Z stack or not.
ip.addParameter('axisOrder', 'x,y,z', @ischar); % stage coordinates axis order
ip.addParameter('dataOrder', 'y,x,z', @ischar); % data axis order, in case data in zyx order
ip.addParameter('objectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('IOScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeZarr', false, @(x) islogical(x) || ischar(x));
ip.addParameter('poolSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('batchSize', [500, 500, 500], @(x) isnumeric(x) || ischar(x));
ip.addParameter('blockSize', [500, 500, 500], @(x) isnumeric(x) || ischar(x));
ip.addParameter('shardSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('resampleType', 'xy_isotropic', @ischar); % by default use xy isotropic
ip.addParameter('resampleFactor', [], @(x) isnumeric(x) || ischar(x)); % user-defined resample factor
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x)); % crop input tile before processing
ip.addParameter('tileOutBbox', [], @(x) isnumeric(x) || ischar(x)); % crop tile after processing 
ip.addParameter('tileOffset', 0, @(x) isnumeric(x) || ischar(x)); % offset added to tile
ip.addParameter('resultDirName', 'matlab_stitch', @ischar);
ip.addParameter('blendMethod', 'none', @ischar);
ip.addParameter('overlapType', '', @ischar); % '', 'none', 'half', or 'full'
ip.addParameter('xcorrShift', true, @(x) islogical(x) || ischar(x));
ip.addParameter('xyMaxOffset', 300, @(x) isvector(x) && numel(x) <= 2 || ischar(x)); % max offsets in xy axes
ip.addParameter('zMaxOffset', 50, @(x) isnumeric(x) || ischar(x)); % max offsets in z axis
ip.addParameter('xcorrDownsample', [2, 2, 1], @(x) isnumeric(x) || ischar(x)); % max offsets in z axis
ip.addParameter('xcorrThresh', 0.25, @(x) isnumeric(x) || ischar(x)); % threshold of of xcorr, ignore shift if xcorr below this threshold.
ip.addParameter('padSize', [], @(x) isnumeric(x) && (isempty(x) || numel(x) == 3) || ischar(x));
ip.addParameter('outBbox', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6) || ischar(x));
ip.addParameter('xcorrMode', 'primaryFirst', @(x) ismember(x, {'primary', 'primaryFirst', 'all'}) || ischar(x));  % 'primary': choose one channel as primary channel, 
ip.addParameter('shiftMethod', 'grid', @ischar); % {'local', 'global', 'grid', 'group', 'test'}
ip.addParameter('axisWeight', [1, 0.1, 10], @(x) isnumeric(x) || ischar(x)); % axis weight for optimization, y, x, z
ip.addParameter('groupFile', '', @ischar); % file to define tile groups
ip.addParameter('primaryCh', '', @(x) isempty(x) || ischar(x)); % format: CamA_ch0. If it is empty, use the first channel as primary channel
ip.addParameter('usePrimaryCoords', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('save16bit', false, @(x) islogical(x) || ischar(x));
ip.addParameter('edgeArtifacts', 0, @(x) isnumeric(x) || ischar(x));
ip.addParameter('distBboxes', [], @(x) isnumeric(x) || ischar(x)); % bounding boxes for distance transform
ip.addParameter('saveMIP', true, @(x) islogical(x) || ischar(x));
ip.addParameter('stitchMIP', [], @(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3)) || ischar(x)); % 1x3 vector or vector, by default, stitch MIP-z
ip.addParameter('onlineStitch', false, @(x) islogical(x) || ischar(x)); % support for online stitch (with partial number of tiles). 
ip.addParameter('processFunPath', '', @(x) isempty(x) || ischar(x) || iscell(x)); % path of user-defined process function handle
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 8, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 0.1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPath, imageListFileName, varargin{:});

pr = ip.Results;
streaming = pr.streaming;
channelPatterns = pr.channelPatterns;
multiLoc = pr.multiLoc;
processedDirStr = pr.processedDirStr;
stitchInfoFullpath = pr.stitchInfoFullpath;
DS = pr.DS;
DSR = pr.DSR;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
skewAngle = pr.skewAngle;
reverse = pr.reverse;
parseSettingFile = pr.parseSettingFile;
axisOrder = pr.axisOrder;
dataOrder = pr.dataOrder;
objectiveScan = pr.objectiveScan;
IOScan = pr.IOScan;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
poolSize = pr.poolSize;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
shardSize = pr.shardSize;
resampleType = pr.resampleType;
resampleFactor = pr.resampleFactor;
inputBbox = pr.inputBbox;
tileOutBbox = pr.tileOutBbox;
tileOffset = pr.tileOffset;
resultDirName = pr.resultDirName;
blendMethod = pr.blendMethod;
overlapType = pr.overlapType;
xcorrShift = pr.xcorrShift;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
xcorrDownsample = pr.xcorrDownsample;
xcorrThresh = pr.xcorrThresh;
padSize = pr.padSize;
outBbox = pr.outBbox;
xcorrMode = pr.xcorrMode;
shiftMethod = pr.shiftMethod;
axisWeight = pr.axisWeight;
groupFile = pr.groupFile;
primaryCh = pr.primaryCh;
usePrimaryCoords = pr.usePrimaryCoords;
save16bit = pr.save16bit;
edgeArtifacts = pr.edgeArtifacts;
distBboxes = pr.distBboxes;
saveMIP = pr.saveMIP;
stitchMIP = pr.stitchMIP;
onlineStitch = pr.onlineStitch;
processFunPath = pr.processFunPath;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPath) && ~isempty(dataPath) && strcmp(dataPath(1), '{')
    dataPath = eval(dataPath);
end
if ischar(imageListFileName) && ~isempty(imageListFileName) && strcmp(imageListFileName(1), '{')
    imageListFileName = eval(imageListFileName);
end
if ischar(streaming)
    streaming = str2num(streaming);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(multiLoc)
    multiLoc = str2num(multiLoc);
end
if ischar(DS)
    DS = str2num(DS);
end
if ischar(DSR)
    DSR = str2num(DSR);
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
if ischar(reverse)
    reverse = str2num(reverse);
end
if ischar(parseSettingFile)
    parseSettingFile = str2num(parseSettingFile);
end
if ischar(objectiveScan)
    objectiveScan = str2num(objectiveScan);
end
if ischar(IOScan)
    IOScan = str2num(IOScan);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(largeZarr)
    largeZarr = str2num(largeZarr);
end
if ischar(poolSize)
    poolSize = str2num(poolSize);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(shardSize)
    shardSize = str2num(shardSize);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(tileOutBbox)
    tileOutBbox = str2num(tileOutBbox);
end
if ischar(tileOffset)
    tileOffset = str2num(tileOffset);
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
if ischar(xcorrDownsample)
    xcorrDownsample = str2num(xcorrDownsample);
end
if ischar(xcorrThresh)
    xcorrThresh = str2num(xcorrThresh);
end
if ischar(padSize)
    padSize = str2num(padSize);
end
if ischar(outBbox)
    outBbox = str2num(outBbox);
end
if ischar(axisWeight)
    axisWeight = str2num(axisWeight);
end
if ischar(usePrimaryCoords)
    usePrimaryCoords = str2num(usePrimaryCoords);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(edgeArtifacts)
    edgeArtifacts = str2num(edgeArtifacts);
end
if ischar(distBboxes)
    distBboxes = str2num(distBboxes);
end
if ischar(saveMIP)
    saveMIP = str2num(saveMIP);
end
if ischar(stitchMIP)
    stitchMIP = str2num(stitchMIP);
end
if ischar(onlineStitch)
    onlineStitch = str2num(onlineStitch);
end
if ischar(processFunPath) && ~isempty(processFunPath) && (strcmp(processFunPath(1), '{') || strcmp(processFunPath(1), '[') || strcmp(processFunPath(1), '@'))
    processFunPath = eval(processFunPath);
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
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_matlab_stitching_wrapper(dataPath, imageListFileName, streaming=streaming, ...
    channelPatterns=channelPatterns, multiLoc=multiLoc, processedDirStr=processedDirStr, ...
    stitchInfoFullpath=stitchInfoFullpath, DS=DS, DSR=DSR, xyPixelSize=xyPixelSize, ...
    dz=dz, skewAngle=skewAngle, reverse=reverse, parseSettingFile=parseSettingFile, ...
    axisOrder=axisOrder, dataOrder=dataOrder, objectiveScan=objectiveScan, ...
    IOScan=IOScan, zarrFile=zarrFile, largeZarr=largeZarr, poolSize=poolSize, ...
    batchSize=batchSize, blockSize=blockSize, shardSize=shardSize, resampleType=resampleType, ...
    resampleFactor=resampleFactor, inputBbox=inputBbox, tileOutBbox=tileOutBbox, ...
    tileOffset=tileOffset, resultDirName=resultDirName, blendMethod=blendMethod, ...
    overlapType=overlapType, xcorrShift=xcorrShift, xyMaxOffset=xyMaxOffset, ...
    zMaxOffset=zMaxOffset, xcorrDownsample=xcorrDownsample, xcorrThresh=xcorrThresh, ...
    padSize=padSize, outBbox=outBbox, xcorrMode=xcorrMode, shiftMethod=shiftMethod, ...
    axisWeight=axisWeight, groupFile=groupFile, primaryCh=primaryCh, usePrimaryCoords=usePrimaryCoords, ...
    save16bit=save16bit, edgeArtifacts=edgeArtifacts, distBboxes=distBboxes, ...
    saveMIP=saveMIP, stitchMIP=stitchMIP, onlineStitch=onlineStitch, processFunPath=processFunPath, ...
    parseCluster=parseCluster, masterCompute=masterCompute, jobLogDir=jobLogDir, ...
    cpusPerTask=cpusPerTask, uuid=uuid, maxTrialNum=maxTrialNum, unitWaitTime=unitWaitTime, ...
    mccMode=mccMode, configFile=configFile);

end

