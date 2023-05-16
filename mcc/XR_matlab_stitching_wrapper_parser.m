function [] = XR_matlab_stitching_wrapper_parser(dataPath, imageListFileName, varargin)


%#function XR_matlab_stitching_wrapper

ip = inputParser;
ip.CaseSensitive = false;

ip.addRequired('dataPath', @(x) ischar(x) || iscell(dataPath));
ip.addRequired('imageListFileName', @(x) ischar(x) || iscell(dataPath));
ip.addParameter('Streaming', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('multiLoc', false, @(x) islogical(x) || ischar(x)); % use subregions from different folders
ip.addParameter('ProcessedDirStr', '', @ischar); % path for using existing exist processed data (i.e., DSR, decon)
ip.addParameter('stitchInfoFullpath', '', @ischar); % use exist stitch info for stitching
ip.addParameter('DS', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSR', false, @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45, @(x) isnumeric(x) || ischar(x));
ip.addParameter('Reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('parseSettingFile', false, @(x) islogical(x) || ischar(x)); % use setting file to decide whether filp Z stack or not.
ip.addParameter('axisOrder', 'x,y,z', @ischar);
ip.addParameter('ObjectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('IOScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeZarr', false, @(x) islogical(x) || ischar(x));
ip.addParameter('poolSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('blockSize', [500, 500, 500], @(x) isnumeric(x) || ischar(x));
ip.addParameter('zarrSubSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('resampleType', 'xy_isotropic', @ischar); % by default use xy isotropic
ip.addParameter('resample', [], @(x) isnumeric(x) || ischar(x)); % user-defined resample factor
ip.addParameter('InputBbox', [], @(x) isnumeric(x) || ischar(x)); % crop input tile before processing
ip.addParameter('tileOutBbox', [], @(x) isnumeric(x) || ischar(x)); % crop tile after processing 
ip.addParameter('TileOffset', 0, @(x) isnumeric(x) || ischar(x)); % offset added to tile
ip.addParameter('Resolution', [0.108, 0.5], @(x) isnumeric(x) || ischar(x));
ip.addParameter('resultDir', 'matlab_stitch', @ischar);
ip.addParameter('BlendMethod', 'none', @ischar);
ip.addParameter('overlapType', '', @ischar); % '', 'none', 'half', or 'full'
ip.addParameter('xcorrShift', true, @(x) islogical(x) || ischar(x));
ip.addParameter('xyMaxOffset', 300, @(x) isnumeric(x) || ischar(x)); % max offsets in xy axes
ip.addParameter('zMaxOffset', 50, @(x) isnumeric(x) || ischar(x)); % max offsets in z axis
ip.addParameter('xcorrDownsample', [2, 2, 1], @(x) isnumeric(x) || ischar(x)); % max offsets in z axis
ip.addParameter('xcorrThresh', 0.25, @(x) isnumeric(x) || ischar(x)); % threshold of of xcorr, ignore shift if xcorr below this threshold.
ip.addParameter('padSize', [], @(x) isnumeric(x) && (isempty(x) || numel(x) == 3) || ischar(x));
ip.addParameter('boundboxCrop', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6) || ischar(x));
ip.addParameter('zNormalize', false, @(x) islogical(x) || ischar(x));
ip.addParameter('onlyFirstTP', false, @(x) islogical(x) || ischar(x)); % only compute first time point (for deciding cropping bouding box)
ip.addParameter('timepoints', [], @(x) isnumeric(x) || ischar(x)); % stitch for given time points, nx1 
ip.addParameter('subtimepoints', [], @(x) isnumeric(x) || ischar(x)); % stitch for given sub time points (subtimepoints), nx1
ip.addParameter('xcorrMode', 'primaryFirst', @(x) ismember(x, {'primary', 'primaryFirst', 'all'}) && ischar(x));  % 'primary': choose one channel as primary channel, 
        % 'all': xcorr shift for each channel;  % 'primaryFirst': the primary channel of first time point
ip.addParameter('shiftMethod', 'grid', @ischar); % {'local', 'global', 'grid', 'group', 'test'}
ip.addParameter('axisWeight', [1, 0.1, 10], @(x) isnumeric(x) || ischar(x)); % axis weight for optimization, y, x, z
ip.addParameter('groupFile', '', @ischar); % file to define tile groups
ip.addParameter('primaryCh', '', @(x) isempty(x) || ischar(x)); % format: CamA_ch0. If it is empty, use the first channel as primary channel
ip.addParameter('usePrimaryCoords', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('Save16bit', false, @(x) islogical(x) || ischar(x));
ip.addParameter('EdgeArtifacts', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('saveMIP', true, @(x) islogical(x) || ischar(x));
ip.addParameter('stitchMIP', [], @(x) islogical(x) && (numel(x) == 1 || numel(x) == 3) || ischar(x)); % 1x3 vector or vector, by default, stitch MIP-z
ip.addParameter('onlineStitch', false, @(x) islogical(x) || ischar(x)); % support for online stitch (with partial number of tiles). 
ip.addParameter('bigStitchData', false, @(x) islogical(x) || ischar(x)); % support for online stitch (with partial number of tiles). 
ip.addParameter('pipeline', 'zarr', @(x) strcmpi(x, 'matlab') || strcmpi(x, 'zarr') && ischar(x));
ip.addParameter('processFunPath', '', @(x) isempty(x) || ischar(x) || iscell(x)); % path of user-defined process function handle
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 8, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 0.1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPath, imageListFileName, varargin{:});

pr = ip.Results;
% Overwrite = pr.Overwrite;
Streaming = pr.Streaming;
ChannelPatterns = pr.ChannelPatterns;
multiLoc = pr.multiLoc;
ProcessedDirStr = pr.ProcessedDirStr;
stitchInfoFullpath = pr.stitchInfoFullpath;
DS = pr.DS;
DSR = pr.DSR;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
parseSettingFile = pr.parseSettingFile;
axisOrder = pr.axisOrder;
ObjectiveScan = pr.ObjectiveScan;
IOScan = pr.IOScan;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
poolSize = pr.poolSize;
blockSize = pr.blockSize;
zarrSubSize = pr.zarrSubSize;
resampleType = pr.resampleType;
resample = pr.resample;
InputBbox = pr.InputBbox;
tileOutBbox = pr.tileOutBbox;
TileOffset = pr.TileOffset;
Resolution = pr.Resolution;
resultDir = pr.resultDir;
BlendMethod = pr.BlendMethod;
overlapType = pr.overlapType;
xcorrShift = pr.xcorrShift;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
xcorrDownsample = pr.xcorrDownsample;
xcorrThresh = pr.xcorrThresh;
padSize = pr.padSize;
boundboxCrop = pr.boundboxCrop;
zNormalize = pr.zNormalize;
onlyFirstTP = pr.onlyFirstTP;
timepoints = pr.timepoints;
subtimepoints = pr.subtimepoints;
xcorrMode = pr.xcorrMode;
shiftMethod = pr.shiftMethod;
axisWeight = pr.axisWeight;
groupFile = pr.groupFile;
primaryCh = pr.primaryCh;
usePrimaryCoords = pr.usePrimaryCoords;
Save16bit = pr.Save16bit;
EdgeArtifacts = pr.EdgeArtifacts;
saveMIP = pr.saveMIP;
stitchMIP = pr.stitchMIP;
onlineStitch = pr.onlineStitch;
bigStitchData = pr.bigStitchData;
pipeline = pr.pipeline;
processFunPath = pr.processFunPath;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataPath) && strcmp(dataPath(1), '{')
    dataPath = eval(dataPath);
end
if ischar(imageListFileName) && strcmp(imageListFileName(1), '{')
    imageListFileName = eval(imageListFileName);
end
if ischar(Streaming)
    Streaming = strcmp(Streaming, 'true');
end
if ischar(ChannelPatterns)
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(multiLoc)
    multiLoc = strcmp(multiLoc, 'true');
end
if ischar(DS)
    DS = strcmp(DS, 'true');
end
if ischar(DSR)
    DSR = strcmp(DSR, 'true');
end
if ischar(SkewAngle)
    SkewAngle = str2double(SkewAngle);
end
if ischar(Reverse)
    Reverse = strcmp(Reverse, 'true');
end
if ischar(parseSettingFile)
    parseSettingFile = strcmp(parseSettingFile, 'true');
end
if ischar(ObjectiveScan)
    ObjectiveScan = strcmp(ObjectiveScan, 'true');
end
if ischar(IOScan)
    IOScan = strcmp(IOScan, 'true');
end
if ischar(zarrFile)
    zarrFile = strcmp(zarrFile, 'true');
end
if ischar(largeZarr)
    largeZarr = strcmp(largeZarr, 'true');
end
if ischar(poolSize)
    poolSize = str2num(poolSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(zarrSubSize)
    zarrSubSize = str2num(zarrSubSize);
end
if ischar(resample)
    resample = str2num(resample);
end
if ischar(InputBbox)
    InputBbox = str2num(InputBbox);
end
if ischar(tileOutBbox)
    tileOutBbox = str2num(tileOutBbox);
end
if ischar(TileOffset)
    TileOffset = strcmp(TileOffset,'true');
end
if ischar(Resolution)
    Resolution = str2num(Resolution);
end
if ischar(xcorrShift)
    xcorrShift = strcmp(xcorrShift, 'true');
end
if ischar(xyMaxOffset)
    xyMaxOffset = str2double(xyMaxOffset);
end
if ischar(zMaxOffset)
    zMaxOffset = str2double(zMaxOffset);
end
if ischar(xcorrDownsample)
    xcorrDownsample = str2num(xcorrDownsample);
end
if ischar(xcorrThresh)
    xcorrThresh = str2double(xcorrThresh);
end
if ischar(padSize)
    padSize = str2num(padSize);
end
if ischar(boundboxCrop)
    boundboxCrop = str2num(boundboxCrop);
end
if ischar(zNormalize)
    zNormalize = strcmp(zNormalize, 'true');
end
if ischar(onlyFirstTP)
    onlyFirstTP = strcmp(onlyFirstTP, 'true');
end
if ischar(timepoints)
    timepoints = str2num(timepoints);
end
if ischar(subtimepoints)
    subtimepoints = str2num(subtimepoints);
end
if ischar(axisWeight)
    axisWeight = str2num(axisWeight);
end
if ischar(usePrimaryCoords)
    usePrimaryCoords = strcmp(usePrimaryCoords, 'true');
end
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit, 'true');
end
if ischar(EdgeArtifacts)
    EdgeArtifacts = str2num(EdgeArtifacts);
end
if ischar(saveMIP)
    saveMIP = str2num(saveMIP);
end
if ischar(stitchMIP)
    stitchMIP = str2num(stitchMIP);
end
if ischar(onlineStitch)
    onlineStitch = strcmp(onlineStitch, 'true');
end
if ischar(bigStitchData)
    bigStitchData = strcmp(bigStitchData, 'true');
end
if ischar(processFunPath) && ~isempty(processFunPath) && strcmp(processFunPath(1), '{')
    processFunPath = eval(processFunPath);
end
if ischar(parseCluster)
    parseCluster = strcmp(parseCluster, 'true');
end
if ischar(masterCompute)
    masterCompute = strcmp(masterCompute, 'true');
end
if ischar(cpusPerTask)
    cpusPerTask = str2double(cpusPerTask);
end
if ischar(maxTrialNum)
    maxTrialNum = str2double(maxTrialNum);
end
if ischar(unitWaitTime)
    unitWaitTime = str2double(unitWaitTime);
end
if ischar(mccMode)
    mccMode = strcmp(mccMode, 'true');
end

XR_matlab_stitching_wrapper(dataPath, imageListFileName, Streaming=Streaming,...
    ChannelPatterns=ChannelPatterns, multiLoc=multiLoc, ProcessedDirStr=ProcessedDirStr,...
    stitchInfoFullpath=stitchInfoFullpath, DS=DS, DSR=DSR, SkewAngle=SkewAngle, ...
    Reverse=Reverse, parseSettingFile=parseSettingFile, axisOrder=axisOrder, ...
    ObjectiveScan=ObjectiveScan, IOScan=IOScan, zarrFile=zarrFile, largeZarr=largeZarr, ...
    poolSize=poolSize, blockSize=blockSize, zarrSubSize=zarrSubSize, resampleType=resampleType, ...
    resample=resample, InputBbox=InputBbox, tileOutBbox=tileOutBbox, TileOffset=TileOffset, ...
    Resolution=Resolution,resultDir=resultDir, BlendMethod=BlendMethod, overlapType=overlapType, ...
    xcorrShift=xcorrShift, xyMaxOffset=xyMaxOffset, zMaxOffset=zMaxOffset, ...
    xcorrDownsample=xcorrDownsample, xcorrThresh=xcorrThresh, padSize=padSize, ...
    boundboxCrop=boundboxCrop, zNormalize=zNormalize, onlyFirstTP=onlyFirstTP, ...
    timepoints=timepoints, subtimepoints=subtimepoints, xcorrMode=xcorrMode, ...
    shiftMethod=shiftMethod, axisWeight=axisWeight, groupFile=groupFile, primaryCh=primaryCh, ...
    usePrimaryCoords=usePrimaryCoords, Save16bit=Save16bit,EdgeArtifacts=EdgeArtifacts, ...
    saveMIP=saveMIP, stitchMIP=stitchMIP, onlineStitch=onlineStitch, bigStitchData=bigStitchData, ...
    pipeline=pipeline, processFunPath=processFunPath, parseCluster=parseCluster, ...
    masterCompute=masterCompute, jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, ...
    uuid=uuid, maxTrialNum=maxTrialNum, unitWaitTime=unitWaitTime, mccMode=mccMode, ...
    ConfigFile=ConfigFile);

end

