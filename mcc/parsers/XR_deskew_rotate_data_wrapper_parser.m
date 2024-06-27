function XR_deskew_rotate_data_wrapper_parser(dataPaths, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('DSDirName', 'DS', @ischar);
ip.addParameter('DSRDirName', 'DSR', @ischar);
ip.addParameter('deskew', true, @(x) islogical(x) || ischar(x));
ip.addParameter('rotate', true, @(x) islogical(x) || ischar(x));
ip.addParameter('overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @(x) iscell(x) || ischar(x));
ip.addParameter('dz', 0.5, @(x) isscalar(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isvector(x) && numel(x) <= 2 || ischar(x));
ip.addParameter('skewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('objectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zStageScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('parseSettingFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('crop', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSRCombined', true, @(x) islogical(x) || ischar(x));
ip.addParameter('FFCorrection', false, @(x) islogical(x) || ischar(x));
ip.addParameter('BKRemoval', false, @(x) islogical(x) || ischar(x));
ip.addParameter('lowerLimit', 0.4, @(x) isnumeric(x) || ischar(x));
ip.addParameter('constOffset', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('FFImagePaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('backgroundPaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('save16bit', true , @(x) islogical(x) || ischar(x));
ip.addParameter('save3DStack', true , @(x) islogical(x) || ischar(x));
ip.addParameter('saveMIP', true , @(x) islogical(x) || ischar(x));
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveZarr', false , @(x) islogical(x) || ischar(x));
ip.addParameter('batchSize', [1024, 1024, 1024] , @(x) isvector(x) || ischar(x));
ip.addParameter('blockSize', [256, 256, 256], @(x) isvector(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isempty(x) || isvector(x) || ischar(x));
ip.addParameter('taskSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('resampleFactor', [], @(x) isempty(x) || isnumeric(x) || ischar(x));
ip.addParameter('interpMethod', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})) || ischar(x));
ip.addParameter('maskFullpaths', {}, @(x) iscell(x) || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
DSDirName = pr.DSDirName;
DSRDirName = pr.DSRDirName;
deskew = pr.deskew;
rotate = pr.rotate;
overwrite = pr.overwrite;
channelPatterns = pr.channelPatterns;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
skewAngle = pr.skewAngle;
objectiveScan = pr.objectiveScan;
zStageScan = pr.zStageScan;
reverse = pr.reverse;
flipZstack = pr.flipZstack;
parseSettingFile = pr.parseSettingFile;
crop = pr.crop;
DSRCombined = pr.DSRCombined;
FFCorrection = pr.FFCorrection;
BKRemoval = pr.BKRemoval;
lowerLimit = pr.lowerLimit;
constOffset = pr.constOffset;
FFImagePaths = pr.FFImagePaths;
backgroundPaths = pr.backgroundPaths;
save16bit = pr.save16bit;
save3DStack = pr.save3DStack;
saveMIP = pr.saveMIP;
largeFile = pr.largeFile;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
inputBbox = pr.inputBbox;
taskSize = pr.taskSize;
resampleFactor = pr.resampleFactor;
interpMethod = pr.interpMethod;
maskFullpaths = pr.maskFullpaths;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
debug = pr.debug;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(deskew)
    deskew = str2num(deskew);
end
if ischar(rotate)
    rotate = str2num(rotate);
end
if ischar(overwrite)
    overwrite = str2num(overwrite);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(skewAngle)
    skewAngle = str2num(skewAngle);
end
if ischar(objectiveScan)
    objectiveScan = str2num(objectiveScan);
end
if ischar(zStageScan)
    zStageScan = str2num(zStageScan);
end
if ischar(reverse)
    reverse = str2num(reverse);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(parseSettingFile)
    parseSettingFile = str2num(parseSettingFile);
end
if ischar(crop)
    crop = str2num(crop);
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
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(save3DStack)
    save3DStack = str2num(save3DStack);
end
if ischar(saveMIP)
    saveMIP = str2num(saveMIP);
end
if ischar(largeFile)
    largeFile = str2num(largeFile);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
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
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(taskSize)
    taskSize = str2num(taskSize);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(maskFullpaths) && ~isempty(maskFullpaths) && strcmp(maskFullpaths(1), '{')
    maskFullpaths = eval(maskFullpaths);
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

XR_deskew_rotate_data_wrapper(dataPaths, DSDirName=DSDirName, DSRDirName=DSRDirName, ...
    deskew=deskew, rotate=rotate, overwrite=overwrite, channelPatterns=channelPatterns, ...
    dz=dz, xyPixelSize=xyPixelSize, skewAngle=skewAngle, objectiveScan=objectiveScan, ...
    zStageScan=zStageScan, reverse=reverse, flipZstack=flipZstack, parseSettingFile=parseSettingFile, ...
    crop=crop, DSRCombined=DSRCombined, FFCorrection=FFCorrection, BKRemoval=BKRemoval, ...
    lowerLimit=lowerLimit, constOffset=constOffset, FFImagePaths=FFImagePaths, ...
    backgroundPaths=backgroundPaths, save16bit=save16bit, save3DStack=save3DStack, ...
    saveMIP=saveMIP, largeFile=largeFile, zarrFile=zarrFile, saveZarr=saveZarr, ...
    batchSize=batchSize, blockSize=blockSize, inputBbox=inputBbox, taskSize=taskSize, ...
    resampleFactor=resampleFactor, interpMethod=interpMethod, maskFullpaths=maskFullpaths, ...
    parseCluster=parseCluster, parseParfor=parseParfor, masterCompute=masterCompute, ...
    jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, uuid=uuid, debug=debug, mccMode=mccMode, ...
    configFile=configFile);

end

