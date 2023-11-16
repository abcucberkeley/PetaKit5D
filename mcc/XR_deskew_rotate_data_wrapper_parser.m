function [] = XR_deskew_rotate_data_wrapper_parser(dataPaths, varargin)
% data-level deskew/rotate wrapper, support small and large scale deskew/rotate
% 
% Author: Xiongtao Ruan (11/15/2023)


ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths');
ip.addParameter('DSDirStr', 'DS/', @ischar);
ip.addParameter('DSRDirStr', 'DSR/', @ischar);
ip.addParameter('Deskew', true, @(x) islogical(x) || ischar(x));
ip.addParameter('Rotate', true, @(x) islogical(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @(x) iscell(x) || ischar(x));
ip.addParameter('dz', 0.5, @(x) isscalar(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isscalar(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('ObjectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ZstageScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('parseSettingFile', false, @(x) islogical(x) || ischar(x)); % use setting file to decide whether filp Z stack or not, it is  poirier over flipZstack
ip.addParameter('Crop', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSRCombined', true, @(x) islogical(x) || ischar(x)); % combined processing 
ip.addParameter('LLFFCorrection', false, @(x) islogical(x) || ischar(x));
ip.addParameter('BKRemoval', false, @(x) islogical(x) || ischar(x));
ip.addParameter('LowerLimit', 0.4, @(x) isnumeric(x) || ischar(x)); % this value is the lowest
ip.addParameter('constOffset', [], @(x) isnumeric(x) || ischar(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('LSImagePaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('BackgroundPaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('Save16bit', false , @(x) islogical(x) || ischar(x)); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('save3DStack', true , @(x) islogical(x) || ischar(x)); % option to save 3D stack or not
ip.addParameter('SaveMIP', true , @(x) islogical(x) || ischar(x)); % save MIP-z for ds and dsr. 
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x)); % use zarr file as input
ip.addParameter('saveZarr', false , @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('BatchSize', [1024, 1024, 1024] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('zarrSubSize', [20, 20, 20], @(x) isnumeric(x) || ischar(x)); % zarr subfolder size
ip.addParameter('inputBbox', [], @(x) isempty(x) || isvector(x) || ischar(x));
ip.addParameter('taskSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x) || ischar(x)); % resampling after rotation 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})) && ischar(x));
ip.addParameter('maskFns', {}, @(x) iscell(x) || ischar(x)); % 2d masks to filter regions to deskew and rotate, in xy, xz, yz order
ip.addParameter('suffix', '', @ischar); % suffix for the folder
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 8, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
DSDirStr = pr.DSDirStr;
DSRDirStr = pr.DSRDirStr;
Deskew = pr.Deskew;
Rotate = pr.Rotate;
Overwrite = pr.Overwrite;
% image related parameters
ChannelPatterns = pr.ChannelPatterns;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
SkewAngle = pr.SkewAngle;
ObjectiveScan = pr.ObjectiveScan;
ZstageScan = pr.ZstageScan;
Reverse = pr.Reverse;
flipZstack = pr.flipZstack;
parseSettingFile = pr.parseSettingFile;
Crop = pr.Crop;
DSRCombined = pr.DSRCombined;
% flat field parameters
LLFFCorrection = pr.LLFFCorrection;
BKRemoval = pr.BKRemoval;
LowerLimit = pr.LowerLimit;
constOffset = pr.constOffset;
LSImagePaths = pr.LSImagePaths;
BackgroundPaths = pr.BackgroundPaths;
% input and output parameters
Save16bit = pr.Save16bit;
save3DStack = pr.save3DStack;
SaveMIP = pr.SaveMIP;
largeFile = pr.largeFile;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
zarrSubSize = pr.zarrSubSize;
inputBbox = pr.inputBbox;
taskSize = pr.taskSize;
resample = pr.resample;
Interp = pr.Interp;
maskFns = pr.maskFns;
suffix = pr.suffix;
% job related parameters
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
debug = pr.debug;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(Deskew)
    Deskew = strcmp(Deskew, 'true');
end
if ischar(Rotate)
    Rotate = strcmp(Rotate, 'true');
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
end
if ischar(ChannelPatterns)
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(ObjectiveScan)
    ObjectiveScan = strcmp(ObjectiveScan, 'true');
end
if ischar(ZstageScan)
    ZstageScan = strcmp(ZstageScan, 'true');
end
if ischar(Reverse)
    Reverse = strcmp(Reverse, 'true');
end
if ischar(flipZstack)
    flipZstack = strcmp(flipZstack, 'true');
end
if ischar(parseSettingFile)
    parseSettingFile = strcmp(parseSettingFile, 'true');
end
if ischar(Crop)
    Crop = strcmp(Crop, 'true');
end
if ischar(DSRCombined)
    DSRCombined = strcmp(DSRCombined, 'true');
end
if ischar(LLFFCorrection)
    LLFFCorrection = strcmp(LLFFCorrection, 'true');
end
if ischar(BKRemoval)
    BKRemoval = strcmp(BKRemoval, 'true');
end
if ischar(LowerLimit)
    LowerLimit = str2num(LowerLimit);
end
if ischar(constOffset)
    constOffset = str2num(constOffset);
end
if ischar(LSImagePaths)
    LSImagePaths = eval(LSImagePaths);
end
if ischar(BackgroundPaths)
    BackgroundPaths = eval(BackgroundPaths);
end
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit, 'true');
end
if ischar(save3DStack)
    save3DStack = strcmp(save3DStack, 'true');
end
if ischar(SaveMIP)
    SaveMIP = strcmp(SaveMIP, 'true');
end
if ischar(largeFile)
    largeFile = strcmp(largeFile, 'true');
end
if ischar(zarrFile)
    zarrFile = strcmp(zarrFile, 'true');
end
if ischar(saveZarr)
    saveZarr = strcmp(saveZarr, 'true');
end
if ischar(BatchSize)
    BatchSize = str2num(BatchSize);
end
if ischar(BlockSize)
    BlockSize = str2num(BlockSize);
end
if ischar(zarrSubSize)
    zarrSubSize = str2num(zarrSubSize);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(taskSize)
    taskSize = str2num(taskSize);
end
if ischar(constOffset)
    resample = str2num(resample);
end
if ischar(maskFns)
    maskFns = eval(maskFns);
end
if ischar(parseCluster)
    parseCluster = strcmp(parseCluster, 'true');
end
if ischar(parseParfor)
    parseParfor = strcmp(parseParfor, 'true');
end
if ischar(masterCompute)
    masterCompute = strcmp(masterCompute, 'true');
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(debug)
    debug = strcmp(debug, 'true');
end
if ischar(mccMode)
    mccMode = strcmp(mccMode, 'true');
end

XR_deskew_rotate_data_wrapper(dataPaths, DSDirStr=DSDirStr, DSRDirStr=DSRDirStr, ...
    Deskew=Deskew, Rotate=Rotate, Overwrite=Overwrite, ChannelPatterns=ChannelPatterns, ...
    dz=dz, xyPixelSize=xyPixelSize, SkewAngle=SkewAngle, ObjectiveScan=ObjectiveScan, ...
    ZstageScan=ZstageScan, Reverse=Reverse, flipZstack=flipZstack, parseSettingFile=parseSettingFile, ...
    Crop=Crop, DSRCombined=DSRCombined, LLFFCorrection=LLFFCorrection, BKRemoval=BKRemoval, ...
    LowerLimit=LowerLimit, constOffset=constOffset, LSImagePaths=LSImagePaths, ...
    BackgroundPaths=BackgroundPaths, Save16bit=Save16bit, save3DStack=save3DStack, ...
    SaveMIP=SaveMIP, largeFile=largeFile, zarrFile=zarrFile, saveZarr=saveZarr, ...
    BatchSize=BatchSize, BlockSize=BlockSize, zarrSubSize=zarrSubSize, inputBbox=inputBbox, ...
    taskSize=taskSize, resample=resample, Interp=Interp, maskFns=maskFns, suffix=suffix, ...
    parseCluster=parseCluster, parseParfor=parseParfor, masterCompute=masterCompute, ...
    jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, uuid=uuid, debug=debug, mccMode=mccMode, ...
    ConfigFile=ConfigFile)


end
