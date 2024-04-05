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
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x)); % data structure from loadConditionData
ip.addParameter('Overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 5) && islogical(x) || ischar(x));
ip.addParameter('Streaming', true,  @(x) islogical(x) || ischar(x)); % if true, check for new files. If false, assume all files transferred completely.
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('Channels', [488, 560, 642], @(x) isnumeric(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('dz', 0.5, @(x) isscalar(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isscalar(x) || ischar(x));
ip.addParameter('Reverse', true, @(x) islogical(x) || ischar(x));
ip.addParameter('ObjectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ZstageScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('sCMOSCameraFlip', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Save16bit', [false, false, false, false], @(x) (numel(x) == 1 || numel(x) == 4) && islogical(x) || ischar(x));
ip.addParameter('onlyFirstTP', false, @(x) islogical(x) || ischar(x));
ip.addParameter('dzFromEncoder', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x)); % use zarr file as input
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); % use zarr file as output
ip.addParameter('save3DStack', true , @(x) islogical(x) || ischar(x)); % option to save 3D stack or not
ip.addParameter('Deskew', true, @(x) islogical(x) || ischar(x));
ip.addParameter('Rotate', true, @(x) islogical(x) || ischar(x));
ip.addParameter('Stitch', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Decon', ~false, @(x) islogical(x) || ischar(x));
ip.addParameter('RotateAfterDecon', false, @(x) islogical(x) || ischar(x));
ip.addParameter('parseSettingFile', false, @(x) islogical(x) || ischar(x)); % use setting file to decide whether filp Z stack or not.
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x)); % 
ip.addParameter('DSRCombined', true, @(x) islogical(x) || ischar(x)); 
ip.addParameter('LLFFCorrection', false, @(x) islogical(x) || ischar(x));
ip.addParameter('BKRemoval', false, @(x) islogical(x) || ischar(x));
ip.addParameter('LowerLimit', 0.4, @(x) isnumeric(x) || ischar(x)); % this value is the lowest
ip.addParameter('constOffset', [], @(x) isnumeric(x) || ischar(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('LSImagePaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('BackgroundPaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('resampleType', 'isotropic', @ischar); % resample type: given, isotropic, xy_isotropic
ip.addParameter('resample', [], @(x) isnumeric(x) || ischar(x)); % resample
ip.addParameter('InputBbox', [], @(x) isnumeric(x) || ischar(x)); % bbox for input in deskew and rotate
ip.addParameter('stitchPipeline', 'zarr', @ischar); % matlab or zarr
ip.addParameter('stitchResultDir', '', @ischar);
ip.addParameter('imageListFullpaths', '', @(x) ischar(x) || iscell(x));
ip.addParameter('axisOrder', 'xyz', @(x) ischar(x));
ip.addParameter('BlendMethod', 'none', @ischar);
ip.addParameter('xcorrShift', false, @(x) islogical(x) || ischar(x));
ip.addParameter('xcorrMode', 'primaryFirst', @(x) ismember(lower(x), {'primary', 'primaryfirst', 'all'}) || ischar(x)); % 'primary': choose one channel as primary channel, 
ip.addParameter('xyMaxOffset', 300, @(x) isnumeric(x) || ischar(x)); % max offsets in xy axes
ip.addParameter('zMaxOffset', 50, @(x) isnumeric(x) || ischar(x)); % max offsets in z axis
ip.addParameter('EdgeArtifacts', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('timepoints', [], @(x) isnumeric(x) || ischar(x)); % stitch for given time points
ip.addParameter('boundboxCrop', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6) || ischar(x));
ip.addParameter('primaryCh', '', @ischar);
ip.addParameter('stitchMIP', [], @(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3)) || ischar(x)); % 1x3 vector or vector, by default, stitch MIP-z
ip.addParameter('onlineStitch', false, @(x) islogical(x) || ischar(x)); % support for online stitch (with partial number of tiles). 
ip.addParameter('generateImageList', '', @(x) ischar(x)); % for real time processing, {'', 'from_encoder', 'from_sqlite'}
ip.addParameter('cudaDecon', false, @(x) islogical(x) || ischar(x));
ip.addParameter('cppDecon', false, @(x) islogical(x) || ischar(x));
ip.addParameter('cppDeconPath', '/global/home/groups/software/sl-7.x86_64/modules/RLDecon_CPU/20200718/build-cluster/cpuDeconv', @ischar);
ip.addParameter('loadModules', 'module load gcc/4.8.5; module load fftw/3.3.6-gcc; module load boost/1.65.1-gcc; module load libtiff/4.1.0; ', @ischar);
ip.addParameter('cudaDeconPath', '/global/home/groups/software/sl-7.x86_64/modules/cudaDecon/bin/cudaDeconv' , @ischar);
ip.addParameter('OTFGENPath', '/global/home/groups/software/sl-7.x86_64/modules/cudaDecon/bin/radialft' , @ischar); % point to radialft file
ip.addParameter('DS', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSR', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Background', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('dzPSF', 0.1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('EdgeErosion', 8, @(x) isnumeric(x) || ischar(x));
ip.addParameter('ErodeByFTP', true, @(x) islogical(x) || ischar(x)); % Edge erosion by the first time point (ranked the first in the inital file list for each dataset).
ip.addParameter('deconRotate', false, @(x) islogical(x) || ischar(x));
ip.addParameter('psfFullpaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('DeconIter', 15 , @(x) isnumeric(x) || ischar(x)); % number of iterations
ip.addParameter('rotatedPSF', false , @(x) islogical(x) || ischar(x)); % psf is rotated (for dsr)
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('fixIter', false, @(x) islogical(x) || ischar(x)); % CPU Memory in Gb
ip.addParameter('errThresh', [], @(x) isnumeric(x) || ischar(x)); % error threshold for simplified code
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x)); % debug mode for simplified code
ip.addParameter('GPUJob', false, @(x) islogical(x) || ischar(x)); % use gpu for chuck deconvolution. 
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('minModifyTime', 1, @(x) isnumeric(x) || ischar(x)); % the minimum duration of last modify time of a file, in minute.
ip.addParameter('maxModifyTime', 10, @(x) isnumeric(x) || ischar(x)); % the maximum duration of last modify time of a file, in minute.
ip.addParameter('maxWaitLoopNum', 10, @(x) isnumeric(x) || ischar(x)); % the max number of loops the loop waits with all existing files processed. 
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);
ip.addParameter('GPUConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
Streaming = pr.Streaming;
ChannelPatterns = pr.ChannelPatterns;
Channels = pr.Channels;
SkewAngle = pr.SkewAngle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
Reverse = pr.Reverse;
ObjectiveScan = pr.ObjectiveScan;
ZstageScan = pr.ZstageScan;
sCMOSCameraFlip = pr.sCMOSCameraFlip;
Save16bit = pr.Save16bit;
onlyFirstTP = pr.onlyFirstTP;
dzFromEncoder = pr.dzFromEncoder;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
save3DStack = pr.save3DStack;
Deskew = pr.Deskew;
Rotate = pr.Rotate;
Stitch = pr.Stitch;
Decon = pr.Decon;
RotateAfterDecon = pr.RotateAfterDecon;
parseSettingFile = pr.parseSettingFile;
flipZstack = pr.flipZstack;
DSRCombined = pr.DSRCombined;
LLFFCorrection = pr.LLFFCorrection;
BKRemoval = pr.BKRemoval;
LowerLimit = pr.LowerLimit;
constOffset = pr.constOffset;
LSImagePaths = pr.LSImagePaths;
BackgroundPaths = pr.BackgroundPaths;
resampleType = pr.resampleType;
resample = pr.resample;
InputBbox = pr.InputBbox;
stitchPipeline = pr.stitchPipeline;
stitchResultDir = pr.stitchResultDir;
imageListFullpaths = pr.imageListFullpaths;
axisOrder = pr.axisOrder;
BlendMethod = pr.BlendMethod;
xcorrShift = pr.xcorrShift;
xcorrMode = pr.xcorrMode;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
EdgeArtifacts = pr.EdgeArtifacts;
timepoints = pr.timepoints;
boundboxCrop = pr.boundboxCrop;
primaryCh = pr.primaryCh;
stitchMIP = pr.stitchMIP;
onlineStitch = pr.onlineStitch;
generateImageList = pr.generateImageList;
cudaDecon = pr.cudaDecon;
cppDecon = pr.cppDecon;
cppDeconPath = pr.cppDeconPath;
loadModules = pr.loadModules;
cudaDeconPath = pr.cudaDeconPath;
OTFGENPath = pr.OTFGENPath;
DS = pr.DS;
DSR = pr.DSR;
Background = pr.Background;
dzPSF = pr.dzPSF;
EdgeErosion = pr.EdgeErosion;
ErodeByFTP = pr.ErodeByFTP;
deconRotate = pr.deconRotate;
psfFullpaths = pr.psfFullpaths;
DeconIter = pr.DeconIter;
rotatedPSF = pr.rotatedPSF;
RLMethod = pr.RLMethod;
fixIter = pr.fixIter;
errThresh = pr.errThresh;
debug = pr.debug;
GPUJob = pr.GPUJob;
largeFile = pr.largeFile;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
minModifyTime = pr.minModifyTime;
maxModifyTime = pr.maxModifyTime;
maxWaitLoopNum = pr.maxWaitLoopNum;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;
GPUConfigFile = pr.GPUConfigFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(Streaming)
    Streaming = str2num(Streaming);
end
if ischar(ChannelPatterns) && ~isempty(ChannelPatterns) && strcmp(ChannelPatterns(1), '{')
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(Channels)
    Channels = str2num(Channels);
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(Reverse)
    Reverse = str2num(Reverse);
end
if ischar(ObjectiveScan)
    ObjectiveScan = str2num(ObjectiveScan);
end
if ischar(ZstageScan)
    ZstageScan = str2num(ZstageScan);
end
if ischar(sCMOSCameraFlip)
    sCMOSCameraFlip = str2num(sCMOSCameraFlip);
end
if ischar(Save16bit)
    Save16bit = str2num(Save16bit);
end
if ischar(onlyFirstTP)
    onlyFirstTP = str2num(onlyFirstTP);
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
if ischar(Deskew)
    Deskew = str2num(Deskew);
end
if ischar(Rotate)
    Rotate = str2num(Rotate);
end
if ischar(Stitch)
    Stitch = str2num(Stitch);
end
if ischar(Decon)
    Decon = str2num(Decon);
end
if ischar(RotateAfterDecon)
    RotateAfterDecon = str2num(RotateAfterDecon);
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
if ischar(LLFFCorrection)
    LLFFCorrection = str2num(LLFFCorrection);
end
if ischar(BKRemoval)
    BKRemoval = str2num(BKRemoval);
end
if ischar(LowerLimit)
    LowerLimit = str2num(LowerLimit);
end
if ischar(constOffset)
    constOffset = str2num(constOffset);
end
if ischar(LSImagePaths) && ~isempty(LSImagePaths) && strcmp(LSImagePaths(1), '{')
    LSImagePaths = eval(LSImagePaths);
end
if ischar(BackgroundPaths) && ~isempty(BackgroundPaths) && strcmp(BackgroundPaths(1), '{')
    BackgroundPaths = eval(BackgroundPaths);
end
if ischar(resample)
    resample = str2num(resample);
end
if ischar(InputBbox)
    InputBbox = str2num(InputBbox);
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
if ischar(EdgeArtifacts)
    EdgeArtifacts = str2num(EdgeArtifacts);
end
if ischar(timepoints)
    timepoints = str2num(timepoints);
end
if ischar(boundboxCrop)
    boundboxCrop = str2num(boundboxCrop);
end
if ischar(stitchMIP)
    stitchMIP = str2num(stitchMIP);
end
if ischar(onlineStitch)
    onlineStitch = str2num(onlineStitch);
end
if ischar(cudaDecon)
    cudaDecon = str2num(cudaDecon);
end
if ischar(cppDecon)
    cppDecon = str2num(cppDecon);
end
if ischar(DS)
    DS = str2num(DS);
end
if ischar(DSR)
    DSR = str2num(DSR);
end
if ischar(Background)
    Background = str2num(Background);
end
if ischar(dzPSF)
    dzPSF = str2num(dzPSF);
end
if ischar(EdgeErosion)
    EdgeErosion = str2num(EdgeErosion);
end
if ischar(ErodeByFTP)
    ErodeByFTP = str2num(ErodeByFTP);
end
if ischar(deconRotate)
    deconRotate = str2num(deconRotate);
end
if ischar(psfFullpaths) && ~isempty(psfFullpaths) && strcmp(psfFullpaths(1), '{')
    psfFullpaths = eval(psfFullpaths);
end
if ischar(DeconIter)
    DeconIter = str2num(DeconIter);
end
if ischar(rotatedPSF)
    rotatedPSF = str2num(rotatedPSF);
end
if ischar(fixIter)
    fixIter = str2num(fixIter);
end
if ischar(errThresh)
    errThresh = str2num(errThresh);
end
if ischar(debug)
    debug = str2num(debug);
end
if ischar(GPUJob)
    GPUJob = str2num(GPUJob);
end
if ischar(largeFile)
    largeFile = str2num(largeFile);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
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

XR_microscopeAutomaticProcessing(dataPaths, Overwrite=Overwrite, Streaming=Streaming, ...
    ChannelPatterns=ChannelPatterns, Channels=Channels, SkewAngle=SkewAngle, ...
    dz=dz, xyPixelSize=xyPixelSize, Reverse=Reverse, ObjectiveScan=ObjectiveScan, ...
    ZstageScan=ZstageScan, sCMOSCameraFlip=sCMOSCameraFlip, Save16bit=Save16bit, ...
    onlyFirstTP=onlyFirstTP, dzFromEncoder=dzFromEncoder, zarrFile=zarrFile, ...
    saveZarr=saveZarr, save3DStack=save3DStack, Deskew=Deskew, Rotate=Rotate, ...
    Stitch=Stitch, Decon=Decon, RotateAfterDecon=RotateAfterDecon, parseSettingFile=parseSettingFile, ...
    flipZstack=flipZstack, DSRCombined=DSRCombined, LLFFCorrection=LLFFCorrection, ...
    BKRemoval=BKRemoval, LowerLimit=LowerLimit, constOffset=constOffset, LSImagePaths=LSImagePaths, ...
    BackgroundPaths=BackgroundPaths, resampleType=resampleType, resample=resample, ...
    InputBbox=InputBbox, stitchPipeline=stitchPipeline, stitchResultDir=stitchResultDir, ...
    imageListFullpaths=imageListFullpaths, axisOrder=axisOrder, BlendMethod=BlendMethod, ...
    xcorrShift=xcorrShift, xcorrMode=xcorrMode, xyMaxOffset=xyMaxOffset, zMaxOffset=zMaxOffset, ...
    EdgeArtifacts=EdgeArtifacts, timepoints=timepoints, boundboxCrop=boundboxCrop, ...
    primaryCh=primaryCh, stitchMIP=stitchMIP, onlineStitch=onlineStitch, generateImageList=generateImageList, ...
    cudaDecon=cudaDecon, cppDecon=cppDecon, cppDeconPath=cppDeconPath, loadModules=loadModules, ...
    cudaDeconPath=cudaDeconPath, OTFGENPath=OTFGENPath, DS=DS, DSR=DSR, Background=Background, ...
    dzPSF=dzPSF, EdgeErosion=EdgeErosion, ErodeByFTP=ErodeByFTP, deconRotate=deconRotate, ...
    psfFullpaths=psfFullpaths, DeconIter=DeconIter, rotatedPSF=rotatedPSF, ...
    RLMethod=RLMethod, fixIter=fixIter, errThresh=errThresh, debug=debug, GPUJob=GPUJob, ...
    largeFile=largeFile, parseCluster=parseCluster, jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, ...
    uuid=uuid, maxTrialNum=maxTrialNum, unitWaitTime=unitWaitTime, minModifyTime=minModifyTime, ...
    maxModifyTime=maxModifyTime, maxWaitLoopNum=maxWaitLoopNum, mccMode=mccMode, ...
    ConfigFile=ConfigFile, GPUConfigFile=GPUConfigFile);

end

