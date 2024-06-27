function [] = XR_decon_data_wrapper_parser(dataPaths, varargin)


%#function XR_cudaDeconFrame3D
%#function XR_cppDeconFrame3D
%#function XR_RLdeconFrame3D
%#function XR_RotateFrame3D

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x)); % data structure from loadConditionData
ip.addParameter('resultDirName', 'matlab_decon',  @(x) ischar(x));
ip.addParameter('overwrite', false,  @(x) islogical(x) || ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('skewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('dz', 0.5, @(x) isscalar(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isvector(x) && numel(x) <= 2 || ischar(x));
ip.addParameter('save16bit', true, @(x) numel(x) == 1 && islogical(x) || ischar(x));
ip.addParameter('parseSettingFile', false, @(x) islogical(x) || ischar(x)); % use setting file to decide whether filp Z stack or not, it is  poirier over flipZstack
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('background', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('dzPSF', 0.1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('edgeErosion', 0, @(x) isnumeric(x) || ischar(x));
ip.addParameter('erodeByFTP', true, @(x) islogical(x) || ischar(x)); % Edge erosion by the first time point (ranked the first in the inital file list for each dataset).
ip.addParameter('psfFullpaths', {'','',''}, @(x) iscell(x) || ischar(x));
ip.addParameter('deconIter', 15 , @(x) isnumeric(x) || ischar(x)); % number of iterations
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('OTFCumThresh', 0.9, @(x) isnumeric(x) || ischar(x)); % OTF cumutative sum threshold
ip.addParameter('hannWinBounds', [0.8, 1.0], @(x) isnumeric(x) || ischar(x)); % apodization range for distance matrix
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x) || ischar(x)); % decon in skewed space
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x)); % debug mode for simplified code
ip.addParameter('saveStep', 5, @(x) isnumeric(x) || ischar(x)); % save intermediate results every given iterations
ip.addParameter('psfGen', true, @(x) islogical(x) || ischar(x)); % psf generation
ip.addParameter('GPUJob', false, @(x) islogical(x) || ischar(x)); % use gpu for chuck deconvolution. 
ip.addParameter('deconRotate', false, @(x) islogical(x) || ischar(x));
ip.addParameter('batchSize', [1024, 1024, 1024] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x)); % block size 
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeMethod', 'inmemory', @ischar); % inmemory, inplace. 
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x)); % use zarr file as input
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('dampFactor', 1, @(x) isnumeric(x) || ischar(x)); % damp factor for decon result
ip.addParameter('scaleFactor', [], @(x) isnumeric(x) || ischar(x)); % scale factor for decon result
ip.addParameter('deconOffset', 0, @(x) isnumeric(x) || ischar(x)); % offset for decon result
ip.addParameter('maskFullpaths', {}, @(x) iscell(x) || ischar(x)); % 2d masks to filter regions to decon, in xy, xz, yz order
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('unitWaitTime', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);
ip.addParameter('GPUConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
overwrite = pr.overwrite;
channelPatterns = pr.channelPatterns;
skewAngle = pr.skewAngle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
save16bit = pr.save16bit;
parseSettingFile = pr.parseSettingFile;
flipZstack = pr.flipZstack;
background = pr.background;
dzPSF = pr.dzPSF;
edgeErosion = pr.edgeErosion;
erodeByFTP = pr.erodeByFTP;
psfFullpaths = pr.psfFullpaths;
deconIter = pr.deconIter;
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
hannWinBounds = pr.hannWinBounds;
skewed = pr.skewed;
debug = pr.debug;
saveStep = pr.saveStep;
psfGen = pr.psfGen;
GPUJob = pr.GPUJob;
deconRotate = pr.deconRotate;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
largeFile = pr.largeFile;
largeMethod = pr.largeMethod;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
dampFactor = pr.dampFactor;
scaleFactor = pr.scaleFactor;
deconOffset = pr.deconOffset;
maskFullpaths = pr.maskFullpaths;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
unitWaitTime = pr.unitWaitTime;
maxTrialNum = pr.maxTrialNum;
mccMode = pr.mccMode;
configFile = pr.configFile;
GPUConfigFile = pr.GPUConfigFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(overwrite)
    overwrite = str2num(overwrite);
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
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(parseSettingFile)
    parseSettingFile = str2num(parseSettingFile);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(background)
    background = str2num(background);
end
if ischar(dzPSF)
    dzPSF = str2num(dzPSF);
end
if ischar(edgeErosion)
    edgeErosion = str2num(edgeErosion);
end
if ischar(erodeByFTP)
    erodeByFTP = str2num(erodeByFTP);
end
if ischar(psfFullpaths) && ~isempty(psfFullpaths) && strcmp(psfFullpaths(1), '{')
    psfFullpaths = eval(psfFullpaths);
end
if ischar(deconIter)
    deconIter = str2num(deconIter);
end
if ischar(wienerAlpha)
    wienerAlpha = str2num(wienerAlpha);
end
if ischar(OTFCumThresh)
    OTFCumThresh = str2num(OTFCumThresh);
end
if ischar(hannWinBounds)
    hannWinBounds = str2num(hannWinBounds);
end
if ischar(skewed)
    skewed = str2num(skewed);
end
if ischar(debug)
    debug = str2num(debug);
end
if ischar(saveStep)
    saveStep = str2num(saveStep);
end
if ischar(psfGen)
    psfGen = str2num(psfGen);
end
if ischar(GPUJob)
    GPUJob = str2num(GPUJob);
end
if ischar(deconRotate)
    deconRotate = str2num(deconRotate);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
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
if ischar(dampFactor)
    dampFactor = str2num(dampFactor);
end
if ischar(scaleFactor)
    scaleFactor = str2num(scaleFactor);
end
if ischar(deconOffset)
    deconOffset = str2num(deconOffset);
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
if ischar(unitWaitTime)
    unitWaitTime = str2num(unitWaitTime);
end
if ischar(maxTrialNum)
    maxTrialNum = str2num(maxTrialNum);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_decon_data_wrapper(dataPaths, resultDirName=resultDirName, overwrite=overwrite, ...
    channelPatterns=channelPatterns, skewAngle=skewAngle, dz=dz, xyPixelSize=xyPixelSize, ...
    save16bit=save16bit, parseSettingFile=parseSettingFile, flipZstack=flipZstack, ...
    background=background, dzPSF=dzPSF, edgeErosion=edgeErosion, erodeByFTP=erodeByFTP, ...
    psfFullpaths=psfFullpaths, deconIter=deconIter, RLMethod=RLMethod, wienerAlpha=wienerAlpha, ...
    OTFCumThresh=OTFCumThresh, hannWinBounds=hannWinBounds, skewed=skewed, ...
    debug=debug, saveStep=saveStep, psfGen=psfGen, GPUJob=GPUJob, deconRotate=deconRotate, ...
    batchSize=batchSize, blockSize=blockSize, largeFile=largeFile, largeMethod=largeMethod, ...
    zarrFile=zarrFile, saveZarr=saveZarr, dampFactor=dampFactor, scaleFactor=scaleFactor, ...
    deconOffset=deconOffset, maskFullpaths=maskFullpaths, parseCluster=parseCluster, ...
    parseParfor=parseParfor, masterCompute=masterCompute, jobLogDir=jobLogDir, ...
    cpusPerTask=cpusPerTask, uuid=uuid, unitWaitTime=unitWaitTime, maxTrialNum=maxTrialNum, ...
    mccMode=mccMode, configFile=configFile, GPUConfigFile=GPUConfigFile);

end

