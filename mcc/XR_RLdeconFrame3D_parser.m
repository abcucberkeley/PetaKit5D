function [] = XR_RLdeconFrame3D_parser(frameFullpaths, xyPixelSize, dz, varargin)

%#function XR_RLdeconFrame3D

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('xyPixelSize', @(x) isnumeric(x) || ischar(x));
ip.addRequired('dz', @(x) isnumeric(x) || ischar(x));
ip.addOptional('deconPath', '', @(x) ischar(x) || isempty(x));
ip.addParameter('PSFfile', '', @ischar);
ip.addParameter('Overwrite', false , @(x) islogical(x) || ischar(x));
ip.addParameter('Save16bit', false , @(x) islogical(x) || ischar(x));
ip.addParameter('Rotate', false , @(x) islogical(x) || ischar(x));
ip.addParameter('Deskew', false , @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', -32.45 , @(x) isnumeric(x) || ischar(x));
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('Background', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('EdgeSoften', 5, @(x) isnumeric(x) || ischar(x)); % # ofxy px to soften
ip.addParameter('zEdgeSoften', 2, @(x) isnumeric(x) || ischar(x)); % # ofxy px to soften
ip.addParameter('Crop', [], @(x) isnumeric(x) || ischar(x)); % requires lower and higher values for cropping
ip.addParameter('dzPSF', 0.1 , @(x) isnumeric(x) || ischar(x)); %in um
ip.addParameter('DeconIter', 15 , @(x) isnumeric(x) || ischar(x)); % number of iterations
ip.addParameter('EdgeErosion', 8 , @(x) isnumeric(x) || ischar(x)); % erode edges for certain size.
ip.addParameter('ErodeMaskfile', '', @ischar); % erode edges file
ip.addParameter('SaveMaskfile', false, @(x) islogical(x) || ischar(x)); % save mask file for common eroded mask
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('OTFCumThresh', 0.9, @(x) isnumeric(x) || ischar(x)); % OTF cumutative sum threshold
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x) || ischar(x)); % decon in skewed space
ip.addParameter('fixIter', false, @(x) islogical(x) || ischar(x)); % CPU Memory in Gb
ip.addParameter('errThresh', [], @(x) isnumeric(x) || ischar(x)); % error threshold for simplified code
ip.addParameter('CPUMaxMem', 500, @(x) isnumeric(x) || ischar(x)); % CPU Memory in Gb
ip.addParameter('BatchSize', [1024, 1024, 1024] , @(x) isnumeric(x) || ischar(x)); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x)); % block overlap
ip.addParameter('zarrSubSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeMethod', 'MemoryJobs', @ischar); % memory jobs, memory single, inplace. 
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('scaleFactor', [], @(x) isnumeric(x) || ischar(x)); % scale factor for result
ip.addParameter('deconMaskFns', {}, @(x) iscell(x) || ischar(x)); % 2d masks to filter regions to decon, in xy, xz, yz order
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('masterCPU', false, @(x) islogical(x) || ischar(x)); % master node is a cpu node, which is just for large file deconvolution. 
ip.addParameter('GPUJob', false, @(x) islogical(x) || ischar(x)); % use gpu for chuck deconvolution. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 4, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveStep', 5, @(x) isnumeric(x) || ischar(x)); % save intermediate results every given iterations
ip.addParameter('psfGen', true, @(x) islogical(x) || ischar(x)); % psf generation
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);
ip.addParameter('GPUConfigFile', '', @ischar);

ip.parse(frameFullpaths, xyPixelSize, dz, varargin{:});

pr = ip.Results;

deconPath = pr.deconPath;
PSFfile = pr.PSFfile;
Overwrite = pr.Overwrite;
% parameters
dzPSF = pr.dzPSF;
DeconIter = pr.DeconIter;
Deskew = pr.Deskew;
SkewAngle = pr.SkewAngle;
flipZstack = pr.flipZstack;
Rotate = pr.Rotate;
Save16bit = pr.Save16bit;
Crop = pr.Crop;
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
skewed = pr.skewed;
GPUJob = pr.GPUJob;
EdgeSoften = pr.EdgeSoften;
zEdgeSoften = pr.zEdgeSoften;
EdgeErosion = pr.EdgeErosion;
ErodeMaskfile = pr.ErodeMaskfile;
SaveMaskfile = pr.SaveMaskfile;
% check if background information available. Currently use 100. 
Background = pr.Background;
% simplified version related options
fixIter = pr.fixIter;
errThresh = pr.errThresh;
debug = pr.debug;
saveStep = pr.saveStep;
psfGen = pr.psfGen;
CPUMaxMem = pr.CPUMaxMem;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
zarrSubSize = pr.zarrSubSize;
largeFile = pr.largeFile;
largeMethod = pr.largeMethod;
saveZarr = pr.saveZarr;
scaleFactor = pr.scaleFactor;
deconMaskFns = pr.deconMaskFns;

parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
masterCompute = pr.masterCompute;
masterCPU = pr.masterCPU;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
uuid = pr.uuid;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;
GPUConfigFile = pr.GPUConfigFile;

if ischar(frameFullpaths) && strcmp(frameFullpaths(1), '{')
    frameFullpaths = eval(frameFullpaths);
end
if ischar(xyPixelSize)
    xyPixelSize = str2double(xyPixelSize);
end
if ischar(dz)
    dz = str2double(dz);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
end
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit, 'true');
end
if ischar(Rotate)
    Rotate = strcmp(Rotate, 'true');
end
if ischar(Deskew)
    Deskew = strcmp(Deskew, 'true');
end
if ischar(SkewAngle)
    SkewAngle = str2double(SkewAngle);
end
if ischar(flipZstack)
    flipZstack = strcmp(flipZstack, 'true');
end
if ischar(Background)
    Background = str2num(Background);
end
if ischar(EdgeSoften)
    EdgeSoften = str2double(EdgeSoften);
end
if ischar(zEdgeSoften)
    zEdgeSoften = str2double(zEdgeSoften);
end
if ischar(Crop)
    Crop = str2num(Crop);
end
if ischar(dzPSF)
    dzPSF = str2double(dzPSF);
end
if ischar(DeconIter)
    DeconIter = str2double(DeconIter);
end
if ischar(EdgeErosion)
    EdgeErosion = str2double(EdgeErosion);
end
if ischar(SaveMaskfile)
    SaveMaskfile = strcmp(SaveMaskfile, 'true');
end
if ischar(wienerAlpha)
    wienerAlpha = str2double(wienerAlpha);
end
if ischar(OTFCumThresh)
    OTFCumThresh = str2double(OTFCumThresh);
end
if ischar(skewed)
    if strcmp(skewed(1), '[')
        skewed = [];
    else
        skewed = strcmp(skewed, 'true');
    end
end
if ischar(fixIter)
    fixIter = strcmp(fixIter, 'true');
end
if ischar(errThresh)
    errThresh = str2double(errThresh);
end
if ischar(CPUMaxMem)
    CPUMaxMem = str2double(CPUMaxMem);
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
if ischar(largeFile)
    largeFile = strcmp(largeFile, 'true');
end
if ischar(saveZarr)
    saveZarr = strcmp(saveZarr, 'true');
end
if ischar(scaleFactor)
    scaleFactor = str2num(scaleFactor);
end
if ischar(deconMaskFns)
    deconMaskFns = eval(deconMaskFns);
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
if ischar(masterCPU)
    masterCPU = strcmp(masterCPU, 'true');
end
if ischar(GPUJob)
    GPUJob = strcmp(GPUJob, 'true');
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
if ischar(debug)
    debug = strcmp(debug, 'true');
end
if ischar(saveStep)
    saveStep = str2double(saveStep);
end
if ischar(psfGen)
    psfGen = strcmp(psfGen, 'true');
end
if ischar(mccMode)
    mccMode = strcmp(mccMode, 'true');
end

XR_RLdeconFrame3D(frameFullpaths, xyPixelSize, dz, deconPath, PSFfile=PSFfile, ...
    Overwrite=Overwrite, Save16bit=Save16bit, Rotate=Rotate, Deskew=Deskew, ...
    SkewAngle=SkewAngle, flipZstack=flipZstack, Background=Background, EdgeSoften=EdgeSoften, ...
    zEdgeSoften=zEdgeSoften, Crop=Crop, dzPSF=dzPSF, DeconIter=DeconIter, EdgeErosion=EdgeErosion, ...
    ErodeMaskfile=ErodeMaskfile, SaveMaskfile=SaveMaskfile, RLMethod=RLMethod, ...
    wienerAlpha=wienerAlpha, OTFCumThresh=OTFCumThresh, skewed=skewed, fixIter=fixIter, ...
    errThresh=errThresh, CPUMaxMem=CPUMaxMem, BatchSize=BatchSize, BlockSize=BlockSize, ...
    zarrSubSize=zarrSubSize, largeFile=largeFile, largeMethod=largeMethod, saveZarr=saveZarr, ...
    scaleFactor=scaleFactor, deconMaskFns=deconMaskFns, parseCluster=parseCluster, ...
    parseParfor=parseParfor, masterCompute=masterCompute, masterCPU=masterCPU, ...
    GPUJob=GPUJob, jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, uuid=uuid, ...
    maxTrialNum=maxTrialNum, unitWaitTime=unitWaitTime, debug=debug, saveStep=saveStep, ...
    psfGen=psfGen, mccMode=mccMode, ConfigFile=ConfigFile, GPUConfigFile=GPUConfigFile);

end

