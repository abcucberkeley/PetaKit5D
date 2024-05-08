function [] = XR_RLdeconFrame3D_parser(frameFullpaths, xyPixelSize, dz, deconPath, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('xyPixelSize', @(x) isnumeric(x) || ischar(x));
ip.addRequired('dz', @(x) isnumeric(x) || ischar(x));
ip.addOptional('deconPath', '', @(x) ischar(x) || isempty(x));
ip.addParameter('PSFfile', '', @ischar);
ip.addParameter('Overwrite', false , @(x) islogical(x) || ischar(x));
ip.addParameter('save16bit', false , @(x) islogical(x) || ischar(x));
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
ip.addParameter('batchSize', [1024, 1024, 1024] , @(x) isnumeric(x) || ischar(x)); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x)); % block overlap
ip.addParameter('zarrSubSize', [20, 20, 20], @(x) isnumeric(x) || ischar(x));
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeMethod', 'inmemory', @ischar); % memory jobs, memory single, inplace. 
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('dampFactor', 1, @(x) isnumeric(x) || ischar(x)); % damp factor for decon result
ip.addParameter('scaleFactor', [], @(x) isnumeric(x) || ischar(x)); % scale factor for decon result
ip.addParameter('deconOffset', 0, @(x) isnumeric(x) || ischar(x)); % offset for decon result
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
ip.addParameter('configFile', '', @ischar);
ip.addParameter('GPUConfigFile', '', @ischar);

ip.parse(frameFullpaths, xyPixelSize, dz, deconPath, varargin{:});

pr = ip.Results;
PSFfile = pr.PSFfile;
Overwrite = pr.Overwrite;
save16bit = pr.save16bit;
Rotate = pr.Rotate;
Deskew = pr.Deskew;
SkewAngle = pr.SkewAngle;
flipZstack = pr.flipZstack;
Background = pr.Background;
EdgeSoften = pr.EdgeSoften;
zEdgeSoften = pr.zEdgeSoften;
Crop = pr.Crop;
dzPSF = pr.dzPSF;
DeconIter = pr.DeconIter;
EdgeErosion = pr.EdgeErosion;
ErodeMaskfile = pr.ErodeMaskfile;
SaveMaskfile = pr.SaveMaskfile;
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
skewed = pr.skewed;
fixIter = pr.fixIter;
errThresh = pr.errThresh;
CPUMaxMem = pr.CPUMaxMem;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
zarrSubSize = pr.zarrSubSize;
largeFile = pr.largeFile;
largeMethod = pr.largeMethod;
saveZarr = pr.saveZarr;
dampFactor = pr.dampFactor;
scaleFactor = pr.scaleFactor;
deconOffset = pr.deconOffset;
deconMaskFns = pr.deconMaskFns;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
masterCPU = pr.masterCPU;
GPUJob = pr.GPUJob;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
debug = pr.debug;
saveStep = pr.saveStep;
psfGen = pr.psfGen;
mccMode = pr.mccMode;
configFile = pr.configFile;
GPUConfigFile = pr.GPUConfigFile;

if ischar(frameFullpaths) && ~isempty(frameFullpaths) && strcmp(frameFullpaths(1), '{')
    frameFullpaths = eval(frameFullpaths);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(Rotate)
    Rotate = str2num(Rotate);
end
if ischar(Deskew)
    Deskew = str2num(Deskew);
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(Background)
    Background = str2num(Background);
end
if ischar(EdgeSoften)
    EdgeSoften = str2num(EdgeSoften);
end
if ischar(zEdgeSoften)
    zEdgeSoften = str2num(zEdgeSoften);
end
if ischar(Crop)
    Crop = str2num(Crop);
end
if ischar(dzPSF)
    dzPSF = str2num(dzPSF);
end
if ischar(DeconIter)
    DeconIter = str2num(DeconIter);
end
if ischar(EdgeErosion)
    EdgeErosion = str2num(EdgeErosion);
end
if ischar(SaveMaskfile)
    SaveMaskfile = str2num(SaveMaskfile);
end
if ischar(wienerAlpha)
    wienerAlpha = str2num(wienerAlpha);
end
if ischar(OTFCumThresh)
    OTFCumThresh = str2num(OTFCumThresh);
end
if ischar(skewed)
    skewed = str2num(skewed);
end
if ischar(fixIter)
    fixIter = str2num(fixIter);
end
if ischar(errThresh)
    errThresh = str2num(errThresh);
end
if ischar(CPUMaxMem)
    CPUMaxMem = str2num(CPUMaxMem);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(zarrSubSize)
    zarrSubSize = str2num(zarrSubSize);
end
if ischar(largeFile)
    largeFile = str2num(largeFile);
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
if ischar(deconMaskFns) && ~isempty(deconMaskFns) && strcmp(deconMaskFns(1), '{')
    deconMaskFns = eval(deconMaskFns);
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
if ischar(masterCPU)
    masterCPU = str2num(masterCPU);
end
if ischar(GPUJob)
    GPUJob = str2num(GPUJob);
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
if ischar(debug)
    debug = str2num(debug);
end
if ischar(saveStep)
    saveStep = str2num(saveStep);
end
if ischar(psfGen)
    psfGen = str2num(psfGen);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_RLdeconFrame3D(frameFullpaths, xyPixelSize, dz, deconPath, PSFfile=PSFfile, ...
    Overwrite=Overwrite, save16bit=save16bit, Rotate=Rotate, Deskew=Deskew, ...
    SkewAngle=SkewAngle, flipZstack=flipZstack, Background=Background, EdgeSoften=EdgeSoften, ...
    zEdgeSoften=zEdgeSoften, Crop=Crop, dzPSF=dzPSF, DeconIter=DeconIter, EdgeErosion=EdgeErosion, ...
    ErodeMaskfile=ErodeMaskfile, SaveMaskfile=SaveMaskfile, RLMethod=RLMethod, ...
    wienerAlpha=wienerAlpha, OTFCumThresh=OTFCumThresh, skewed=skewed, fixIter=fixIter, ...
    errThresh=errThresh, CPUMaxMem=CPUMaxMem, batchSize=batchSize, blockSize=blockSize, ...
    zarrSubSize=zarrSubSize, largeFile=largeFile, largeMethod=largeMethod, ...
    saveZarr=saveZarr, dampFactor=dampFactor, scaleFactor=scaleFactor, deconOffset=deconOffset, ...
    deconMaskFns=deconMaskFns, parseCluster=parseCluster, parseParfor=parseParfor, ...
    masterCompute=masterCompute, masterCPU=masterCPU, GPUJob=GPUJob, jobLogDir=jobLogDir, ...
    cpusPerTask=cpusPerTask, uuid=uuid, maxTrialNum=maxTrialNum, unitWaitTime=unitWaitTime, ...
    debug=debug, saveStep=saveStep, psfGen=psfGen, mccMode=mccMode, configFile=configFile, ...
    GPUConfigFile=GPUConfigFile);

end

