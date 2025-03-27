function XR_stitching_frame_zarr_dev_v1_parser(tileFullpaths, coordinates, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tileFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('coordinates', @(x) isnumeric(x) || ischar(x));
ip.addParameter('ResultPath', '', @ischar);
ip.addParameter('tileInfoFullpath', '', @ischar); % matfile contains tileFullpaths and coordinates for too many tiles
ip.addParameter('stitchInfoDir', 'stitchInfo', @ischar);
ip.addParameter('stitchInfoFullpath', '', @ischar); % filename that contain stitch information for secondrary channels
ip.addParameter('ProcessedDirstr', '', @ischar); % path str for processed data 
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('axisOrder', 'x,y,z', @ischar);
ip.addParameter('dataOrder', 'y,x,z', @ischar);
ip.addParameter('flippedTile', [], @(x) isempty(x) || all(islogical(x) | isnumeric(x)) || ischar(x));
ip.addParameter('dz', 0.5, @(x) isscalar(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isscalar(x) || ischar(x));
ip.addParameter('objectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('IOScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('sCMOSCameraFlip', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Crop', false, @(x) islogical(x) || ischar(x));
ip.addParameter('InputBbox', [], @(x) isnumeric(x) || ischar(x)); % crop input tile before processing
ip.addParameter('tileOutBbox', [], @(x) isnumeric(x) || ischar(x)); % crop tile after processing
ip.addParameter('TileOffset', 0, @(x) isnumeric(x) || ischar(x)); % offset added to the tile
ip.addParameter('df', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('save16bit', true , @(x) islogical(x) || ischar(x)); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('EdgeArtifacts', 0, @(x) isnumeric(x) || ischar(x));
ip.addParameter('Decon', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DS', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSR', false, @(x) islogical(x) || ischar(x));
ip.addParameter('resampleType', 'xy_isotropic', @ischar);
ip.addParameter('resampleFactor', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('deconRotate', false, @(x) islogical(x) || ischar(x));
ip.addParameter('BlendMethod', 'none', @ischar);
ip.addParameter('blendWeightDegree', 10, @(x) isnumeric(x) || ischar(x));
ip.addParameter('halfOrder', [3, 2, 1], @(x) isnumeric(x) || ischar(x));
ip.addParameter('overlapType', '', @ischar); % '', 'none', 'half', or 'full'
ip.addParameter('xcorrShift', true, @(x) islogical(x) || ischar(x));
ip.addParameter('isPrimaryCh', true, @(x) islogical(x) || ischar(x));
ip.addParameter('usePrimaryCoords', false, @(x) islogical(x) || ischar(x)); % use primary coordinates for secondary channels/tps
ip.addParameter('stitchPadSize', [2, 2, 2], @(x) isnumeric(x) && numel(x) == 3 || ischar(x));
ip.addParameter('outBbox', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6) || ischar(x));
ip.addParameter('zNormalize', false, @(x) islogical(x) || ischar(x));
ip.addParameter('xcorrDownsample', [2, 2, 1], @(x) isnumeric(x) || ischar(x)); % y,x,z
ip.addParameter('xcorrThresh', 0.25, @(x) isnumeric(x) || ischar(x)); % threshold of of xcorr, ignore shift if xcorr below this threshold.
ip.addParameter('xyMaxOffset', 300, @(x) isnumeric(x) || ischar(x)); % max offsets in xy axes
ip.addParameter('zMaxOffset', 50, @(x) isnumeric(x) || ischar(x)); % max offsets in z axis
ip.addParameter('shiftMethod', 'grid', @ischar); % {'local', 'global', 'grid', 'test', 'group'}
ip.addParameter('axisWeight', [1, 0.1, 10], @(x) isnumeric(x) || ischar(x)); % axis weight for optimization, y, x, z
ip.addParameter('groupFile', '', @ischar); % file to define tile groups
ip.addParameter('singleDistMap', ~false, @(x) islogical(x) || ischar(x)); % compute distance map for the first tile and apply to all other tiles
ip.addParameter('distBboxes', [], @(x) isnumeric(x) || ischar(x)); % bounding boxes for distance transform
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('poolSize', [], @(x) isnumeric(x) || ischar(x)); % max pooling size for large zarr MIPs
ip.addParameter('blockSize', [500, 500, 500], @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('batchSize', [500, 500, 500], @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('shardSize', [], @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('saveMultires', false, @(x) islogical(x) || ischar(x)); % save as multi resolution dataset
ip.addParameter('resLevel', 4, @(x) isnumeric(x) || ischar(x)); % downsample to 2^1-2^resLevel
ip.addParameter('BorderSize', [0, 0, 0], @(x) isnumeric(x) || ischar(x));
ip.addParameter('saveMIP', true , @(x) islogical(x) || ischar(x)); % save MIP-z for stitch. 
ip.addParameter('tileIdx', [] , @(x) isnumeric(x) || ischar(x)); % tile indices 
ip.addParameter('processFunPath', '', @(x) isempty(x) || iscell(x) || ischar(x)); % path of user-defined process function handle
ip.addParameter('stitchMIP', [], @(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3)) || ischar(x)); % 1x3 vector or vector, by default, stitch MIP-z
ip.addParameter('stitch2D', false, @(x)islogical(x) || ischar(x));  
ip.addParameter('bigStitchData', false, @(x)islogical(x) || ischar(x));  
ip.addParameter('maxFileNumPerFolder', 20000, @(x)isscalar(x) || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 30, @(x) isnumeric(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(tileFullpaths, coordinates, varargin{:});

pr = ip.Results;
ResultPath = pr.ResultPath;
tileInfoFullpath = pr.tileInfoFullpath;
stitchInfoDir = pr.stitchInfoDir;
stitchInfoFullpath = pr.stitchInfoFullpath;
ProcessedDirstr = pr.ProcessedDirstr;
Overwrite = pr.Overwrite;
SkewAngle = pr.SkewAngle;
axisOrder = pr.axisOrder;
dataOrder = pr.dataOrder;
flippedTile = pr.flippedTile;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
objectiveScan = pr.objectiveScan;
IOScan = pr.IOScan;
sCMOSCameraFlip = pr.sCMOSCameraFlip;
Reverse = pr.Reverse;
Crop = pr.Crop;
InputBbox = pr.InputBbox;
tileOutBbox = pr.tileOutBbox;
TileOffset = pr.TileOffset;
df = pr.df;
save16bit = pr.save16bit;
EdgeArtifacts = pr.EdgeArtifacts;
Decon = pr.Decon;
DS = pr.DS;
DSR = pr.DSR;
resampleType = pr.resampleType;
resampleFactor = pr.resampleFactor;
deconRotate = pr.deconRotate;
BlendMethod = pr.BlendMethod;
blendWeightDegree = pr.blendWeightDegree;
halfOrder = pr.halfOrder;
overlapType = pr.overlapType;
xcorrShift = pr.xcorrShift;
isPrimaryCh = pr.isPrimaryCh;
usePrimaryCoords = pr.usePrimaryCoords;
stitchPadSize = pr.stitchPadSize;
outBbox = pr.outBbox;
zNormalize = pr.zNormalize;
xcorrDownsample = pr.xcorrDownsample;
xcorrThresh = pr.xcorrThresh;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
shiftMethod = pr.shiftMethod;
axisWeight = pr.axisWeight;
groupFile = pr.groupFile;
singleDistMap = pr.singleDistMap;
distBboxes = pr.distBboxes;
zarrFile = pr.zarrFile;
largeFile = pr.largeFile;
poolSize = pr.poolSize;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
shardSize = pr.shardSize;
saveMultires = pr.saveMultires;
resLevel = pr.resLevel;
BorderSize = pr.BorderSize;
saveMIP = pr.saveMIP;
tileIdx = pr.tileIdx;
processFunPath = pr.processFunPath;
stitchMIP = pr.stitchMIP;
stitch2D = pr.stitch2D;
bigStitchData = pr.bigStitchData;
maxFileNumPerFolder = pr.maxFileNumPerFolder;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
mccMode = pr.mccMode;
configFile = pr.configFile;
debug = pr.debug;

if ischar(tileFullpaths) && ~isempty(tileFullpaths) && strcmp(tileFullpaths(1), '{')
    tileFullpaths = eval(tileFullpaths);
end
if ischar(coordinates)
    coordinates = str2num(coordinates);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(flippedTile)
    flippedTile = str2num(flippedTile);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(objectiveScan)
    objectiveScan = str2num(objectiveScan);
end
if ischar(IOScan)
    IOScan = str2num(IOScan);
end
if ischar(sCMOSCameraFlip)
    sCMOSCameraFlip = str2num(sCMOSCameraFlip);
end
if ischar(Reverse)
    Reverse = str2num(Reverse);
end
if ischar(Crop)
    Crop = str2num(Crop);
end
if ischar(InputBbox)
    InputBbox = str2num(InputBbox);
end
if ischar(tileOutBbox)
    tileOutBbox = str2num(tileOutBbox);
end
if ischar(TileOffset)
    TileOffset = str2num(TileOffset);
end
if ischar(df)
    df = str2num(df);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(EdgeArtifacts)
    EdgeArtifacts = str2num(EdgeArtifacts);
end
if ischar(Decon)
    Decon = str2num(Decon);
end
if ischar(DS)
    DS = str2num(DS);
end
if ischar(DSR)
    DSR = str2num(DSR);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(deconRotate)
    deconRotate = str2num(deconRotate);
end
if ischar(blendWeightDegree)
    blendWeightDegree = str2num(blendWeightDegree);
end
if ischar(halfOrder)
    halfOrder = str2num(halfOrder);
end
if ischar(xcorrShift)
    xcorrShift = str2num(xcorrShift);
end
if ischar(isPrimaryCh)
    isPrimaryCh = str2num(isPrimaryCh);
end
if ischar(usePrimaryCoords)
    usePrimaryCoords = str2num(usePrimaryCoords);
end
if ischar(stitchPadSize)
    stitchPadSize = str2num(stitchPadSize);
end
if ischar(outBbox)
    outBbox = str2num(outBbox);
end
if ischar(zNormalize)
    zNormalize = str2num(zNormalize);
end
if ischar(xcorrDownsample)
    xcorrDownsample = str2num(xcorrDownsample);
end
if ischar(xcorrThresh)
    xcorrThresh = str2num(xcorrThresh);
end
if ischar(xyMaxOffset)
    xyMaxOffset = str2num(xyMaxOffset);
end
if ischar(zMaxOffset)
    zMaxOffset = str2num(zMaxOffset);
end
if ischar(axisWeight)
    axisWeight = str2num(axisWeight);
end
if ischar(singleDistMap)
    singleDistMap = str2num(singleDistMap);
end
if ischar(distBboxes)
    distBboxes = str2num(distBboxes);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(largeFile)
    largeFile = str2num(largeFile);
end
if ischar(poolSize)
    poolSize = str2num(poolSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(shardSize)
    shardSize = str2num(shardSize);
end
if ischar(saveMultires)
    saveMultires = str2num(saveMultires);
end
if ischar(resLevel)
    resLevel = str2num(resLevel);
end
if ischar(BorderSize)
    BorderSize = str2num(BorderSize);
end
if ischar(saveMIP)
    saveMIP = str2num(saveMIP);
end
if ischar(tileIdx)
    tileIdx = str2num(tileIdx);
end
if ischar(processFunPath) && ~isempty(processFunPath) && (strcmp(processFunPath(1), '{') || strcmp(processFunPath(1), '[') || strcmp(processFunPath(1), '@'))
    processFunPath = eval(processFunPath);
end
if ischar(stitchMIP)
    stitchMIP = str2num(stitchMIP);
end
if ischar(stitch2D)
    stitch2D = str2num(stitch2D);
end
if ischar(bigStitchData)
    bigStitchData = str2num(bigStitchData);
end
if ischar(maxFileNumPerFolder)
    maxFileNumPerFolder = str2num(maxFileNumPerFolder);
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
if ischar(debug)
    debug = str2num(debug);
end

XR_stitching_frame_zarr_dev_v1(tileFullpaths, coordinates, ResultPath=ResultPath, ...
    tileInfoFullpath=tileInfoFullpath, stitchInfoDir=stitchInfoDir, stitchInfoFullpath=stitchInfoFullpath, ...
    ProcessedDirstr=ProcessedDirstr, Overwrite=Overwrite, SkewAngle=SkewAngle, ...
    axisOrder=axisOrder, dataOrder=dataOrder, flippedTile=flippedTile, dz=dz, ...
    xyPixelSize=xyPixelSize, objectiveScan=objectiveScan, IOScan=IOScan, sCMOSCameraFlip=sCMOSCameraFlip, ...
    Reverse=Reverse, Crop=Crop, InputBbox=InputBbox, tileOutBbox=tileOutBbox, ...
    TileOffset=TileOffset, df=df, save16bit=save16bit, EdgeArtifacts=EdgeArtifacts, ...
    Decon=Decon, DS=DS, DSR=DSR, resampleType=resampleType, resampleFactor=resampleFactor, ...
    deconRotate=deconRotate, BlendMethod=BlendMethod, blendWeightDegree=blendWeightDegree, ...
    halfOrder=halfOrder, overlapType=overlapType, xcorrShift=xcorrShift, isPrimaryCh=isPrimaryCh, ...
    usePrimaryCoords=usePrimaryCoords, stitchPadSize=stitchPadSize, outBbox=outBbox, ...
    zNormalize=zNormalize, xcorrDownsample=xcorrDownsample, xcorrThresh=xcorrThresh, ...
    xyMaxOffset=xyMaxOffset, zMaxOffset=zMaxOffset, shiftMethod=shiftMethod, ...
    axisWeight=axisWeight, groupFile=groupFile, singleDistMap=singleDistMap, ...
    distBboxes=distBboxes, zarrFile=zarrFile, largeFile=largeFile, poolSize=poolSize, ...
    blockSize=blockSize, batchSize=batchSize, shardSize=shardSize, saveMultires=saveMultires, ...
    resLevel=resLevel, BorderSize=BorderSize, saveMIP=saveMIP, tileIdx=tileIdx, ...
    processFunPath=processFunPath, stitchMIP=stitchMIP, stitch2D=stitch2D, ...
    bigStitchData=bigStitchData, maxFileNumPerFolder=maxFileNumPerFolder, parseCluster=parseCluster, ...
    masterCompute=masterCompute, jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, ...
    uuid=uuid, maxTrialNum=maxTrialNum, unitWaitTime=unitWaitTime, mccMode=mccMode, ...
    configFile=configFile, debug=debug);

end

