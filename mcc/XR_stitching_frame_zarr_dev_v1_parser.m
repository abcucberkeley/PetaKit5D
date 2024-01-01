function [] = XR_stitching_frame_zarr_dev_v1_parser(tileFullpaths, coordinates, varargin)


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
ip.addParameter('dataOrder', 'y,x,z', @ischar); % data axis order, in case data in zyx order
ip.addParameter('flippedTile', [], @(x) all(islogical(x) | isnumeric(x)) || ischar(x));
ip.addParameter('dz', 0.5, @(x) isscalar(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isscalar(x) || ischar(x));
ip.addParameter('ObjectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('IOScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('sCMOSCameraFlip', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Crop', false, @(x) islogical(x) || ischar(x));
ip.addParameter('InputBbox', [], @(x) isnumeric(x) || ischar(x)); % crop input tile before processing
ip.addParameter('tileOutBbox', [], @(x) isnumeric(x) || ischar(x)); % crop tile after processing
ip.addParameter('TileOffset', 0, @(x) isnumeric(x) || ischar(x)); % offset added to the tile
ip.addParameter('df', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('Save16bit', false , @(x) islogical(x) || ischar(x)); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('EdgeArtifacts', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('Decon', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DS', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSR', false, @(x) islogical(x) || ischar(x));
ip.addParameter('resampleType', 'xy_isotropic', @ischar);
ip.addParameter('resample', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('deconRotate', false, @(x) islogical(x) || ischar(x));
ip.addParameter('BlendMethod', 'none', @ischar);
ip.addParameter('blendWeightDegree', 10, @(x) isnumeric(x) || ischar(x));
ip.addParameter('halfOrder', [3, 2, 1], @(x) isnumeric(x) || ischar(x));
ip.addParameter('overlapType', '', @ischar); % '', 'none', 'half', or 'full'
ip.addParameter('xcorrShift', true, @(x) islogical(x) || ischar(x));
ip.addParameter('isPrimaryCh', true, @(x) islogical(x) || ischar(x));
ip.addParameter('usePrimaryCoords', false, @(x) islogical(x) || ischar(x)); % use primary coordinates for secondary channels/tps
ip.addParameter('stitchPadSize', [2, 2, 1], @(x) isnumeric(x) && numel(x) == 3 || ischar(x));
ip.addParameter('padSize', [], @(x) isnumeric(x) && (isempty(x) || numel(x) == 3) || ischar(x));
ip.addParameter('boundboxCrop', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6) || ischar(x));
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
ip.addParameter('largeZarr', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('poolSize', [], @(x) isnumeric(x) || ischar(x)); % max pooling size for large zarr MIPs
ip.addParameter('batchSize', [500, 500, 500], @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('blockSize', [500, 500, 500], @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('shardSize', [], @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('zarrSubSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('saveMultires', false, @(x) islogical(x) || ischar(x)); % save as multi resolution dataset
ip.addParameter('resLevel', 4, @(x) isnumeric(x) || ischar(x)); % downsample to 2^1-2^resLevel
ip.addParameter('BorderSize', [0, 0, 0], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BlurSigma', 10, @(x) isnumeric(x) || ischar(x));
ip.addParameter('SaveMIP', true , @(x) islogical(x) || ischar(x)); % save MIP-z for stitch. 
ip.addParameter('tileIdx', [] , @(x) isnumeric(x) || ischar(x)); % tile indices 
ip.addParameter('processFunPath', '', @(x) isempty(x) || iscell(x) || ischar(x)); % path of user-defined process function handle
ip.addParameter('stitchMIP', [], @(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3)) || ischar(x)); % 1x3 vector or vector, by default, stitch MIP-z
ip.addParameter('stitch2D', false, @(x) islogical(x) || ischar(x));  
ip.addParameter('bigStitchData', false, @(x) islogical(x) || ischar(x));  
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 30, @(x) isnumeric(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(tileFullpaths, coordinates, varargin{:});

pr = ip.Results;
ResultPath = pr.ResultPath;
tileInfoFullpath = pr.tileInfoFullpath;
stitchInfoDir = pr.stitchInfoDir;
stitchInfoFullpath = pr.stitchInfoFullpath;
SkewAngle = pr.SkewAngle;
flippedTile = pr.flippedTile;
axisOrder = pr.axisOrder;
dataOrder = pr.dataOrder;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
ObjectiveScan = pr.ObjectiveScan;
IOScan = pr.IOScan;
sCMOSCameraFlip = pr.sCMOSCameraFlip;
Reverse = pr.Reverse;
Crop = pr.Crop;
InputBbox = pr.InputBbox;
tileOutBbox = pr.tileOutBbox;
TileOffset = pr.TileOffset;
df = pr.df;
ProcessedDirstr = pr.ProcessedDirstr;
Overwrite = pr.Overwrite;
Decon = pr.Decon;
DS = pr.DS;
DSR = pr.DSR;
resample = pr.resample;
resampleType = pr.resampleType;
deconRotate = pr.deconRotate;
BlendMethod = pr.BlendMethod;
blendWeightDegree = pr.blendWeightDegree;
halfOrder = pr.halfOrder;
overlapType = pr.overlapType;
xcorrShift = pr.xcorrShift;
shiftMethod = pr.shiftMethod;
axisWeight = pr.axisWeight;
groupFile = pr.groupFile;
isPrimaryCh = pr.isPrimaryCh;
usePrimaryCoords = pr.usePrimaryCoords;
stitchPadSize = pr.stitchPadSize;
padSize = pr.padSize;
boundboxCrop = pr.boundboxCrop;
zNormalize = pr.zNormalize;
xcorrDownsample = pr.xcorrDownsample;
xcorrThresh = pr.xcorrThresh;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
singleDistMap = pr.singleDistMap;
distBboxes = pr.distBboxes;
saveMultires = pr.saveMultires;
resLevel = pr.resLevel;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
Save16bit = pr.Save16bit;
EdgeArtifacts = pr.EdgeArtifacts;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
poolSize = pr.poolSize;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
shardSize = pr.shardSize;
zarrSubSize = pr.zarrSubSize;
BorderSize = pr.BorderSize;
BlurSigma = pr.BlurSigma;
SaveMIP = pr.SaveMIP;
tileIdx = pr.tileIdx;
processFunPath = pr.processFunPath;
stitchMIP = pr.stitchMIP;
stitch2D = pr.stitch2D;
bigStitchData = pr.bigStitchData;
uuid = pr.uuid;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;
debug = pr.debug;

if ischar(tileFullpaths) && strcmp(tileFullpaths(1), '{')
    tileFullpaths = eval(tileFullpaths);
end
if ischar(coordinates)
    coordinates = str2num(coordinates);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
end
if ischar(SkewAngle)
    SkewAngle = str2double(SkewAngle);
end
if ischar(flippedTile)
    flippedTile = str2num(flippedTile);
end
if ischar(dz)
    dz = str2double(dz);
end
if ischar(xyPixelSize)
    xyPixelSize = str2double(xyPixelSize);
end
if ischar(ObjectiveScan)
    ObjectiveScan = strcmp(ObjectiveScan, 'true');
end
if ischar(IOScan)
    IOScan = strcmp(IOScan, 'true');
end
if ischar(sCMOSCameraFlip)
    sCMOSCameraFlip = strcmp(sCMOSCameraFlip, 'true');
end
if ischar(Reverse)
    Reverse = strcmp(Reverse, 'true');
end
if ischar(Crop)
    Crop = strcmp(Crop, 'true');
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
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit, 'true');
end
if ischar(EdgeArtifacts)
    EdgeArtifacts = str2num(EdgeArtifacts);
end
if ischar(Decon)
    Decon = strcmp(Decon, 'true');
end
if ischar(DS)
    DS = strcmp(DS, 'true');
end
if ischar(DSR)
    DSR = strcmp(DSR, 'true');
end
if ischar(resample)
    resample = str2num(resample);
end
if ischar(deconRotate)
    deconRotate = str2num(deconRotate);
end
if ischar(blendWeightDegree)
    blendWeightDegree = str2double(blendWeightDegree);
end
if ischar(halfOrder)
    halfOrder = str2num(halfOrder);
end
if ischar(xcorrShift)
    xcorrShift = strcmp(xcorrShift, 'true');
end
if ischar(isPrimaryCh)
    isPrimaryCh = strcmp(isPrimaryCh, 'true');
end
if ischar(usePrimaryCoords)
    usePrimaryCoords = strcmp(usePrimaryCoords, 'true');
end
if ischar(stitchPadSize)
    stitchPadSize = str2num(stitchPadSize);
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
if ischar(xcorrDownsample)
    xcorrDownsample = str2num(xcorrDownsample);
end
if ischar(xcorrThresh)
    xcorrThresh = str2double(xcorrThresh);
end
if ischar(xyMaxOffset)
    xyMaxOffset = str2double(xyMaxOffset);
end
if ischar(zMaxOffset)
    zMaxOffset = str2double(zMaxOffset);
end
if ischar(axisWeight)
    axisWeight = str2num(axisWeight);
end
if ischar(singleDistMap)
    singleDistMap = strcmp(singleDistMap, 'true');
end
if ischar(distBboxes)
    distBboxes = str2num(distBboxes);
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
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(shardSize)
    shardSize = str2num(shardSize);
end
if ischar(zarrSubSize)
    zarrSubSize = str2num(zarrSubSize);
end
if ischar(saveMultires)
    saveMultires = strcmp(saveMultires, 'true');
end
if ischar(resLevel)
    resLevel = str2double(resLevel);
end
if ischar(BorderSize)
    BorderSize = str2num(BorderSize);
end
if ischar(BlurSigma)
    BlurSigma = str2double(BlurSigma);
end
if ischar(SaveMIP)
    SaveMIP = strcmp(SaveMIP, 'true');
end
if ischar(tileIdx)
    tileIdx = str2num(tileIdx);
end
if ischar(processFunPath)
    processFunPath = eval(processFunPath);
end
if ischar(stitchMIP)
    stitchMIP = str2num(stitchMIP);
end
if ischar(stitch2D)
    stitch2D = strcmp(stitch2D, 'true');
end
if ischar(bigStitchData)
    bigStitchData = strcmp(bigStitchData, 'true');
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
if ischar(debug)
    debug = strcmp(debug, 'true');
end

XR_stitching_frame_zarr_dev_v1(tileFullpaths, coordinates, ResultPath=ResultPath,...
    tileInfoFullpath=tileInfoFullpath, stitchInfoDir=stitchInfoDir, stitchInfoFullpath=stitchInfoFullpath,...
    ProcessedDirstr=ProcessedDirstr, Overwrite=Overwrite, SkewAngle=SkewAngle, ...
    axisOrder=axisOrder, dataOrder=dataOrder, flippedTile=flippedTile, dz=dz, ...
    xyPixelSize=xyPixelSize, ObjectiveScan=ObjectiveScan, IOScan=IOScan, sCMOSCameraFlip=sCMOSCameraFlip, ...
    Reverse=Reverse, Crop=Crop, InputBbox=InputBbox, tileOutBbox=tileOutBbox, ...
    TileOffset=TileOffset, df=df, Save16bit=Save16bit, EdgeArtifacts=EdgeArtifacts, ...
    Decon=Decon, DS=DS, DSR=DSR, resampleType=resampleType, resample=resample, ...
    deconRotate=deconRotate, BlendMethod=BlendMethod, blendWeightDegree=blendWeightDegree, ...
    halfOrder=halfOrder, overlapType=overlapType, xcorrShift=xcorrShift, isPrimaryCh=isPrimaryCh, ...
    usePrimaryCoords=usePrimaryCoords, stitchPadSize=stitchPadSize, padSize=padSize, ...
    boundboxCrop=boundboxCrop, zNormalize=zNormalize, xcorrDownsample=xcorrDownsample, ...
    xcorrThresh=xcorrThresh, xyMaxOffset=xyMaxOffset, zMaxOffset=zMaxOffset, ...
    shiftMethod=shiftMethod, axisWeight=axisWeight, groupFile=groupFile, singleDistMap=singleDistMap, ...
    distBboxes=distBboxes, zarrFile=zarrFile, largeZarr=largeZarr, poolSize=poolSize, ...
    blockSize=blockSize, batchSize=batchSize, shardSize=shardSize, zarrSubSize=zarrSubSize, ...
    saveMultires=saveMultires, resLevel=resLevel, BorderSize=BorderSize, BlurSigma=BlurSigma, ...
    SaveMIP=SaveMIP, tileIdx=tileIdx, processFunPath=processFunPath, stitchMIP=stitchMIP, ...
    stitch2D=stitch2D, bigStitchData=bigStitchData, parseCluster=parseCluster, ...
    masterCompute=masterCompute, jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, ...
    uuid=uuid, maxTrialNum=maxTrialNum, unitWaitTime=unitWaitTime, mccMode=mccMode, ...
    ConfigFile=ConfigFile, debug=debug);

end

