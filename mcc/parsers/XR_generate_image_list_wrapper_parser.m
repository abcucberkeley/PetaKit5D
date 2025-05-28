function [] = XR_generate_image_list_wrapper_parser(dataPaths, generationMethod, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('generationMethod', @(x) ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('mapTilename', true, @(x) islogical(x) || ischar(x));
ip.addParameter('tilePatterns', {'0000t', 'ch0', '000x', '000y', '000z'}, @(x) iscell(x) || ischar(x));
ip.addParameter('tileFilenames', {}, @(x) iscell(x) || ischar(x));
ip.addParameter('tileIndices', [], @(x) isempty(x) || (isnumeric(x) && size(x, 2) == 5) || ischar(x));
ip.addParameter('tileInterval', [], @(x) isempty(x) || (isnumeric(x) && size(x, 1) == 1 && size(x, 2) == 3) || ischar(x));
ip.addParameter('DS', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSR', false, @(x) islogical(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dz', 0.2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('skewAngle', 32.45, @(x) isnumeric(x) || ischar(x));
ip.addParameter('axisOrder', 'x,y,z', @ischar);
ip.addParameter('dataOrder', 'y,x,z', @ischar);
ip.addParameter('objectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('IOScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('overlapSize', [10, 10, 10], @(x) isnumeric(x) || ischar(x));
ip.addParameter('overlapSizeType', 'pixel', @(x) ischar(x) && ismember(lower(x), {'pixel', 'um'}));
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPaths, generationMethod, varargin{:});

pr = ip.Results;
channelPatterns = pr.channelPatterns;
mapTilename = pr.mapTilename;
tilePatterns = pr.tilePatterns;
tileFilenames = pr.tileFilenames;
tileIndices = pr.tileIndices;
tileInterval = pr.tileInterval;
DS = pr.DS;
DSR = pr.DSR;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
skewAngle = pr.skewAngle;
axisOrder = pr.axisOrder;
dataOrder = pr.dataOrder;
objectiveScan = pr.objectiveScan;
IOScan = pr.IOScan;
zarrFile = pr.zarrFile;
overlapSize = pr.overlapSize;
overlapSizeType = pr.overlapSizeType;
uuid = pr.uuid;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(mapTilename)
    mapTilename = str2num(mapTilename);
end
if ischar(tilePatterns) && ~isempty(tilePatterns) && strcmp(tilePatterns(1), '{')
    tilePatterns = eval(tilePatterns);
end
if ischar(tileFilenames) && ~isempty(tileFilenames) && strcmp(tileFilenames(1), '{')
    tileFilenames = eval(tileFilenames);
end
if ischar(tileIndices)
    tileIndices = str2num(tileIndices);
end
if ischar(tileInterval)
    tileInterval = str2num(tileInterval);
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
if ischar(objectiveScan)
    objectiveScan = str2num(objectiveScan);
end
if ischar(IOScan)
    IOScan = str2num(IOScan);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(overlapSize)
    overlapSize = str2num(overlapSize);
end

XR_generate_image_list_wrapper(dataPaths, generationMethod, channelPatterns=channelPatterns, ...
    mapTilename=mapTilename, tilePatterns=tilePatterns, tileFilenames=tileFilenames, ...
    tileIndices=tileIndices, tileInterval=tileInterval, DS=DS, DSR=DSR, xyPixelSize=xyPixelSize, ...
    dz=dz, skewAngle=skewAngle, axisOrder=axisOrder, dataOrder=dataOrder, objectiveScan=objectiveScan, ...
    IOScan=IOScan, zarrFile=zarrFile, overlapSize=overlapSize, overlapSizeType=overlapSizeType, ...
    uuid=uuid);

end

