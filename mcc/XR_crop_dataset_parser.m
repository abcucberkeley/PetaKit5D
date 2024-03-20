function [] = XR_crop_dataset_parser(dataPaths, resultPaths, bbox, varargin)

%#function XR_crop_frame

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('resultPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('bbox', @(x) isnumeric(x) || ischar(x));
ip.addParameter('cropType', 'fixed', @ischar); % fixed or moving or center
ip.addParameter('pad', false, @(x) islogical(x) || ischar(x)); % pad region that is outside the bbox
ip.addParameter('lastStart', [], @(x) isnumeric(x) || ischar(x)); % start coordinate of the last time point
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('zarrFile', false , @(x) islogical(x) || ischar(x)); % read zarr
ip.addParameter('largeZarr', false, @(x) islogical(x) || ischar(x)); % use zarr file as input
ip.addParameter('saveZarr', false , @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('BlockSize', [500, 500, 500] , @(x) isnumeric(x) || ischar(x)); % save as zarr
ip.addParameter('Save16bit', false, @(x) islogical(x) || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPaths, resultPaths, bbox, varargin{:});

pr = ip.Results;
% Overwrite = pr.Overwrite;
cropType = pr.cropType;
pad = pr.pad;
lastStart = pr.lastStart;
ChannelPatterns = pr.ChannelPatterns;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
saveZarr = pr.saveZarr;
BlockSize = pr.BlockSize;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataPaths)
    dataPaths = eval(dataPaths);
end
if ischar(resultPaths)
    resultPaths = eval(resultPaths);
end
if ischar(bbox)
    bbox = str2num(bbox);
end

if ischar(pad)
    pad = strcmp(pad,'true');
end
if ischar(lastStart)
    lastStart = str2num(lastStart);
end
if ischar(ChannelPatterns)
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(zarrFile)
    zarrFile = strcmp(zarrFile,'true');
end
if ischar(largeZarr)
    largeZarr = strcmp(largeZarr,'true');
end
if ischar(saveZarr)
    saveZarr = strcmp(saveZarr,'true');
end
if ischar(BlockSize)
    BlockSize = str2num(BlockSize);
end
if ischar(parseCluster)
    parseCluster = strcmp(parseCluster,'true');
end
if ischar(masterCompute)
    masterCompute = strcmp(masterCompute,'true');
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(mccMode)
    mccMode = strcmp(mccMode, 'true');
end

XR_crop_dataset(dataPaths, resultPaths, bbox, 'cropType', cropType, 'pad', pad,...
    'lastStart', lastStart, 'ChannelPatterns', ChannelPatterns, 'zarrFile', zarrFile, ...
    'largeZarr',largeZarr,'saveZarr',saveZarr, 'BlockSize', BlockSize, 'parseCluster', parseCluster, ...
    'masterCompute', masterCompute, 'jobLogDir', jobLogDir, 'cpusPerTask', cpusPerTask, ...
    uuid=uuid, mccMode=mccMode, ConfigFile=ConfigFile);

end
