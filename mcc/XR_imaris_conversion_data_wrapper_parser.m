function [] = XR_imaris_conversion_data_wrapper_parser(dataPaths, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths');
ip.addParameter('imsPathstr', 'imaris',  @(x) ischar(x));
ip.addParameter('Overwrite', false,  @(x) islogical(x) || ischar(x));
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('pixelSizes', [0.108, 0.108, 0.108], @(x) isnumeric(x) || ischar(x)); % y, x, z
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x)); % use zarr file as input
ip.addParameter('blockSize', [64, 64, 64], @(x) isnumeric(x) || ischar(x)); % y, x, z
ip.addParameter('bbox', [], @(x) isnumeric(x) || ischar(x)); % ymin, xmin, zmin, ymax, xmax, zmax
ip.addParameter('timepoints', [], @(x) isnumeric(x) || ischar(x)); % number of time points included
ip.addParameter('ImsConverter', '/clusterfs/fiona/matthewmueller/imarisWriter/writeImarisParallel', @(x) ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 24, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

% parameters
pr = ip.Results;
imsPathstr = pr.imsPathstr;
Overwrite = pr.Overwrite;
ChannelPatterns = pr.ChannelPatterns;
pixelSizes = pr.pixelSizes;
zarrFile = pr.zarrFile;
blockSize = pr.blockSize;
bbox = pr.bbox;
timepoints = pr.timepoints;
ImsConverter = pr.ImsConverter;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
end
if ischar(ChannelPatterns)
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(pixelSizes)
    pixelSizes = str2num(pixelSizes);
end
if ischar(zarrFile)
    zarrFile = strcmp(zarrFile, 'true');
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(bbox)
    bbox = str2num(bbox);
end
if ischar(timepoints)
    timepoints = str2num(timepoints);
end
if ischar(parseCluster)
    parseCluster = strcmp(parseCluster, 'true');
end
if ischar(cpusPerTask)
    cpusPerTask = str2double(cpusPerTask);
end
if ischar(mccMode)
    mccMode = strcmp(mccMode,'true');
end

XR_imaris_conversion_data_wrapper(dataPaths, imsPathstr=imsPathstr, Overwrite=Overwrite, ...
    ChannelPatterns=ChannelPatterns, pixelSizes=pixelSizes, zarrFile=zarrFile, ...
    blockSize=blockSize, bbox=bbox, timepoints=timepoints, ImsConverter=ImsConverter, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, ...
    uuid=uuid, mccMode=mccMode, ConfigFile=ConfigFile);

end

