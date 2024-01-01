function [] = XR_MIP_wrapper_parser(dataPaths, varargin)

%#function XR_MIP_zarr
%#function saveMIP_zarr
%#function saveMIP_tiff

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('axis', [0, 0, 1], @(x) isnumeric(x) || ischar(x)); % y, x, z
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @(x) iscell(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x)); % use zarr file as input
ip.addParameter('largeZarr', false, @(x) islogical(x) || ischar(x)); % use zarr file as input
ip.addParameter('Save16bit', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs/', @ischar);
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing.
ip.addParameter('cpusPerTask', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;

axis = pr.axis;
ChannelPatterns =  pr.ChannelPatterns;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
Save16bit = pr.Save16bit;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
jobLogDir = pr.jobLogDir;
uuid = pr.uuid;
debug = pr.debug;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataPaths)
    dataPaths = eval(dataPaths);
end

if ischar(axis)
    axis = str2num(axis);
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
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit,'true');
end
if ischar(parseCluster)
    parseCluster = strcmp(parseCluster,'true');
end
if ischar(parseParfor)
    parseParfor = strcmp(parseParfor,'true');
end
if ischar(masterCompute)
    masterCompute = strcmp(masterCompute,'true');
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(debug)
    debug = strcmp(debug,'true');
end
if ischar(mccMode)
    mccMode = strcmp(mccMode, 'true');
end

XR_MIP_wrapper(dataPaths,'axis',axis,'ChannelPatterns',ChannelPatterns,'zarrFile',zarrFile, ...
    'largeZarr',largeZarr,'Save16bit',Save16bit,'parseCluster',parseCluster, ...
    'parseParfor',parseParfor,'masterCompute',masterCompute,'cpusPerTask',cpusPerTask, ...
    'jobLogDir',jobLogDir,'uuid',uuid,'debug',debug, mccMode=mccMode, ConfigFile=ConfigFile);

end