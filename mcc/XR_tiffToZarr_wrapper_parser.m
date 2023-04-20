function [] = XR_tiffToZarr_wrapper_parser(tiffFullpaths, varargin)

%#function tiffToZarr

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tiffFullpaths', @(x) iscell(x) || ischar(x));
ip.addParameter('zarrPathstr', 'zarr', @ischar);
ip.addParameter('blockSize', [500, 500, 250], @(x) isnumeric(x) || ischar(x));
ip.addParameter('flippedTile', [], @(x) isempty(x) || islogical(x) || ischar(x));
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x) || ischar(x));
ip.addParameter('partialFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ChannelPatterns', {'tif'}, @(x) iscell(x) || ischar(x));
ip.addParameter('InputBbox', [], @(x) isnumeric(x) || ischar(x)); % crop input tile before processing
ip.addParameter('tileOutBbox', [], @(x) isnumeric(x) || ischar(x)); % crop output tile after processing
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x) || isstring(x) || iscell(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('bigData', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 30, @(x) isnumeric(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);


ip.parse(tiffFullpaths, varargin{:});

pr = ip.Results;
 % Resolution = pr.Resolution;
zarrPathstr = pr.zarrPathstr;
blockSize = pr.blockSize;
flippedTile = pr.flippedTile;
resample = pr.resample;
partialFile = pr.partialFile;
ChannelPatterns = pr.ChannelPatterns;
InputBbox = pr.InputBbox;
tileOutBbox = pr.tileOutBbox;
usrFcn = pr.usrFcn;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
bigData = pr.bigData;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(tiffFullpaths)
    tiffFullpaths = eval(tiffFullpaths);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(flippedTile)
    flippedTile = eval(flippedTile);
end
if ischar(resample)
    resample = eval(resample);
end
if ischar(partialFile)
    partialFile = strcmp(partialFile,'true');
end
if ischar(ChannelPatterns)
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(InputBbox)
    InputBbox = str2num(InputBbox);
end
if ischar(tileOutBbox)
    tileOutBbox = str2num(tileOutBbox);
end
if ischar(parseCluster)
    parseCluster = strcmp(parseCluster,'true');
end
if ischar(bigData)
    bigData = strcmp(bigData,'true');
end
if ischar(masterCompute)
    masterCompute = strcmp(masterCompute,'true');
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(mccMode)
    mccMode = strcmp(mccMode,'true');
end

XR_tiffToZarr_wrapper(tiffFullpaths,'zarrPathstr',zarrPathstr,'blockSize',...
    blockSize,'flippedTile',flippedTile,'resample',resample,'partialFile',partialFile,...
    'ChannelPatterns',ChannelPatterns,'InputBbox',InputBbox,'tileOutBbox',tileOutBbox,...
    'usrFcn',usrFcn,'parseCluster',parseCluster,'bigData',bigData,'masterCompute',masterCompute,...
    'jobLogDir',jobLogDir,'cpusPerTask',cpusPerTask,'uuid',uuid,'maxTrialNum',maxTrialNum,...
    'unitWaitTime',unitWaitTime,'mccMode',mccMode,'ConfigFile',ConfigFile);

end