function [] = XR_zarrToTiff_wrapper_parser(dataPaths, varargin)

%#function zarrToTiff

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) iscell(x) || ischar(x));
ip.addParameter('ChannelPatterns', {'CamA', 'CamB'}, @(x) iscell(x) || ischar(x));
ip.addParameter('resultDirStr', 'tiffs/', @ischar);
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 30, @(x) isnumeric(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
% Resolution = pr.Resolution;
ChannelPatterns = pr.ChannelPatterns;
resultDirStr = pr.resultDirStr;
usrFcn = pr.usrFcn;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataPaths)
    dataPaths = eval(dataPaths);
end
if ischar(ChannelPatterns)
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(mccMode)
    mccMode = strcmp(mccMode,'true');
end

XR_zarrToTiff_wrapper(dataPaths,'ChannelPatterns',ChannelPatterns,'resultDirStr',resultDirStr,...
    'usrFcn',usrFcn,'mccMode',mccMode,'ConfigFile',ConfigFile)

end