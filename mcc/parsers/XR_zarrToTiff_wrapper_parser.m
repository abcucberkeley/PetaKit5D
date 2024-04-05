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
ChannelPatterns = pr.ChannelPatterns;
resultDirStr = pr.resultDirStr;
usrFcn = pr.usrFcn;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(ChannelPatterns) && strcmp(ChannelPatterns(1), '{')
    ChannelPatterns = eval(ChannelPatterns);
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

XR_zarrToTiff_wrapper(dataPaths, ChannelPatterns=ChannelPatterns, resultDirStr=resultDirStr, ...
    usrFcn=usrFcn, parseCluster=parseCluster, masterCompute=masterCompute, ...
    jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, uuid=uuid, maxTrialNum=maxTrialNum, ...
    unitWaitTime=unitWaitTime, mccMode=mccMode, ConfigFile=ConfigFile);

end

