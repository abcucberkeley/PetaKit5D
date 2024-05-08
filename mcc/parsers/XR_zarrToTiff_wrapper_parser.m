function [] = XR_zarrToTiff_wrapper_parser(dataPaths, varargin)


%#function zarrToTiff

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) iscell(x) || ischar(x));
ip.addParameter('resultDirName', 'tiffs', @ischar);
ip.addParameter('channelPatterns', {'CamA', 'CamB'}, @(x) iscell(x) || ischar(x));
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 30, @(x) isnumeric(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
channelPatterns = pr.channelPatterns;
usrFcn = pr.usrFcn;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(usrFcn) && ~isempty(usrFcn) && (strcmp(usrFcn(1), '{') || strcmp(usrFcn(1), '[') || strcmp(usrFcn(1), '@'))
    usrFcn = eval(usrFcn);
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

XR_zarrToTiff_wrapper(dataPaths, resultDirName=resultDirName, channelPatterns=channelPatterns, ...
    usrFcn=usrFcn, parseCluster=parseCluster, masterCompute=masterCompute, ...
    jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, uuid=uuid, maxTrialNum=maxTrialNum, ...
    unitWaitTime=unitWaitTime, mccMode=mccMode, configFile=configFile);

end

