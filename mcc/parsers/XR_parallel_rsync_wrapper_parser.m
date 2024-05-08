function [] = XR_parallel_rsync_wrapper_parser(source, dest, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('source', @(x) ischar(x));
ip.addRequired('dest', @(x) ischar(x));
ip.addParameter('cpusPerTask', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('numStream', inf, @(x) isnumeric(x) || ischar(x)); % number of streams for the transfer
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar); % slurm parameters 
ip.addParameter('niceFactor', 0, @(x) isnumeric(x) || ischar(x)); % nice factor for priority 
ip.addParameter('includeSubdir', true, @(x) islogical(x) || ischar(x)); % also transfer subfolders
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(source, dest, varargin{:});

pr = ip.Results;
cpusPerTask = pr.cpusPerTask;
numStream = pr.numStream;
SlurmParam = pr.SlurmParam;
niceFactor = pr.niceFactor;
includeSubdir = pr.includeSubdir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(numStream)
    numStream = str2num(numStream);
end
if ischar(niceFactor)
    niceFactor = str2num(niceFactor);
end
if ischar(includeSubdir)
    includeSubdir = str2num(includeSubdir);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(masterCompute)
    masterCompute = str2num(masterCompute);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_parallel_rsync_wrapper(source, dest, cpusPerTask=cpusPerTask, numStream=numStream, ...
    SlurmParam=SlurmParam, niceFactor=niceFactor, includeSubdir=includeSubdir, ...
    parseCluster=parseCluster, masterCompute=masterCompute, mccMode=mccMode, ...
    configFile=configFile);

end

