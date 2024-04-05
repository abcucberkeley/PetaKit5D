function XR_MIP_zarr_parser(zarrFullpath, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addParameter('mipDirStr', 'MIPs', @ischar); % y, x, z
ip.addParameter('axis', [1, 1, 1], @(x) isnumeric(x) || ischar(x)); % y, x, z
ip.addParameter('BatchSize', [2048, 2048, 2048] , @(x) isnumeric(x) || ischar(x)); % in y, x, z
ip.addParameter('poolSize', [] , @(x) isnumeric(x) || ischar(x)); % pooling size for mips
ip.addParameter('zarrSubSize', [20, 20, 20] , @(x) isnumeric(x) || ischar(x)); % in y, x, z
ip.addParameter('mipSlab', false, @(x) islogical(x) || ischar(x)); % compute MIP slabs (without the final MIPs)
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs/', @ischar);
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(zarrFullpath, varargin{:});

pr = ip.Results;
mipDirStr = pr.mipDirStr;
axis = pr.axis;
BatchSize = pr.BatchSize;
poolSize = pr.poolSize;
zarrSubSize = pr.zarrSubSize;
mipSlab = pr.mipSlab;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;
uuid = pr.uuid;
debug = pr.debug;

if ischar(axis)
    axis = str2num(axis);
end
if ischar(BatchSize)
    BatchSize = str2num(BatchSize);
end
if ischar(poolSize)
    poolSize = str2num(poolSize);
end
if ischar(zarrSubSize)
    zarrSubSize = str2num(zarrSubSize);
end
if ischar(mipSlab)
    mipSlab = str2num(mipSlab);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(parseParfor)
    parseParfor = str2num(parseParfor);
end
if ischar(masterCompute)
    masterCompute = str2num(masterCompute);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_MIP_zarr(zarrFullpath, mipDirStr=mipDirStr, axis=axis, BatchSize=BatchSize, ...
    poolSize=poolSize, zarrSubSize=zarrSubSize, mipSlab=mipSlab, parseCluster=parseCluster, ...
    parseParfor=parseParfor, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    mccMode=mccMode, ConfigFile=ConfigFile, uuid=uuid, debug=debug);

end

