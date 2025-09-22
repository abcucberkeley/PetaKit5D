function [] = XR_chromatic_shift_correction_data_wrapper_parser(dataPaths, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('chromaticOffset', [], @(x) isnumeric(x) || ischar(x)); % y, x, z in voxels
ip.addParameter('resultDirName', 'Chromatic_Shift_Corrected', @ischar); 
ip.addParameter('mode', 'valid', @ischar); % same, valid or full
ip.addParameter('padValue', 0, @(x) isnumeric(x) || ischar(x));
ip.addParameter('newOrigin', [], @(x) isnumeric(x) || ischar(x)); % reset origin (only for same mode).
ip.addParameter('channelPatterns', {'CamA', 'CamB'}, @(x) iscell(x) || ischar(x));
ip.addParameter('psfFullpaths', {}, @(x) iscell(x) || ischar(x));
ip.addParameter('maxOffset', [20, 20, 20], @(x) isnumeric(x) || ischar(x)); % max offset across channels in voxel in y, x, z
ip.addParameter('cropLength', [0, 0, 0], @(x) isnumeric(x) || ischar(x)); % max offset across channels in voxel in y, x, z
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x));
ip.addParameter('batchSize', [1024, 1024, 1024] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
chromaticOffset = pr.chromaticOffset;
resultDirName = pr.resultDirName;
mode = pr.mode;
padValue = pr.padValue;
newOrigin = pr.newOrigin;
channelPatterns = pr.channelPatterns;
psfFullpaths = pr.psfFullpaths;
maxOffset = pr.maxOffset;
cropLength = pr.cropLength;
zarrFile = pr.zarrFile;
largeFile = pr.largeFile;
saveZarr = pr.saveZarr;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
configFile = pr.configFile;
mccMode = pr.mccMode;
uuid = pr.uuid;
debug = pr.debug;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(chromaticOffset)
    chromaticOffset = str2num(chromaticOffset);
end
if ischar(padValue)
    padValue = str2num(padValue);
end
if ischar(newOrigin)
    newOrigin = str2num(newOrigin);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(psfFullpaths) && ~isempty(psfFullpaths) && strcmp(psfFullpaths(1), '{')
    psfFullpaths = eval(psfFullpaths);
end
if ischar(maxOffset)
    maxOffset = str2num(maxOffset);
end
if ischar(cropLength)
    cropLength = str2num(cropLength);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(largeFile)
    largeFile = str2num(largeFile);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
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
if ischar(mccMode)
    mccMode = str2num(mccMode);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_chromatic_shift_correction_data_wrapper(dataPaths, chromaticOffset=chromaticOffset, ...
    resultDirName=resultDirName, mode=mode, padValue=padValue, newOrigin=newOrigin, ...
    channelPatterns=channelPatterns, psfFullpaths=psfFullpaths, maxOffset=maxOffset, ...
    cropLength=cropLength, zarrFile=zarrFile, largeFile=largeFile, saveZarr=saveZarr, ...
    batchSize=batchSize, blockSize=blockSize, parseCluster=parseCluster, masterCompute=masterCompute, ...
    jobLogDir=jobLogDir, cpusPerTask=cpusPerTask, configFile=configFile, mccMode=mccMode, ...
    uuid=uuid, debug=debug);

end

