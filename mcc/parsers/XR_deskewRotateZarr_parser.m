function [is_done_flag] = XR_deskewRotateZarr_parser(frameFullpath, xyPixelSize, dz, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpath', @(x) ischar(x) || iscell(x));
ip.addRequired('xyPixelSize', @(x) isscalar(x) || ischar(x)); 
ip.addRequired('dz', @(x) isscalar(x) || ischar(x)); 
ip.addParameter('resultDirStr', 'DSR/', @ischar);
ip.addParameter('objectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('crop', false, @(x) islogical(x) || ischar(x));
ip.addParameter('skewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSRCombined', true, @(x) islogical(x) || ischar(x)); % combined processing 
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('save16bit', true , @(x) islogical(x) || ischar(x)); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('saveMIP', true , @(x) islogical(x) || ischar(x)); % save MIP-z for ds and dsr. 
ip.addParameter('saveZarr', false , @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('batchSize', [1024, 1024, 1024] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256], @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('inputBbox', [], @(x) isempty(x) || isvector(x) || ischar(x));
ip.addParameter('taskSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('resampleFactor', [], @(x) isempty(x) || isnumeric(x) || ischar(x)); % resampling after rotation 
ip.addParameter('interpMethod', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})) || ischar(x));
ip.addParameter('maskFullpaths', {}, @(x) iscell(x) || ischar(x)); % 2d masks to filter regions to deskew and rotate, in xy, xz, yz order
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 8, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(frameFullpath, xyPixelSize, dz, varargin{:});

pr = ip.Results;
resultDirStr = pr.resultDirStr;
objectiveScan = pr.objectiveScan;
overwrite = pr.overwrite;
crop = pr.crop;
skewAngle = pr.skewAngle;
reverse = pr.reverse;
DSRCombined = pr.DSRCombined;
flipZstack = pr.flipZstack;
save16bit = pr.save16bit;
saveMIP = pr.saveMIP;
saveZarr = pr.saveZarr;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
inputBbox = pr.inputBbox;
taskSize = pr.taskSize;
resampleFactor = pr.resampleFactor;
interpMethod = pr.interpMethod;
maskFullpaths = pr.maskFullpaths;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
debug = pr.debug;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(frameFullpath) && ~isempty(frameFullpath) && strcmp(frameFullpath(1), '{')
    frameFullpath = eval(frameFullpath);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(objectiveScan)
    objectiveScan = str2num(objectiveScan);
end
if ischar(overwrite)
    overwrite = str2num(overwrite);
end
if ischar(crop)
    crop = str2num(crop);
end
if ischar(skewAngle)
    skewAngle = str2num(skewAngle);
end
if ischar(reverse)
    reverse = str2num(reverse);
end
if ischar(DSRCombined)
    DSRCombined = str2num(DSRCombined);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(saveMIP)
    saveMIP = str2num(saveMIP);
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
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(taskSize)
    taskSize = str2num(taskSize);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(maskFullpaths) && ~isempty(maskFullpaths) && strcmp(maskFullpaths(1), '{')
    maskFullpaths = eval(maskFullpaths);
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
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(debug)
    debug = str2num(debug);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_deskewRotateZarr(frameFullpath, xyPixelSize, dz, resultDirStr=resultDirStr, ...
    objectiveScan=objectiveScan, overwrite=overwrite, crop=crop, skewAngle=skewAngle, ...
    reverse=reverse, DSRCombined=DSRCombined, flipZstack=flipZstack, save16bit=save16bit, ...
    saveMIP=saveMIP, saveZarr=saveZarr, batchSize=batchSize, blockSize=blockSize, ...
    inputBbox=inputBbox, taskSize=taskSize, resampleFactor=resampleFactor, ...
    interpMethod=interpMethod, maskFullpaths=maskFullpaths, parseCluster=parseCluster, ...
    parseParfor=parseParfor, masterCompute=masterCompute, jobLogDir=jobLogDir, ...
    cpusPerTask=cpusPerTask, uuid=uuid, debug=debug, mccMode=mccMode, configFile=configFile);

end

