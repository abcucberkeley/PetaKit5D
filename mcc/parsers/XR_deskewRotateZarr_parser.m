function [is_done_flag] = XR_deskewRotateZarr_parser(frameFullpath, xyPixelSize, dz, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpath', @(x) ischar(x) || iscell(x));
ip.addRequired('xyPixelSize', @(x) isscalar(x) || ischar(x)); 
ip.addRequired('dz', @(x) isscalar(x) || ischar(x)); 
ip.addParameter('resultDirStr', 'DSR/', @ischar);
ip.addParameter('ObjectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Crop', false, @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('Reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('DSRCombined', true, @(x) islogical(x) || ischar(x)); % combined processing 
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Save16bit', false , @(x) islogical(x) || ischar(x)); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('SaveMIP', true , @(x) islogical(x) || ischar(x)); % save MIP-z for ds and dsr. 
ip.addParameter('saveZarr', false , @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('BatchSize', [1024, 1024, 1024] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('zarrSubSize', [20, 20, 20], @(x) isnumeric(x) || ischar(x)); % zarr subfolder size
ip.addParameter('inputBbox', [], @(x) isempty(x) || isvector(x) || ischar(x));
ip.addParameter('taskSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x) || ischar(x)); % resampling after rotation 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})) || ischar(x));
ip.addParameter('maskFns', {}, @(x) iscell(x) || ischar(x)); % 2d masks to filter regions to deskew and rotate, in xy, xz, yz order
ip.addParameter('suffix', '', @ischar); % suffix for the folder
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('parseParfor', false, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 8, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(frameFullpath, xyPixelSize, dz, varargin{:});

pr = ip.Results;
resultDirStr = pr.resultDirStr;
ObjectiveScan = pr.ObjectiveScan;
Overwrite = pr.Overwrite;
Crop = pr.Crop;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
DSRCombined = pr.DSRCombined;
flipZstack = pr.flipZstack;
Save16bit = pr.Save16bit;
SaveMIP = pr.SaveMIP;
saveZarr = pr.saveZarr;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
zarrSubSize = pr.zarrSubSize;
inputBbox = pr.inputBbox;
taskSize = pr.taskSize;
resample = pr.resample;
Interp = pr.Interp;
maskFns = pr.maskFns;
suffix = pr.suffix;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
debug = pr.debug;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(frameFullpath) && strcmp(frameFullpath(1), '{')
    frameFullpath = eval(frameFullpath);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(ObjectiveScan)
    ObjectiveScan = str2num(ObjectiveScan);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(Crop)
    Crop = str2num(Crop);
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(Reverse)
    Reverse = str2num(Reverse);
end
if ischar(DSRCombined)
    DSRCombined = str2num(DSRCombined);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(Save16bit)
    Save16bit = str2num(Save16bit);
end
if ischar(SaveMIP)
    SaveMIP = str2num(SaveMIP);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(BatchSize)
    BatchSize = str2num(BatchSize);
end
if ischar(BlockSize)
    BlockSize = str2num(BlockSize);
end
if ischar(zarrSubSize)
    zarrSubSize = str2num(zarrSubSize);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(taskSize)
    taskSize = str2num(taskSize);
end
if ischar(resample)
    resample = str2num(resample);
end
if ischar(maskFns) && strcmp(maskFns(1), '{')
    maskFns = eval(maskFns);
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
    ObjectiveScan=ObjectiveScan, Overwrite=Overwrite, Crop=Crop, SkewAngle=SkewAngle, ...
    Reverse=Reverse, DSRCombined=DSRCombined, flipZstack=flipZstack, Save16bit=Save16bit, ...
    SaveMIP=SaveMIP, saveZarr=saveZarr, BatchSize=BatchSize, BlockSize=BlockSize, ...
    zarrSubSize=zarrSubSize, inputBbox=inputBbox, taskSize=taskSize, resample=resample, ...
    Interp=Interp, maskFns=maskFns, suffix=suffix, parseCluster=parseCluster, ...
    parseParfor=parseParfor, masterCompute=masterCompute, jobLogDir=jobLogDir, ...
    cpusPerTask=cpusPerTask, uuid=uuid, debug=debug, mccMode=mccMode, ConfigFile=ConfigFile);

end

