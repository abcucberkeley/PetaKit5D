function XR_deskewRotateZarr_parser(frameFullpath, xyPixelSize, dz, varargin)
% Deskew and/or rotate data for a single zarr file that cannot be fitted to
% memory
% 
%
% Author: Xiongtao Ruan (02/16/2022)
%
% Based on XR_deskewRotateFrame.m
%
% xruan (12/14/2022): add support for input bbox, that is, crop the data before dsr


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('xyPixelSize'); 
ip.addRequired('dz'); 
ip.addParameter('ObjectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Crop', false, @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('Reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Save16bit', false , @(x) islogical(x) || ischar(x)); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('SaveMIP', true , @(x) islogical(x) || ischar(x)); % save MIP-z for ds and dsr. 
ip.addParameter('saveZarr', false , @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('BatchSize', [1024, 1024, 1024] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('zarrSubSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) (isempty(x) || isvector(x)) || ischar(x));
ip.addParameter('taskSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('DSRCombined', true, @(x) islogical(x) || ischar(x)); % combined processing 
ip.addParameter('resample', [], @(x) (isempty(x) || isnumeric(x)) || ischar(x)); % resampling after rotation 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})) && ischar(x));
ip.addParameter('maskFns', {}, @(x) iscell(x) || ischar(x)); % 2d masks to filter regions to deskew and rotate, in xy, xz, yz orde
ip.addParameter('surffix', '', @ischar); % suffix for the folder
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
Crop = pr.Crop;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
Overwrite = pr.Overwrite;
ObjectiveScan = pr.ObjectiveScan;
flipZstack = pr.flipZstack;
Save16bit = pr.Save16bit;
SaveMIP = pr.SaveMIP;
DSRCombined = pr.DSRCombined;
resample = pr.resample;
saveZarr = pr.saveZarr;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
zarrSubSize = pr.zarrSubSize;
inputBbox = pr.inputBbox;
taskSize = pr.taskSize;
Interp = pr.Interp;
maskFns = pr.maskFns;
surffix = pr.surffix;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;

uuid = pr.uuid;
debug = pr.debug;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(frameFullpath)
    frameFullpath = eval(frameFullpath);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(ObjectiveScan)
    ObjectiveScan = strcmp(ObjectiveScan,'true');
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite,'true');
end
if ischar(Crop)
    Crop = strcmp(Crop,'true');
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(Reverse)
    Reverse = strcmp(Reverse,'true');
end
if ischar(flipZstack)
    flipZstack = strcmp(flipZstack,'true');
end
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit,'true');
end
if ischar(SaveMIP)
    SaveMIP = strcmp(SaveMIP,'true');
end
if ischar(saveZarr)
    saveZarr = strcmp(saveZarr,'true');
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
if ischar(DSRCombined)
    DSRCombined = strcmp(DSRCombined,'true');
end
if ischar(resample)
    resample = str2num(resample);
end
if ischar(maskFns)
    maskFns = eval(maskFns);
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

XR_deskewRotateZarr(frameFullpath, xyPixelSize, dz, 'ObjectiveScan',ObjectiveScan,...
    'Overwrite',Overwrite,'Crop',Crop,'SkewAngle',SkewAngle,'Reverse',Reverse,...
    'flipZstack',flipZstack,'Save16bit',Save16bit,'SaveMIP',SaveMIP,'saveZarr',saveZarr,...
    'BatchSize',BatchSize,'BlockSize',BlockSize,'inputBbox',inputBbox,'zarrSubSize',zarrSubSize, ...
    'taskSize',taskSize,'DSRCombined',DSRCombined,'resample',resample,'Interp',Interp, ...
    'maskFns',maskFns,'surffix',surffix,'parseCluster',parseCluster,'parseParfor',parseParfor, ...
    'masterCompute',masterCompute,'jobLogDir',jobLogDir,'cpusPerTask',cpusPerTask,...
    'uuid',uuid,'debug',debug,'mccMode',mccMode,'ConfigFile',ConfigFile);

end
