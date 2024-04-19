function [] = XR_resample_dataset_parser(dataPaths, rsfactor, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('rsfactor', @(x) isnumeric(x) || ischar(x));
ip.addParameter('outDirStr', 'resampled', @ischar);
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @(x) iscell(x) || ischar(x));
ip.addParameter('bbox', [], @(x) isnumeric(x) || ischar(x)); % bbox for input
ip.addParameter('Interp', 'linear', @(x) ischar(x) && any(strcmpi(x, {'cubic', 'linear', 'nearest'})));
ip.addParameter('Save16bit', true, @(x) islogical(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('largeZarr', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); % use zarr file as output
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x)); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @(x) isnumeric(x) || ischar(x)); % size to process in one batch
ip.addParameter('BorderSize', [5, 5, 5], @(x) isnumeric(x) || ischar(x)); % padded boarder for each batch
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 2, @(x) isscalar(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPaths, rsfactor, varargin{:});

pr = ip.Results;
outDirStr = pr.outDirStr;
ChannelPatterns = pr.ChannelPatterns;
bbox = pr.bbox;
Interp = pr.Interp;
Save16bit = pr.Save16bit;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
BorderSize = pr.BorderSize;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(rsfactor)
    rsfactor = str2num(rsfactor);
end
if ischar(ChannelPatterns) && ~isempty(ChannelPatterns) && strcmp(ChannelPatterns(1), '{')
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(bbox)
    bbox = str2num(bbox);
end
if ischar(Save16bit)
    Save16bit = str2num(Save16bit);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(largeZarr)
    largeZarr = str2num(largeZarr);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(BorderSize)
    BorderSize = str2num(BorderSize);
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

XR_resample_dataset(dataPaths, rsfactor, outDirStr=outDirStr, ChannelPatterns=ChannelPatterns, ...
    bbox=bbox, Interp=Interp, Save16bit=Save16bit, zarrFile=zarrFile, largeZarr=largeZarr, ...
    saveZarr=saveZarr, blockSize=blockSize, batchSize=batchSize, BorderSize=BorderSize, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, uuid=uuid, mccMode=mccMode, ConfigFile=ConfigFile);

end

