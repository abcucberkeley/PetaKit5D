function [] = XR_resampleSingleZarr_parser(zarrFullpath, rsFullpath, rsFactor, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath', @ischar); 
ip.addRequired('rsFullpath', @ischar); 
ip.addRequired('rsFactor', @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('bbox', [], @(x) isnumeric(x) || ischar(x)); % bbox for input
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x)); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @(x) isnumeric(x) || ischar(x)); % size to process in one batch 
ip.addParameter('zarrSubSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BorderSize', [5, 5, 5], @(x) isnumeric(x) || ischar(x)); % padded boarder for each batch
ip.addParameter('Interp', 'linear', @(x) ischar(x) && any(strcmpi(x, {'cubic', 'linear', 'nearest'})));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('cpusPerTask', 1, @(x) isscalar(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(zarrFullpath, rsFullpath, rsFactor, varargin{:});

pr = ip.Results;
bbox = pr.bbox;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
zarrSubSize = pr.zarrSubSize;
BorderSize = pr.BorderSize;
Interp = pr.Interp;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(rsFactor)
    rsFactor = str2num(rsFactor);
end
if ischar(bbox)
    bbox = str2num(bbox);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(zarrSubSize)
    zarrSubSize = str2num(zarrSubSize);
end
if ischar(BorderSize)
    BorderSize = str2num(BorderSize);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_resampleSingleZarr(zarrFullpath, rsFullpath, rsFactor, bbox=bbox, blockSize=blockSize, ...
    batchSize=batchSize, zarrSubSize=zarrSubSize, BorderSize=BorderSize, Interp=Interp, ...
    parseCluster=parseCluster, cpusPerTask=cpusPerTask, uuid=uuid, mccMode=mccMode, ...
    ConfigFile=ConfigFile);

end

