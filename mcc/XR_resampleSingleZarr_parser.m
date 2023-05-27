function [] = XR_resampleSingleZarr_parser(zarrFullpath, dsFullpath, dsFactor, varargin)
% 
% Author: Xiongtao Ruan (12/14/2020)
% xruan (11/09/2021): enable arbitray blockSize

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath'); 
ip.addRequired('dsFullpath'); 
ip.addRequired('dsFactors'); 
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric || ischar(x)); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @(x) isnumeric || ischar(x)); % size to process in one batch
ip.addParameter('zarrSubSize', [20, 20, 20], @(x) isnumeric || ischar(x));
ip.addParameter('BorderSize', [5, 5, 5], @(x) isnumeric|| ischar(x)); % padded boarder for each batch
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear', 'nearest'})) && ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical || ischar(x));
ip.addParameter('cpusPerTask', 1, @(x) islogical || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(zarrFullpath, dsFullpath, dsFactor, varargin{:});

pr = ip.Results;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
BorderSize = pr.BorderSize;
Interp = pr.Interp;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dsFactor)
    dsFactor = str2num(dsFactor);
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
    parseCluster = strcmp(parseCluster,'true');
end
if ischar(cpusPerTask)
    cpusPerTask = strcmp(cpusPerTask,'true');
end
if ischar(mccMode)
    mccMode = strcmp(mccMode,'true');
end

XR_resampleSingleZarr(zarrFullpath,dsFullpath,dsFactor,'blockSize',blockSize,...
    'batchSize',batchSize,'zarrSubSize',zarrSubSize,'BorderSize',BorderSize,...
    'Interp',Interp,'parseCluster',parseCluster,'cpusPerTask',cpusPerTask,...
    'uuid',uuid, mccMode=mccMode,ConfigFile=ConfigFile);

end
