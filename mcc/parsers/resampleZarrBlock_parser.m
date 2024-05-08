function [] = resampleZarrBlock_parser(batchInds, zarrFullpath, dsFullpath, flagFullname, resampleFactor, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('dsFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('resampleFactor', @(x) isnumeric(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('batchSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('blockSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('borderSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('interpMethod', 'linear', @ischar);

ip.parse(batchInds, zarrFullpath, dsFullpath, flagFullname, resampleFactor, varargin{:});

pr = ip.Results;
inputBbox = pr.inputBbox;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
borderSize = pr.borderSize;
overwrite = pr.overwrite;
interpMethod = pr.interpMethod;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(borderSize)
    borderSize = str2num(borderSize);
end
if ischar(overwrite)
    overwrite = str2num(overwrite);
end

resampleZarrBlock(batchInds, zarrFullpath, dsFullpath, flagFullname, resampleFactor, ...
    inputBbox=inputBbox, batchSize=batchSize, blockSize=blockSize, borderSize=borderSize, ...
    overwrite=overwrite, interpMethod=interpMethod);

end

