function [] = processStitchBlock_parser(batchInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchFullname, stitchBlockInfo, tileFns, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('BlockInfoFullname', @(x) ischar(x));
ip.addRequired('PerBlockInfoFullname', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('stitchFullname', @(x) ischar(x));
ip.addOptional('stitchBlockInfo', [], @(x) isempty(x) || isstruct(x) || ischar(x));
ip.addOptional('tileFns', [], @(x) isempty(x) || iscell(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('imSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('batchSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('dtype', 'uint16', @ischar);
ip.addParameter('BlendMethod', 'feather', @ischar);
ip.addParameter('BorderSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('imdistFullpaths', {}, @(x) iscell(x) || ischar(x)); % image distance paths
ip.addParameter('imdistFileIdx', [], @(x) isnumeric(x) || ischar(x)); % image distance paths indices
ip.addParameter('poolSize', [], @(x) isnumeric(x) || ischar(x)); % distance matrix with max pooling factors
ip.addParameter('weightDegree', 10, @(x) isnumeric(x) || ischar(x)); % weight degree for image distances

ip.parse(batchInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchFullname, stitchBlockInfo, tileFns, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
imSize = pr.imSize;
batchSize = pr.batchSize;
dtype = pr.dtype;
BlendMethod = pr.BlendMethod;
BorderSize = pr.BorderSize;
imdistFullpaths = pr.imdistFullpaths;
imdistFileIdx = pr.imdistFileIdx;
poolSize = pr.poolSize;
weightDegree = pr.weightDegree;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(stitchBlockInfo) && ~isempty(stitchBlockInfo) && (strcmp(stitchBlockInfo(1), '{') || strcmp(stitchBlockInfo(1), '[') || strcmp(stitchBlockInfo(1), '@'))
    stitchBlockInfo = eval(stitchBlockInfo);
end
if ischar(tileFns) && ~isempty(tileFns) && (strcmp(tileFns(1), '{') || strcmp(tileFns(1), '[') || strcmp(tileFns(1), '@'))
    tileFns = eval(tileFns);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(imSize)
    imSize = str2num(imSize);
end
if ischar(batchSize)
    batchSize = str2num(batchSize);
end
if ischar(BorderSize)
    BorderSize = str2num(BorderSize);
end
if ischar(imdistFullpaths) && ~isempty(imdistFullpaths) && strcmp(imdistFullpaths(1), '{')
    imdistFullpaths = eval(imdistFullpaths);
end
if ischar(imdistFileIdx)
    imdistFileIdx = str2num(imdistFileIdx);
end
if ischar(poolSize)
    poolSize = str2num(poolSize);
end
if ischar(weightDegree)
    weightDegree = str2num(weightDegree);
end

processStitchBlock(batchInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, ...
    stitchFullname, stitchBlockInfo, tileFns, Overwrite=Overwrite, imSize=imSize, ...
    batchSize=batchSize, dtype=dtype, BlendMethod=BlendMethod, BorderSize=BorderSize, ...
    imdistFullpaths=imdistFullpaths, imdistFileIdx=imdistFileIdx, poolSize=poolSize, ...
    weightDegree=weightDegree);

end

