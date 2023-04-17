function [] = processStitchBlock_parser(blockInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchFullname, stitchBlockInfo, tileFns, varargin)

%#function processStitchBlock

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('BlockInfoFullname', @(x) ischar(x));
ip.addRequired('PerBlockInfoFullname', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('stitchFullnanme', @(x) ischar(x));
ip.addOptional('stitchBlockInfo', [], @(x) isnumeric(x) || ischar(x));
ip.addOptional('tileFns', [], @(x) iscell(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('imSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('blockSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('dtype', [], @(x) ischar(x));
ip.addParameter('BlendMethod', 'mean', @ischar);
ip.addParameter('BorderSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BlurSigma', 5, @(x) isnumeric(x) || ischar(x)); % blurred sigma for blurred blend
ip.addParameter('imdistFullpaths', {}, @(x) iscell(x) || ischar(x)); % image distance paths
ip.addParameter('weightDegree', 10, @(x) isnumeric(x) || ischar(x)); % weight degree for image distances

ip.parse(blockInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchFullname, stitchBlockInfo, tileFns, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
imSize = pr.imSize;
blockSize = pr.blockSize;
dtype = pr.dtype;
BlendMethod = pr.BlendMethod;
BorderSize = pr.BorderSize;
BlurSigma = pr.BlurSigma;
imdistFullpaths = pr.imdistFullpaths;
weightDegree = pr.weightDegree;

if ischar(blockInds)
    blockInds = str2num(blockInds);
end
if ischar(stitchBlockInfo)
    stitchBlockInfo = str2num(stitchBlockInfo);
end
if ischar(tileFns)
    tileFns = eval(tileFns);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
end
if ischar(imSize)
    imSize = str2num(imSize);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(BorderSize)
    BorderSize = str2num(BorderSize);
end
if ischar(BlurSigma)
    BlurSigma = str2double(BlurSigma);
end
if ischar(imdistFullpaths)
    imdistFullpaths = eval(imdistFullpaths);
end
if ischar(weightDegree)
    weightDegree = str2double(weightDegree);
end

processStitchBlock(blockInds, BlockInfoFullname, PerBlockInfoFullname, ...
    flagFullname, stitchFullname, stitchBlockInfo, tileFns, Overwrite=Overwrite, ...
    imSize=imSize, blockSize=blockSize, dtype=dtype, BlendMethod=BlendMethod, ...
    BorderSize=BorderSize, BlurSigma=BlurSigma,imdistFullpaths=imdistFullpaths, ...
    weightDegree=weightDegree);

end

