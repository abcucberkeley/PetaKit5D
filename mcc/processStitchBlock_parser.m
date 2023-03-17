function [] = processStitchBlock_parser(blockInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchBlockInfo, zarrHeaders, nv_bim, varargin)

%#function processStitchBlock

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('BlockInfoFullname', @(x) ischar(x));
ip.addRequired('PerBlockInfoFullname', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addOptional('stitchBlockInfo', [], @(x) isnumeric(x) || ischar(x));
ip.addOptional('zarrHeaders', [], @(x) isnumeric(x) || ischar(x));
ip.addOptional('nv_bim', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('BlendMethod', 'mean', @ischar);
ip.addParameter('BorderSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BlurSigma', 5, @(x) isnumeric(x) || ischar(x)); % blurred sigma for blurred blend
ip.addParameter('imdistFullpaths', {}, @(x) iscell(x) || ischar(x)); % image distance paths
ip.addParameter('weightDegree', 10, @(x) isnumeric(x) || ischar(x)); % weight degree for image distances

ip.parse(blockInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchBlockInfo, zarrHeaders, nv_bim, varargin{:});

Overwrite = ip.Results.Overwrite;
BlendMethod = ip.Results.BlendMethod;
BorderSize = ip.Results.BorderSize;
BlurSigma = ip.Results.BlurSigma;
imdistFullpaths = ip.Results.imdistFullpaths;
weightDegree = ip.Results.weightDegree;

if ischar(blockInds)
    blockInds = str2num(blockInds);
end
if ischar(stitchBlockInfo)
    stitchBlockInfo = str2num(stitchBlockInfo);
end
if ischar(zarrHeaders)
    zarrHeaders = str2num(zarrHeaders);
end
if ischar(nv_bim)
    nv_bim = str2num(nv_bim);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
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
    flagFullname, stitchBlockInfo, zarrHeaders, nv_bim, Overwrite=Overwrite, ...
    BlendMethod=BlendMethod, BorderSize=BorderSize, BlurSigma=BlurSigma, ...
    imdistFullpaths=imdistFullpaths, weightDegree=weightDegree);

end

