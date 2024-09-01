function [] = tiffToZarr_parser(tifFilename, zarrFilename, frame, varargin)


%#function processFFCorrectionFrame
%#function erodeVolumeBy2DProjection

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tifFilename', @(x) ischar(x) || iscell(x));
ip.addRequired('zarrFilename', @ischar);
ip.addOptional('frame', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('blockSize', [500, 500, 250], @(x) isnumeric(x) || ischar(x));
ip.addParameter('shardSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('expand2dDim', true, @(x) islogical(x) || ischar(x)); % expand the z dimension for 2d data
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('resampleFactor', [], @(x) isempty(x) || isnumeric(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('tileOutBbox', [], @(x) isempty(x) || isnumeric(x) || ischar(x));
ip.addParameter('readWholeTiff', true, @(x) islogical(x) || ischar(x));
ip.addParameter('compressor', 'zstd', @ischar);
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x) || isstring(x));
ip.addParameter('uuid', '', @ischar);

ip.parse(tifFilename, zarrFilename, frame, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
blockSize = pr.blockSize;
shardSize = pr.shardSize;
expand2dDim = pr.expand2dDim;
flipZstack = pr.flipZstack;
resampleFactor = pr.resampleFactor;
inputBbox = pr.inputBbox;
tileOutBbox = pr.tileOutBbox;
readWholeTiff = pr.readWholeTiff;
compressor = pr.compressor;
usrFcn = pr.usrFcn;
uuid = pr.uuid;

if ischar(tifFilename) && ~isempty(tifFilename) && strcmp(tifFilename(1), '{')
    tifFilename = eval(tifFilename);
end
if ischar(frame)
    frame = str2num(frame);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(shardSize)
    shardSize = str2num(shardSize);
end
if ischar(expand2dDim)
    expand2dDim = str2num(expand2dDim);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(tileOutBbox)
    tileOutBbox = str2num(tileOutBbox);
end
if ischar(readWholeTiff)
    readWholeTiff = str2num(readWholeTiff);
end
if ischar(usrFcn) && ~isempty(usrFcn) && (strcmp(usrFcn(1), '{') || strcmp(usrFcn(1), '[') || strcmp(usrFcn(1), '@'))
    usrFcn = eval(usrFcn);
end

tiffToZarr(tifFilename, zarrFilename, frame, Overwrite=Overwrite, blockSize=blockSize, ...
    shardSize=shardSize, expand2dDim=expand2dDim, flipZstack=flipZstack, resampleFactor=resampleFactor, ...
    inputBbox=inputBbox, tileOutBbox=tileOutBbox, readWholeTiff=readWholeTiff, ...
    compressor=compressor, usrFcn=usrFcn, uuid=uuid);

end

