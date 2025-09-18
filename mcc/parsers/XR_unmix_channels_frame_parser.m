function XR_unmix_channels_frame_parser(frameFullpaths, unmixFactors, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('unmixFactors', @(x) isnumeric(x) || ischar(x));
ip.addParameter('mode', 'linear', @ischar); % linear vs gaussian
ip.addParameter('unmixSigmas', [], @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('resultDirName', 'Unmixed', @ischar); 
ip.addParameter('channelInd', 1, @(x) isnumeric(x) || ischar(x)); % unmix for which channel
ip.addParameter('save16bit', true, @(x) islogical(x) || ischar(x)); 
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('blockSize', [256, 256, 256] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(frameFullpaths, unmixFactors, varargin{:});

pr = ip.Results;
mode = pr.mode;
unmixSigmas = pr.unmixSigmas;
resultDirName = pr.resultDirName;
channelInd = pr.channelInd;
save16bit = pr.save16bit;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
uuid = pr.uuid;
debug = pr.debug;

if ischar(frameFullpaths) && ~isempty(frameFullpaths) && strcmp(frameFullpaths(1), '{')
    frameFullpaths = eval(frameFullpaths);
end
if ischar(unmixFactors)
    unmixFactors = str2num(unmixFactors);
end
if ischar(unmixSigmas)
    unmixSigmas = str2num(unmixSigmas);
end
if ischar(channelInd)
    channelInd = str2num(channelInd);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_unmix_channels_frame(frameFullpaths, unmixFactors, mode=mode, unmixSigmas=unmixSigmas, ...
    resultDirName=resultDirName, channelInd=channelInd, save16bit=save16bit, ...
    zarrFile=zarrFile, saveZarr=saveZarr, blockSize=blockSize, uuid=uuid, debug=debug);

end

