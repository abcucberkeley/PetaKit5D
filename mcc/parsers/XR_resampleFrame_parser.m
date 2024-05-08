function [] = XR_resampleFrame_parser(fn, fnout, resampleFactor, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fn', @ischar);
ip.addRequired('fnout', @ischar);
ip.addRequired('resampleFactor', @(x) isnumeric(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x)); % bbox for input
ip.addParameter('interpMethod', 'linear', @ischar);
ip.addParameter('save16bit', true ,@(x) islogical(x) || ischar(x)); % saves 16bit, else single
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); % use zarr file as output
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x)); % blcoksize
ip.addParameter('uuid', '', @ischar);

ip.parse(fn, fnout, resampleFactor, varargin{:});

pr = ip.Results;
inputBbox = pr.inputBbox;
interpMethod = pr.interpMethod;
save16bit = pr.save16bit;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
uuid = pr.uuid;

if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
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

XR_resampleFrame(fn, fnout, resampleFactor, inputBbox=inputBbox, interpMethod=interpMethod, ...
    save16bit=save16bit, zarrFile=zarrFile, saveZarr=saveZarr, blockSize=blockSize, ...
    uuid=uuid);

end

