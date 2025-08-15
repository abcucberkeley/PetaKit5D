function [] = XR_bioformats_to_tiff_or_zarr_wrapper_parser(dataPaths, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'tiffs', @ischar);
ip.addParameter('channelPatterns', {'.nd2', '.czi'}, @(x) iscell(x) || ischar(x));
ip.addParameter('nChannels', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dataFormats', {'.nd2', '.czi'}, @(x) iscell(x) || ischar(x));
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x));
ip.addParameter('blockSize', [256, 256, 256], @(x) isnumeric(x) || ischar(x));
ip.addParameter('overWrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
channelPatterns = pr.channelPatterns;
nChannels = pr.nChannels;
dataFormats = pr.dataFormats;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
overWrite = pr.overWrite;
uuid = pr.uuid;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(nChannels)
    nChannels = str2num(nChannels);
end
if ischar(dataFormats) && ~isempty(dataFormats) && strcmp(dataFormats(1), '{')
    dataFormats = eval(dataFormats);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(overWrite)
    overWrite = str2num(overWrite);
end

XR_bioformats_to_tiff_or_zarr_wrapper(dataPaths, resultDirName=resultDirName, ...
    channelPatterns=channelPatterns, nChannels=nChannels, dataFormats=dataFormats, ...
    saveZarr=saveZarr, blockSize=blockSize, overWrite=overWrite, uuid=uuid);

end

