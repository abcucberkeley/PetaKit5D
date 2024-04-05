function [] = XR_crop_frame_parser(dataFullpath, saveFullpath, bbox, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataFullpath', @(x) ischar(x) || iscell(x));
ip.addRequired('saveFullpath', @(x) ischar(x) || iscell(x));
ip.addRequired('bbox', @(x) isnumeric(x) || ischar(x));
ip.addParameter('overwrite', false, @(x) islogical(x) || ischar(x)); % start coordinate of the last time point
ip.addParameter('pad', false, @(x) islogical(x) || ischar(x)); % pad region that is outside the bbox
ip.addParameter('zarrFile', false , @(x) islogical(x) || ischar(x)); % read zarr
ip.addParameter('largeZarr', false, @(x) islogical(x) || ischar(x)); % use zarr file as input
ip.addParameter('saveZarr', false , @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('BlockSize', [500, 500, 500] , @(x) isnumeric(x) || ischar(x)); % save as zarr
ip.addParameter('uuid', '', @ischar);
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataFullpath, saveFullpath, bbox, varargin{:});

pr = ip.Results;
overwrite = pr.overwrite;
pad = pr.pad;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
saveZarr = pr.saveZarr;
BlockSize = pr.BlockSize;
uuid = pr.uuid;
parseCluster = pr.parseCluster;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataFullpath) && strcmp(dataFullpath(1), '{')
    dataFullpath = eval(dataFullpath);
end
if ischar(saveFullpath) && strcmp(saveFullpath(1), '{')
    saveFullpath = eval(saveFullpath);
end
if ischar(bbox)
    bbox = str2num(bbox);
end
if ischar(overwrite)
    overwrite = str2num(overwrite);
end
if ischar(pad)
    pad = str2num(pad);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(largeZarr)
    largeZarr = str2num(largeZarr);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(BlockSize)
    BlockSize = str2num(BlockSize);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_crop_frame(dataFullpath, saveFullpath, bbox, overwrite=overwrite, pad=pad, ...
    zarrFile=zarrFile, largeZarr=largeZarr, saveZarr=saveZarr, BlockSize=BlockSize, ...
    uuid=uuid, parseCluster=parseCluster, mccMode=mccMode, ConfigFile=ConfigFile);

end

