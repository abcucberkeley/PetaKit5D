function [] = XR_crop_frame_parser(dataFullpath, saveFullpath, bbox, varargin)
% crop a frame 
% bbox: ymin, xmin, zmin, ymax, xmax, zmax
% If ymax, xmax or zmax is larger than image size, use image size as upper bounds. 
% 
% Author: Xiongtao Ruan (03/11/2020)
% 
% xruan (08/22/2020): update function for writing results
% xruan (07/13/2021): add option to pad data if it is outside of the bbox
% xruan (01/25/2022): add support for zarr read and write
% xruan (06/03/2022): add support for large zarr files


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
    overwrite = strcmp(overwrite,'true');
end
if ischar(pad)
    pad = strcmp(pad,'true');
end
if ischar(zarrFile)
    zarrFile = strcmp(zarrFile,'true');
end
if ischar(largeZarr)
    largeZarr = strcmp(largeZarr,'true');
end
if ischar(saveZarr)
    saveZarr = strcmp(saveZarr,'true');
end
if ischar(BlockSize)
    BlockSize = str2num(BlockSize);
end
if ischar(mccMode)
    mccMode = strcmp(mccMode, 'true');
end

XR_crop_frame(dataFullpath,saveFullpath,bbox,'overwrite',overwrite,'pad',pad, ...
    'zarrFile',zarrFile,'largeZarr',largeZarr,'saveZarr',saveZarr,'BlockSize',BlockSize, ...
    mccMode=mccMode, ConfigFile=ConfigFile);

end
