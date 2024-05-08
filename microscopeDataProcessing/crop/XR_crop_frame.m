function [] = XR_crop_frame(dataFullpath, saveFullpath, bbox, varargin)
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
ip.addRequired('bbox', @isnumeric);
ip.addParameter('overwrite', false, @islogical); % start coordinate of the last time point
ip.addParameter('pad', false, @islogical); % pad region that is outside the bbox
ip.addParameter('zarrFile', false , @islogical); % read zarr
ip.addParameter('largeZarr', false, @islogical); % use zarr file as input
ip.addParameter('saveZarr', false , @islogical); % save as zarr
ip.addParameter('blockSize', [500, 500, 500] , @isnumeric); % save as zarr
ip.addParameter('uuid', '', @ischar);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataFullpath, saveFullpath, bbox, varargin{:});

pr = ip.Results;
overwrite = pr.overwrite;
pad = pr.pad;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
uuid = pr.uuid;
parseCluster = pr.parseCluster;
mccMode = pr.mccMode;
configFile = pr.configFile;

if isempty(uuid)
    uuid = get_uuid();
end

if ~exist(dataFullpath, 'file')
    warning('The file %s does not exist!', dataFullpath);
    return;
end

if exist(saveFullpath, 'file') && ~overwrite
    fprintf('The cropped file %s is already exist, skip it!\n', saveFullpath);
    return;
end

fprintf('Crop frame %s with bounding box [%s]...\n', dataFullpath, num2str(bbox));

imSize = getImageSize(dataFullpath);
if pad && any(bbox(1 : 3) < 1 | bbox(4 : 6) > imSize)
    bbox_1 = [max(1, bbox(1 : 3)), min(imSize, bbox(4 : 6))];
else
    bbox_1 = bbox;
    bbox_1(4 : 6) = min(bbox_1(4 : 6), imSize);
    % only read cropped slices. 
end

% read data
if zarrFile
    if largeZarr
        batchSize = min([1024, 1024, 1024], blockSize * 4);
        parseParfor = false;
        XR_crop_zarr(dataFullpath, saveFullpath, bbox, 'pad', pad, 'batchSize', batchSize, ...
            'blockSize', blockSize, 'uuid', uuid, 'parseCluster', parseCluster, ...
            'parseParfor', parseParfor, mccMode=mccMode, configFile=configFile);
        return;
    end
    im = readzarr(dataFullpath, 'inputBbox', bbox_1);
else
    im = readtiff(dataFullpath, 'range', [bbox_1(3), bbox_1(6)]);
    try
        im = crop3d_mex(im, [bbox_1(1 : 2), 1, bbox_1(4 : 5), size(im, 3)]);
    catch ME
        disp(ME)
        im = im(bbox_1(1) : bbox_1(4), bbox_1(2) : bbox_1(5), :);
    end
end

% pad cropped data
if pad 
    if any(bbox(1 : 3) < 1)
        im = padarray(im, max(0, 1 - bbox(1 : 3)), 0, 'pre');
    end
    if any(bbox(4 : 6) > imSize)
        im = padarray(im, max(0, bbox(4 : 6) - imSize), 0, 'post');        
    end
end

% save data
if saveZarr
    tmpPath = sprintf('%s_%s.zarr', saveFullpath(1 : end - 5), uuid);
    writezarr(im, tmpPath, 'blockSize', blockSize);    
else
    tmpPath = sprintf('%s_%s.tif', saveFullpath(1 : end - 4), uuid);
    writetiff(im, tmpPath);
end
movefile(tmpPath, saveFullpath);

% fprintf('Done!\n\n');

end
