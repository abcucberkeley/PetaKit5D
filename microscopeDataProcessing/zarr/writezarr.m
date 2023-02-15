function writezarr(data, filepath, options)
% wrapper for zarr writer 
% 
% Author: Xiongtao Ruan (01/25/2022)
% 
% xruan (05/23/2022): add support for bounding box write; also add
% parallelWriteZarr as default method


arguments
    data {mustBeNumeric}
    filepath char
    options.overwrite (1, 1) logical = false
    options.blockSize (1, :) {mustBeNumeric} = [500, 500, 500]
    options.expand2dDim (1, 1) logical = true  % expand the z dimension for 2d data
    options.groupWrite (1, 1) logical = true
    options.bbox (1, :) {mustBeNumeric} = []
end

overwrite = options.overwrite;
expand2dDim = options.expand2dDim;
blockSize = options.blockSize;
groupWrite = options.groupWrite;
bbox = options.bbox;

dtype = class(data);
sz = size(data);
init_val = zeros(1, dtype);
if ismatrix(data) 
    if expand2dDim
        sz(3) = 1;
        blockSize(3) = 1;
    else
        blockSize = blockSize(1 : 2);
    end
end
blockSize = min(sz, blockSize);
% overwrite the zarr file
if overwrite && exist(filepath, 'dir') && isempty(bbox)
    rmdir(filepath, 's');
end

% if not overwrite data, the file exist and bbox is empty, define bbox as
% the whole data size. 
if isempty(bbox)
    bbox = [1, 1, 1, size(data, 1 : 3)];
    if ~expand2dDim
        bbox = [1, 1, size(data, 1 : 2)];
    end
end

try 
    if ~exist(filepath, 'dir') || (exist(filepath, 'dir') && isempty(bbox))
        if exist(filepath, 'dir')
            rmdir(filepath, 's');
        end
        bim = blockedImage(filepath, sz, blockSize, init_val, "Adapter", CZarrAdapter, 'Mode', 'w');
        bim.Adapter.close();
    else
        bim = blockedImage(filepath, "Adapter", CZarrAdapter);
        bim.Adapter.close();
    end

    parallelWriteZarr(filepath, data, 1, bbox);
catch ME
    disp(ME);
    disp('Use the alternative zarr writer (ZarrAdapter)...');    
    if ~exist(filepath, 'dir') || (exist(filepath, 'dir') && isempty(bbox))
        if exist(filepath, 'dir')
            rmdir(filepath, 's');
        end
        bim = blockedImage(filepath, sz, blockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');             
    else
        bim = blockedImage(filepath, "Adapter", ZarrAdapter);
    end

    % for data greater than 2GB, use multiprocessing
    if isempty(bbox)
        if ~ispc && prod(sz) * 2 / 1024^3 > 2 && ~ismatrix(data) 
            bim.Adapter.setData(data);
        else
            bim.Adapter.setRegion(ones(1, numel(bim.Size)), bim.Size, data)
        end
    else
        bim.Adapter.setRegion(bbox(1 : 3), bbox(4 : 6), data)
    end
    bim.Adapter.close();
end

if groupWrite
    try
        fileattrib(filepath, '+w', 'g');
    catch
        warning('Unable to change file attribe for group write!');
    end
end

end
