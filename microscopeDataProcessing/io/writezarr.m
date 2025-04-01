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
    options.shardSize (1, :) {mustBeNumeric} = []
    options.zarrSubSize (1, :) {mustBeNumeric} = []
    options.expand2dDim (1, 1) logical = true  % expand the z dimension for 2d data
    options.groupWrite (1, 1) logical = true
    options.bbox (1, :) {mustBeNumeric} = []
    options.sparseData (1, :) logical = false
    options.dimSeparator char = '.'
    options.create logical = true
end

overwrite = options.overwrite;
expand2dDim = options.expand2dDim;
blockSize = options.blockSize;
shardSize = options.shardSize;
zarrSubSize = options.zarrSubSize;
groupWrite = options.groupWrite;
bbox = options.bbox;
sparseData = options.sparseData;
dimSeparator = options.dimSeparator;
create = options.create;

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

newFile = false;
try 
    if ~exist(filepath, 'dir')
        if create
            if exist(filepath, 'dir')
                rmdir(filepath, 's');
            end
            createzarr(filepath, dataSize=sz, blockSize=blockSize, shardSize=shardSize, ...
                dtype=dtype, zarrSubSize=zarrSubSize, dimSeparator=dimSeparator);
            newFile = true;
        else
            error('The zarr file %s does not exist!', filepath);
        end
    end

    parallelWriteZarr(filepath, data, 'bbox', bbox, 'sparse', sparseData);
catch ME
    pe = pyenv;    
    if pe.Status ~= "Loaded"
        throw(ME);
    end

    disp(ME);
    disp('Use the alternative zarr writer (ZarrAdapter)...');    
    if ~exist(filepath, 'dir') || (exist(filepath, 'dir') && isempty(bbox))
        if create
            if exist(filepath, 'dir')
                rmdir(filepath, 's');
            end
            bim = blockedImage(filepath, sz, blockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');
            newFile = true;
        else
            error('The zarr file %s does not exist!', filepath);
        end
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

if isunix && groupWrite && newFile
    try
        fileattrib(filepath, '+w', 'g');
    catch
        warning('Unable to change file attribe for group write!');
    end
end

end
