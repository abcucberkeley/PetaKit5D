function createzarr(filepath, options)
% wrapper for initializing zarr file
% 
% Author: Xiongtao Ruan (03/12/2023)


arguments
    filepath char
    options.overwrite (1, 1) logical = false
    options.dataSize (1, :) {mustBeNumeric} = [500, 500, 500]
    options.blockSize (1, :) {mustBeNumeric} = [500, 500, 500]
    options.dtype char = 'uint16'
    options.order char = 'F'
    options.expand2dDim (1, 1) logical = true  % expand the z dimension for 2d data
    options.groupWrite (1, 1) logical = true
    options.compressor char = 'zstd'
    options.zarrSubSize (1, :) {mustBeNumeric} = []
end

dataSize = options.dataSize;
blockSize = options.blockSize;
dtype = options.dtype;
order = options.order;
expand2dDim = options.expand2dDim;
groupWrite = options.groupWrite;
compressor = options.compressor;
zarrSubSize = options.zarrSubSize;

if numel(dataSize) == 2 
    if expand2dDim
        dataSize(3) = 1;
        blockSize(3) = 1;
    else
        blockSize = blockSize(1 : 2);
    end
end
blockSize = min(dataSize, blockSize);
% overwrite the zarr file
if exist(filepath, 'dir')
    rmdir(filepath, 's');
end

switch dtype
    case 'single'
        ddtype = 'f4';
    case 'uint16'
        ddtype = 'u2';
    otherwise
        error('Unsupported data type');
end

try
    if ~isempty(zarrSubSize) && prod(ceil(dataSize ./ blockSize)) < prod(zarrSubSize)
        zarrSubSize = [];
    end
    if isempty(zarrSubSize)
        createZarrFile(filepath, 'chunks', blockSize, 'dtype', ddtype, 'order', order, ...
            'shape', dataSize, 'cname', compressor, 'level', 1); 
    else
        createZarrFile(filepath, 'chunks', blockSize, 'dtype', ddtype, 'order', order, ...
            'shape', dataSize, 'cname', compressor, 'level', 1, 'subfolders', zarrSubSize); 
    end
catch ME
    disp(ME);
    disp('Use the alternative zarr initializer (ZarrAdapter)...');    
    init_val = zeros(1, dtype);
    bim = blockedImage(filepath, dataSize, blockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');             
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
