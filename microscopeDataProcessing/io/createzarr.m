function createzarr(filepath, options)
% wrapper for initializing zarr file
% 
% Author: Xiongtao Ruan (03/12/2023)


arguments
    filepath char
    options.overwrite (1, 1) logical = false
    options.dataSize (1, :) {mustBeNumeric} = [500, 500, 500]
    options.blockSize (1, :) {mustBeNumeric} = [500, 500, 500]
    options.shardSize (1, :) {mustBeNumeric} = []    
    options.dtype char = 'uint16'
    options.order char = 'F'
    options.expand2dDim (1, 1) logical = true  % expand the z dimension for 2d data
    options.groupWrite (1, 1) logical = true
    options.compressor char = 'zstd'
    options.clevel (1, 1) {mustBeNumeric} = 1
    options.zarrSubSize (1, :) {mustBeNumeric} = []
    options.dimSeparator char = '.'
end

dataSize = options.dataSize;
blockSize = options.blockSize;
shardSize = options.shardSize;
dtype = options.dtype;
order = options.order;
expand2dDim = options.expand2dDim;
groupWrite = options.groupWrite;
compressor = options.compressor;
clevel = options.clevel;
zarrSubSize = options.zarrSubSize;
dimSeparator = options.dimSeparator;

if numel(dataSize) == 2 
    if expand2dDim
        dataSize(3) = 1;
        blockSize(3) = 1;
    else
        blockSize = blockSize(1 : 2);
    end
end
blockSize = min(dataSize, blockSize);
if ~isempty(shardSize)
    shardSize = min(blockSize, shardSize);
    ind = rem(blockSize, shardSize) ~= 0;
    if any(ind)
        warning(['The shard size needs to be divisors of the block size, set', ...
            'the size of the axis that is not divisor to be the block size.']);
        shardSize(ind) = blockSize(ind);
    end
end
% overwrite the zarr file
if exist(filepath, 'dir')
    rmdir(filepath, 's');
end

switch dtype
    case 'double'
        ddtype = 'f8';
    case 'single'
        ddtype = 'f4';
    case 'uint16'
        ddtype = 'u2';
    case 'uint8'
        ddtype = 'u1';
    otherwise
        error('Unsupported data type %s', dtype);
end

try
    if ~isempty(zarrSubSize) && prod(ceil(dataSize ./ blockSize)) < prod(zarrSubSize)
        zarrSubSize = [];
    end
    if isempty(zarrSubSize)
        if isempty(shardSize)
            createZarrFile(filepath, 'chunks', blockSize, 'dtype', ddtype, 'order', order, ...
                'shape', dataSize, 'cname', compressor, 'clevel', clevel, 'dimension_separator', dimSeparator);
        else
            createZarrFile(filepath, 'chunks', blockSize, 'chunk_shape', shardSize, ...
                'dtype', ddtype, 'order', order, 'shape', dataSize, 'cname', compressor, ...
                'clevel', clevel, 'dimension_separator', dimSeparator);
        end
    else
        if isempty(shardSize) 
            createZarrFile(filepath, 'chunks', blockSize, 'dtype', ddtype, 'order', order, ...
                'shape', dataSize, 'cname', compressor, 'clevel', clevel, 'subfolders', zarrSubSize, ...
                'dimension_separator', dimSeparator);
        else
            createZarrFile(filepath, 'chunks', blockSize, 'chunk_shape', shardSize, ...
                'dtype', ddtype, 'order', order, 'shape', dataSize, 'cname', compressor, ...
                'clevel', clevel, 'subfolders', zarrSubSize, 'dimension_separator', dimSeparator);
        end
    end
catch ME
    pe = pyenv;    
    if pe.Status ~= "Loaded"
        throw(ME);
    end
    
    disp(ME);
    disp('Use the alternative zarr initializer (ZarrAdapter)...');    
    init_val = zeros(1, dtype);
    bim = blockedImage(filepath, dataSize, blockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');             
    bim.Adapter.close();
end

if isunix && groupWrite
    try
        fileattrib(filepath, '+w', 'g');
    catch
        warning('Unable to change file attribe for group write!');
    end
end

end
