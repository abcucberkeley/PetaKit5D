function writezarr(data, filepath, varargin)
% wrapper for zarr writer 
% 
% Author: Xiongtao Ruan (01/25/2022)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isnumeric(x));
ip.addRequired('filepath', @(x) ischar(x));
ip.addParameter('blockSize', [500, 500, 500], @isnumeric);
ip.addParameter('groupWrite', true, @islogical);
ip.parse(data, filepath, varargin{:});

pr = ip.Results;
blockSize = pr.blockSize;
groupWrite = pr.groupWrite;

dtype = class(data);
sz = size(data);
init_val = zeros(1, dtype);
blockSize = min(sz, blockSize);
% overwrite the zarr file
if exist(filepath, 'dir')
    rmdir(filepath, 's');
end
bim = blockedImage(filepath, sz, blockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');
% for data greater than 2GB, use multiprocessing
if ~ispc && prod(sz) * 2 / 1024^3 > 2
    bim.Adapter.setData(data);
else
    bim.Adapter.setRegion([1, 1, 1], bim.Size, data)
end
bim.Adapter.close();

if groupWrite
    try
        fileattrib(filepath, '+w', 'g');
    catch
        warning('Unable to change file attribe for group write!');
    end
end

end