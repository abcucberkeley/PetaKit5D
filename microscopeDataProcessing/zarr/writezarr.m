function writezarr(data, filepath, varargin)
% wrapper for zarr writer 
% 
% Author: Xiongtao Ruan (01/25/2022)
% 
% xruan (05/23/2022): add support for bounding box write; also add
% parallelWriteZarr as default method


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isnumeric(x));
ip.addRequired('filepath', @(x) ischar(x) || isstring(x));
ip.addParameter('overwrite', false, @islogical);
ip.addParameter('blockSize', [500, 500, 500], @isnumeric);
ip.addParameter('expand2dDim', true, @islogical); % expand the z dimension for 2d data
ip.addParameter('groupWrite', true, @islogical);
ip.addParameter('bbox', [], @isnumeric);
ip.parse(data, filepath, varargin{:});

pr = ip.Results;
overwrite = pr.overwrite;
expand2dDim = pr.expand2dDim;
blockSize = pr.blockSize;
groupWrite = pr.groupWrite;
bbox = pr.bbox;

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
if ~overwrite && exist(filepath, 'dir')
    if isempty(bbox)
        bbox = [1, 1, 1, size(data, 1 : 3)];
        if ~expand2dDim
            bbox = [1, 1, size(data, 1 : 2)];
        end
    end
end

if isstring(filepath)
    filepath = char(filepath);
end
try 
    if ismatrix(data)
        % error('No support for 2d data for now!')
    end
    if isempty(bbox)
        parallelWriteZarr(filepath, data, 1, blockSize);
    else
        parallelWriteZarr(filepath, data, 1, bbox);
    end
catch ME
    disp(ME);
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
