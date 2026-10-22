function saveMIP_zarr(zarrFullname, MIPFullname, varargin)
% save MIP from zarr use dask for efficient processing

% Author: Xiongtao Ruan (11/19/2020)
% xruan (07/11/2022): add support for user defined axis MIPs
% xruan (03/03/2023): simplify the function and use readzarr as default
% method other than dask


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullname', @(x) ischar(x));
ip.addRequired('MIPFullname', @(x) ischar(x));
ip.addParameter('dtype',  '', @(x) ischar(x) || isstring(x));
ip.addParameter('axis',  [0, 0, 1], @(x) isvector(x) || numel(x) == 3);
ip.addParameter('inputBbox', [] , @(x) isempty(x) || isvector(x));

ip.parse(zarrFullname, MIPFullname, varargin{:});

pr = ip.Results;
dtype = pr.dtype;
axis = pr.axis;
inputBbox = pr.inputBbox;

if isempty(dtype)
    dtype = getImageDataType(zarrFullname);
end

axis_strs = {'y', 'x', 'z'};
uuid = get_uuid();

try
    MIPs = saveMIP_zarr_slabs(zarrFullname, axis, inputBbox);
catch ME_1
    disp(ME_1)
    MIPs = cell(1, 3);
end
for i = 1 : 3
    if axis(i) == 0
        continue;
    end
    axis_i = i;
    fprintf('Generate MIP %s... ', axis_strs{i});
    MIP = MIPs{i};
    MIPs{i} = [];
    if isempty(MIP)
        MIP = saveMIP_zarr_axis(zarrFullname, axis_i, inputBbox);
    end
    MIP = squeeze(MIP);
    disp('Done!');
    MIPFullname = sprintf('%s_MIP_%s.tif', MIPFullname(1 : end - 10), axis_strs{i});

    MIP = cast(MIP, dtype);
    MIP = squeeze(MIP);
    % write to a temporary file and rename it: the job framework counts a MIP
    % task as done once its MIP exists
    MIPTmpname = sprintf('%s_%s.tif', MIPFullname(1 : end - 4), uuid);
    writetiff(MIP, MIPTmpname);
    if ~exist(MIPTmpname, 'file')
        imwrite(MIP, MIPTmpname);
    end
    movefile(MIPTmpname, MIPFullname);
end

end


function [MIPs] = saveMIP_zarr_slabs(zarrFullname, axis, inputBbox)
% MIPs of the image (or its inputBbox region) read in z slabs of whole chunk
% layers, so that the image is not held in memory at once: a running max
% along z (max keeps the first of equal values, as for the whole image), and
% the y and x MIPs of each slab placed along z. A slab has at least 64
% chunks, so that the reads stay parallel; one read for zarr files without
% .zarray.

zs = [];
try
    zInfo = getZarrInfo(zarrFullname);
    bbox = inputBbox;
    if isempty(bbox)
        bbox = [1, 1, 1, zInfo.size];
    end
    if numel(zInfo.size) == 3 && numel(bbox) == 6
        % starts of the slabs along z, at multiples of whole chunk layers
        cSize = zInfo.blockSize;
        layerChunkNum = prod(ceil(bbox(4 : 5) ./ cSize(1 : 2)) - floor((bbox(1 : 2) - 1) ./ cSize(1 : 2)));
        step = ceil(64 / layerChunkNum) * cSize(3);
        zs = [bbox(3), (floor((bbox(3) - 1) / step) + 1) * step + 1 : step : bbox(6)];
    end
catch ME
    disp(ME);
end

% range_max_3d_mex is the same as max, 2-3 times faster for uint16 (max if
% it isn't built for the platform)
useMex = true;
MIPs = cell(1, 3);
if numel(zs) <= 1
    im = readzarr(zarrFullname, 'inputBbox', inputBbox);
    for i = 1 : 3
        if axis(i) ~= 0
            [MIPs{i}, useMex] = axis_max(im, i, useMex);
        end
    end
    return;
end

ze = [zs(2 : end) - 1, bbox(6)];
for s = 1 : numel(zs)
    im = readzarr(zarrFullname, 'inputBbox', [bbox(1 : 2), zs(s), bbox(4 : 5), ze(s)]);
    for i = 1 : 3
        if axis(i) == 0
            continue;
        end
        [mip, useMex] = axis_max(im, i, useMex);
        if i == 3
            if s == 1
                MIPs{i} = mip;
            else
                MIPs{i} = max(MIPs{i}, mip);
            end
        else
            if s == 1
                mipSize = bbox(4 : 6) - bbox(1 : 3) + 1;
                mipSize(i) = 1;
                MIPs{i} = zeros(mipSize, class(mip));
            end
            MIPs{i}(:, :, zs(s) - bbox(3) + 1 : ze(s) - bbox(3) + 1) = mip;
        end
    end
    clear im;
end

end


function [mip, useMex] = axis_max(im, axis_ind, useMex)
% max along an axis, with range_max_3d_mex if it is available

if useMex
    try
        mip = range_max_3d_mex(im, axis_ind, 1);
        return;
    catch ME
        disp(ME);
        useMex = false;
    end
end
mip = max(im, [], axis_ind);

end


function [MIP] = saveMIP_zarr_axis(zarrFullname, axis_ind, inputBbox)
% MIP of one axis computed blockwise, for images that can't be read

if ~isempty(inputBbox)
    error('Batch MIP is not supported with input BBox yet!');
end

% this step is pretty slow in a single node, takes ~15min for 313 GB data
nv_bim = blockedImage(zarrFullname, 'Adapter', CZarrAdapter);

blockSize = min(5000, nv_bim.BlockSize * 10);
blockSize(axis_ind) = nv_bim.Size(axis_ind);

bmip = apply(nv_bim, @(bs) max(bs.Data, [], axis_ind), 'blockSize', blockSize);
MIP = gather(bmip);

end

