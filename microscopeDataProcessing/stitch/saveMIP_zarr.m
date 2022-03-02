function saveMIP_zarr(zarrFullname, MIPFullname, dtype, axis)
% save MIP from zarr use dask for efficient processing

% Author: Xiongtao Ruan (11/19/2020)

if nargin < 3
    nv_bim = blockedImage(zarrFullname, 'Adapter', ZarrAdapter);
    dtype = nv_bim.ClassUnderlying;
end
if nargin < 4
    axis = 3;
end

uuid = get_uuid();
try
    % for matlab R2020b, it may fail to write tiff file if this step
    % runs before using Tiff write, to avoid the issue, we first write
    % a rand image to tmp folder to initialize Tiff in matlab
    writetiff(uint8(rand(10) > 0.5), sprintf('/tmp/%s.tif', uuid));
    MIP = py.daskAPI.daskZarrMaxProjection(zarrFullname, axis-1);    
catch ME
    disp(ME);
    % this step is pretty slow in a single node, takes ~15min for 313 GB data
    nv_bim = blockedImage(zarrFullname, 'Adapter', ZarrAdapter);
    
    blockSize = min(5000, nv_bim.BlockSize * 10);
    blockSize(axis) = nv_bim.Size(axis);
        
    bmip = apply(nv_bim, @(bs) max(bs.Data, [], axis), 'blockSize', blockSize);
    MIP = gather(bmip);        
end

MIP = cast(MIP, dtype);
MIP = squeeze(MIP);
writetiff(MIP, MIPFullname);
if ~exist(MIPFullname, 'file')
    imwrite(MIP, MIPFullname);
end

end
