function saveMIP_zarr(zarrFullname, MIPFullname, dtype, axis)
% save MIP from zarr use dask for efficient processing

% Author: Xiongtao Ruan (11/19/2020)
% xruan (07/11/2022): add support for user defined axis MIPs

if nargin < 3
    nv_bim = blockedImage(zarrFullname, 'Adapter', ZarrAdapter);
    dtype = nv_bim.ClassUnderlying;
end
if nargin < 4
    axis = [0, 0, 1];
end

uuid = get_uuid();

for i = 1 : 3
    if axis(i) == 0
        continue;
    end
    axis_i = i;
    switch i
        case 1    
            try
                % for matlab R2020b, it may fail to write tiff file if this step
                % runs before using Tiff write, to avoid the issue, we first write
                % a rand image to tmp folder to initialize Tiff in matlab
                writetiff(uint8(rand(10) > 0.5), sprintf('/tmp/%s.tif', uuid));
                MIP = py.daskAPI.daskZarrMaxProjection(zarrFullname, axis_i-1);    
            catch ME
                disp(ME);
                % this step is pretty slow in a single node, takes ~15min for 313 GB data
                nv_bim = blockedImage(zarrFullname, 'Adapter', CZarrAdapter);
                
                blockSize = min(5000, nv_bim.BlockSize * 10);
                blockSize(axis_i) = nv_bim.Size(axis_i);
                    
                bmip = apply(nv_bim, @(bs) max(bs.Data, [], axis_i), 'blockSize', blockSize);
                MIP = gather(bmip);        
            end
            MIPFullname = sprintf('%s_MIP_y.tif', MIPFullname(1 : end - 10));
        case 2 
            try
                writetiff(uint8(rand(10) > 0.5), sprintf('/tmp/%s.tif', uuid));
                MIP = py.daskAPI.daskZarrMaxProjection(zarrFullname, axis_i-1);    
            catch ME
                disp(ME);
                % this step is pretty slow in a single node, takes ~15min for 313 GB data
                nv_bim = blockedImage(zarrFullname, 'Adapter', CZarrAdapter);
                
                blockSize = min(5000, nv_bim.BlockSize * 10);
                blockSize(axis_i) = nv_bim.Size(axis_i);
                    
                bmip = apply(nv_bim, @(bs) max(bs.Data, [], axis_i), 'blockSize', blockSize);
                MIP = gather(bmip);        
            end
            MIPFullname = sprintf('%s_MIP_x.tif', MIPFullname(1 : end - 10));
        case 3 
            try
                writetiff(uint8(rand(10) > 0.5), sprintf('/tmp/%s.tif', uuid));
                MIP = py.daskAPI.daskZarrMaxProjection(zarrFullname, axis_i-1);    
            catch ME
                disp(ME);
                % this step is pretty slow in a single node, takes ~15min for 313 GB data
                try
                    nv_bim = blockedImage(zarrFullname, 'Adapter', CZarrAdapter);
                catch ME
                    disp(ME)
                    nv_bim = blockedImage(zarrFullname, 'Adapter', ZarrAdapter);
                end
                
                blockSize = min(nv_bim.Size, min(5000, nv_bim.BlockSize * 10));
                blockSize(axis_i) = nv_bim.Size(axis_i);
                    
                bmip = apply(nv_bim, @(bs) max(bs.Data, [], axis_i), 'blockSize', blockSize);
                MIP = gather(bmip);        
            end
            MIPFullname = sprintf('%s_MIP_z.tif', MIPFullname(1 : end - 10));
    end
    
    MIP = cast(MIP, dtype);
    MIP = squeeze(MIP);
    writetiff(MIP, MIPFullname);
    if ~exist(MIPFullname, 'file')
        imwrite(MIP, MIPFullname);
    end
end


end
