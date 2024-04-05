function saveMIP_zarr(zarrFullname, MIPFullname, dtype, axis)
% save MIP from zarr use dask for efficient processing

% Author: Xiongtao Ruan (11/19/2020)
% xruan (07/11/2022): add support for user defined axis MIPs
% xruan (03/03/2023): simplify the function and use readzarr as default
% method other than dask


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullname', @(x) ischar(x));
ip.addRequired('MIPFullname', @(x) ischar(x));
ip.addOptional('dtype',  '', @(x) ischar(x) || isstring(x));
ip.addOptional('axis',  [0, 0, 1], @(x) isvector(x) || numel(x) == 3);

ip.parse(zarrFullname, MIPFullname, dtype, axis);

pr = ip.Results;
dtype = pr.dtype;
axis = pr.axis;

if isempty(dtype)
    dtype = getImageDataType(zarrFullname);
end

axis_strs = {'y', 'x', 'z'};

for i = 1 : 3
    if axis(i) == 0
        continue;
    end
    axis_i = i;
    MIP = saveMIP_zarr_axis(zarrFullname, axis_i);
    MIPFullname = sprintf('%s_MIP_%s.tif', MIPFullname(1 : end - 10), axis_strs{i});
    
    MIP = cast(MIP, dtype);
    MIP = squeeze(MIP);
    writetiff(MIP, MIPFullname);
    if ~exist(MIPFullname, 'file')
        imwrite(MIP, MIPFullname);
    end
end

end


function [MIP] = saveMIP_zarr_axis(zarrFullname, axis_ind)

% try
%     uuid = get_uuid();    
%     % for matlab R2020b, it may fail to write tiff file if this step
%     % runs before using Tiff write, to avoid the issue, we first write
%     % a rand image to tmp folder to initialize Tiff in matlab
%     writetiff(uint8(rand(10) > 0.5), sprintf('/tmp/%s.tif', uuid));
%     MIP = py.daskAPI.daskZarrMaxProjection(zarrFullname, axis_ind - 1);    
% catch ME
%     disp(ME);
%     disp('Use matlab image block method for MIP computing...');
%     try
%         im = readzarr(zarrFullname);
%         MIP = squeeze(max(im, [], axis_ind));
%     catch ME_1
%         disp(ME_1)
% 
%         % this step is pretty slow in a single node, takes ~15min for 313 GB data
%         nv_bim = blockedImage(zarrFullname, 'Adapter', CZarrAdapter);
%         
%         blockSize = min(5000, nv_bim.BlockSize * 10);
%         blockSize(axis_ind) = nv_bim.Size(axis_ind);
%             
%         bmip = apply(nv_bim, @(bs) max(bs.Data, [], axis_ind), 'blockSize', blockSize);
%         MIP = gather(bmip);
%     end
%     disp('Done!');    
% end

try
    im = readzarr(zarrFullname);
    MIP = squeeze(max(im, [], axis_ind));
catch ME_1
    disp(ME_1)

    % this step is pretty slow in a single node, takes ~15min for 313 GB data
    nv_bim = blockedImage(zarrFullname, 'Adapter', CZarrAdapter);
    
    blockSize = min(5000, nv_bim.BlockSize * 10);
    blockSize(axis_ind) = nv_bim.Size(axis_ind);
        
    bmip = apply(nv_bim, @(bs) max(bs.Data, [], axis_ind), 'blockSize', blockSize);
    MIP = gather(bmip);
end
disp('Done!');    


end

