function saveMIP_zarr_parser(zarrFullname, MIPFullname, dtype, axis)
% save MIP from zarr use dask for efficient processing

% Author: Xiongtao Ruan (11/19/2020)
% xruan (07/11/2022): add support for user defined axis MIPs
% xruan (03/03/2023): simplify the function and use readzarr as default
% method other than dask

if ischar(axis)
    axis = str2num(axis);
end

saveMIP_zarr(zarrFullname, MIPFullname, dtype, axis);