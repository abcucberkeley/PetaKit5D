function saveMIP_zarr_parser(zarrFullname, MIPFullname, dtype, axis)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullname', @(x) ischar(x));
ip.addRequired('MIPFullname', @(x) ischar(x));
ip.addOptional('dtype',  '', @(x) ischar(x) || isstring(x));
ip.addOptional('axis',  [0, 0, 1], @(x) isvector(x) || numel(x) == 3 || ischar(x));

ip.parse(zarrFullname, MIPFullname, dtype, axis);

pr = ip.Results;

if ischar(axis)
    axis = str2num(axis);
end

saveMIP_zarr(zarrFullname, MIPFullname, dtype, axis);

end

