function saveMIP_zarr_parser(zarrFullname, MIPFullname, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullname', @(x) ischar(x));
ip.addRequired('MIPFullname', @(x) ischar(x));
ip.addParameter('dtype',  '', @(x) ischar(x) || isstring(x));
ip.addParameter('axis',  [0, 0, 1], @(x) isvector(x) || numel(x) == 3 || ischar(x));
ip.addParameter('inputBbox', [] , @(x) isempty(x) || isvector(x) || ischar(x));

ip.parse(zarrFullname, MIPFullname, varargin{:});

pr = ip.Results;
dtype = pr.dtype;
axis = pr.axis;
inputBbox = pr.inputBbox;

if ischar(axis)
    axis = str2num(axis);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end

saveMIP_zarr(zarrFullname, MIPFullname, dtype=dtype, axis=axis, inputBbox=inputBbox);

end

