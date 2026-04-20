function [] = saveMIP_tiff_parser(frameFullname, MIPFullname, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullname', @ischar); 
ip.addRequired('MIPFullname', @ischar); 
ip.addParameter('dtype', 'uint16', @ischar);
ip.addParameter('axis', [0, 0, 1], @(x) isnumeric(x) || ischar(x));
ip.addParameter('inputBbox', [] , @(x) isempty(x) || isvector(x) || ischar(x));

ip.parse(frameFullname, MIPFullname, varargin{:});

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

saveMIP_tiff(frameFullname, MIPFullname, dtype=dtype, axis=axis, inputBbox=inputBbox);

end

