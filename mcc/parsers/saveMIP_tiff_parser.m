function [] = saveMIP_tiff_parser(frameFullname, MIPFullname, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullname', @ischar); 
ip.addRequired('MIPFullname', @ischar); 
ip.addParameter('dtype', 'uint16', @ischar); % suffix for the folder
ip.addParameter('axis', [0, 0, 1], @(x) isnumeric(x) || ischar(x)); % suffix for the folder

ip.parse(frameFullname, MIPFullname, varargin{:});

pr = ip.Results;
dtype = pr.dtype;
axis = pr.axis;

if ischar(axis)
    axis = str2num(axis);
end

saveMIP_tiff(frameFullname, MIPFullname, dtype=dtype, axis=axis);

end

