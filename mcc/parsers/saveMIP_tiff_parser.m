function [] = saveMIP_tiff_parser(frameFullname, MIPFullname, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpath', @ischar); 
ip.addRequired('mipFullpath', @ischar); 
ip.addParameter('dtype', 'uint16', @ischar); % suffix for the folder
ip.addParameter('axis', [0, 0, 1], @(x) isnumeric(x) || ischar(x)); % suffix for the folder

ip.parse(frameFullname, MIPFullname, varargin{:});

pr = ip.Results;
dtype = pr.dtype;
axis = pr.axis;

if ischar(axis)
    axis = str2num(axis);
end

saveMIP_tiff(frameFullpath, mipFullpath, dtype=dtype, axis=axis);

end

