function [] = saveMIP_tiff_parser(frameFullname, MIPFullname, varargin)
% get the mip for a given image with image path
% 
% Author: Xiongtao Ruan (03/30/2021)
% 
% xruan (05/05/2022): add support for user defined axis MIPs

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpath'); 
ip.addRequired('mipFullpath'); 
ip.addParameter('dtype', 'uint16', @ischar); % suffix for the folder
ip.addParameter('axis', [0, 0, 1], @(x) isnumeric(x) || ischar(x)); % suffix for the folder

ip.parse(frameFullname, MIPFullname, varargin{:});

if exist(MIPFullname, 'file')
    fprintf('The mip file %s already exists, skip it!\n', MIPFullname)
    return;
end

dtype = ip.Results.dtype;
axis = ip.Results.axis;

if ischar(axis)
    axis = str2num(axis);
end

saveMIP_tiff(frameFullname,MIPFullname,'dtype',dtype,'axis',axis);