function writetiff(img, filepath, varargin)
% wrapper for saveastiff.m to support big tiff format, and keep
% compatibility of previous script for the decrapted writetiff. 
% 
% Author: Xiongtao Ruan (07/30/2020)
% xruan (08/07/2020): add option for group write

ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('Compression', 'lzw', @(x) any(strcmpi(x, {'none', 'lzw'})));
ip.addParameter('Mode', 'libtiff', @(x) any(strcmpi(x, {'libtiff', 'imwrite'})));
ip.addParameter('groupWrite', true, @islogical);
ip.parse(varargin{:});

pr = ip.Results;

options = struct();
switch pr.Compression
    case 'none'
        options.compress = 'no';
    case 'lzw'
        options.compress = 'lzw';
end

options.message = false;
options.overwrite = true;
saveastiff(img, filepath, options);
if pr.groupWrite
    try
        fileattrib(filepath, '+w', 'g');
    catch
        warning('Unable to change file attribe for group write!');
    end
end

end

