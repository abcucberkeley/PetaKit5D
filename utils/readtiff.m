function [data] = readtiff(filepath, varargin)
% wrapper for tiff reader with both mex parallel version and matlab version code
% 
%
% 
% Author: Xiongtao Ruan (09/28/2022)


ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('filepath');
ip.addParameter('range', [], @isnumeric); % z range (start and end z)
ip.parse(filepath, varargin{:});

pr = ip.Results;
range = pr.range;

try 
    if isstring(filepath)
        filepath = convertStringsToChars(filepath);
    end
    if isempty(range)
        data = parallelReadTiff(filepath);
    else
        data = parallelReadTiff(filepath, range);
    end
catch ME
    disp(ME);
    if isempty(range)
        data = readtiff_matlab(filepath);
    else
        data = readtiff_matlab(filepath, range(1) : range(2));
    end
end

end