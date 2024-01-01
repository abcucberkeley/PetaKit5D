function [data] = readtiff(filepath, options)
% wrapper for tiff reader with both mex parallel version and matlab version code
% 
%
% 
% Author: Xiongtao Ruan (09/28/2022)


arguments
    filepath char
    options.range (1, :) {mustBeNumeric} = []
end

range = options.range;

try 
    if isempty(range)
        data = parallelReadTiff(filepath);
    else
        data = parallelReadTiff(filepath, range);
    end
catch ME
    disp(ME);
    disp('Use the alternative tiff reader (matlab libtiff)...');
    if isempty(range)
        data = readtiff_matlab(filepath);
    else
        data = readtiff_matlab(filepath, range(1) : range(2));
    end
end

end