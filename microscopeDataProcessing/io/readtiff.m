function [data] = readtiff(filepath, options)
% wrapper for tiff reader with both mex parallel version and matlab version code
% 
% range: [first, last] slice (finite values: the mex ended MATLAB for an Inf
% end and read one slice for a NaN end)
% 
% Author: Xiongtao Ruan (09/28/2022)


arguments
    filepath char
    options.range (1, :) {mustBeNumeric, mustBeFinite} = []
end

range = options.range;

try 
    if isempty(range)
        data = parallelReadTiff(filepath);
    else
        data = parallelReadTiff(filepath, range);
    end
catch ME
    mex_fallback('parallelReadTiff', ME, {'tiff:dataTypeError'});
    if isempty(range)
        data = readtiff_matlab(filepath);
    else
        data = readtiff_matlab(filepath, range(1) : range(2));
    end
end

end