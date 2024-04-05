function [] = stitch_global_grid_assignment_wrapper_parser(inputFullpath, outputFullpath, xcorr_thresh, axisWeight, uuid)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpath', @(x) ischar(x));
ip.addRequired('outputFullpath', @(x) ischar(x));
ip.addRequired('xcorr_thresh', @(x) isscalar(x) || ischar(x));
ip.addRequired('axisWeight', @(x) isvector(x) || ischar(x));
ip.addRequired('uuid', @(x) ischar(x));

ip.parse(inputFullpath, outputFullpath, xcorr_thresh, axisWeight, uuid);

pr = ip.Results;

if ischar(xcorr_thresh)
    xcorr_thresh = str2num(xcorr_thresh);
end
if ischar(axisWeight)
    axisWeight = str2num(axisWeight);
end

stitch_global_grid_assignment_wrapper(inputFullpath, outputFullpath, xcorr_thresh, ...
    axisWeight, uuid);

end

