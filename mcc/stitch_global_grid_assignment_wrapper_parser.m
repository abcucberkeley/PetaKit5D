function [] = stitch_global_grid_assignment_wrapper_parser(inputFullpath, outputFullpath, xcorr_thresh, axisWeight, uuid)
% wrapper for grid assigment optimization used for parallel optimization
% for multiple groups


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpath', @(x) ischar(x));
ip.addRequired('outputFullpath', @(x) ischar(x));
ip.addRequired('xcorr_thresh', @(x) isnumeric(x) || ischar(x));
ip.addRequired('axisWeight', @(x) isnumeric(x) || ischar(x));
ip.addRequired('uuid', @(x) ischar(x));

ip.parse(inputFullpath, outputFullpath, xcorr_thresh, axisWeight, uuid);

if ischar(xcorr_thresh)
    xcorr_thresh = str2num(xcorr_thresh);
end
if ischar(axisWeight)
    axisWeight = str2num(axisWeight);
end

stitch_global_grid_assignment_wrapper(inputFullpath, outputFullpath, xcorr_thresh, axisWeight, uuid);

end
