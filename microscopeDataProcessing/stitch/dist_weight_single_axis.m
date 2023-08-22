function [dist_weight] = dist_weight_single_axis(sz, endPoints, bufferSize, winType)
% calculate distance weight based on the endpoints and buffer size for
% certain window type

if nargin < 3
    bufferSize = 10;
end
if nargin < 4
    winType = 'hann';
end

dist_weight = ones(sz, 1);

s = endPoints(1);
t = endPoints(2);

if s == 1 && t == sz
    return;
end

dist_weight(1 : s) = 0;
dist_weight(t : end) = 0;

win_dist = 0 : bufferSize;

switch lower(winType)
    case 'hann'
        win_func = @(x, y) cos(0.5 * pi * (y - x) ./ y) .^ 2;
    otherwise
        error('Unsupported window type %s', winType);
end

win_weight = win_func(win_dist, bufferSize);

dist_weight(max(1, s - bufferSize) : s) = win_weight(bufferSize - (s - max(1, s - bufferSize)) + 1 : end);
win_weight_f = flip(win_weight);
dist_weight(t : min(sz, t + bufferSize)) = win_weight_f(1 : min(sz, t + bufferSize) - t + 1);

end

