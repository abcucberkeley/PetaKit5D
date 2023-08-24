function [dist_weight] = distance_weight_single_axis(sz, endPoints, bufferSize, dfactor, winType)
% calculate distance weight based on the endpoints and buffer size for
% certain window type
%
% output as single

if nargin < 3
    bufferSize = 10;
end
% decay factor
if nargin < 4
    dfactor = 0.99;
end
if nargin < 5
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
win_weight = max(win_weight, 1e-3);

s0 = max(1, s - bufferSize);
dist_weight(s0 : s) = win_weight(bufferSize - (s - max(1, s - bufferSize)) + 1 : end);

% lower bound in case of the end weight is zero. 
lb = 1e-4 .* (dfactor > 0);
ub = 1e-3 .* (dfactor > 0);
if s0 > 1
    dist_weight(1 : s0 - 1) = dfactor .^ (s0 - 1 : -1 : 1) .* min(max(lb, dist_weight(s0)), ub);
end

win_weight_f = flip(win_weight);
t1 = min(sz, t + bufferSize);
dist_weight(t : t1) = win_weight_f(1 : min(sz, t + bufferSize) - t + 1);
if t1 < sz
    dist_weight(t1 + 1 : sz) = dfactor .^ (1 : sz - t1) .* min(max(lb, dist_weight(t1)), ub);
end

if dfactor > 0
    dist_weight = max(dist_weight, 1e-30);
end

dist_weight = single(dist_weight);

end

