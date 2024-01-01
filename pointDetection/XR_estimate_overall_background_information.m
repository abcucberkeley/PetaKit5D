function [backgroundInfo] = XR_estimate_overall_background_information(frame, varargin)
% The function is used to estimate some statistics for a 3D image: mean,
% median, std. 
% The idea is to crop some window from the boundary/corners of the image,
% and then take the median of the statistics as the overall background
% information. 
% 
% backgroundInfo contains estimated [mean, median, std];
% 
% Author: Xiongtao Ruan 11/08/2019
% 


ip = inputParser;
ip.addRequired('frame', @isnumeric);
ip.addParameter('WindowSize', [10, 2], @isnumeric);
ip.addParameter('WindowNum', 200, @isnumeric);
ip.addParameter('BoundaryDistanceRatio', 0.3, @isnumeric); % percentage of distance to the boundary
ip.addParameter('minValidRatio', 0.3, @isnumeric); % percentage of none-NaN points
ip.addParameter('randseed', 1, @isnumeric);
ip.parse(frame, varargin{:});

ws = ip.Results.WindowSize;
win_num = ip.Results.WindowNum;
randseed = ip.Results.randseed;

rng(randseed);

backgroundInfo_mat = zeros(win_num, 3);
% corp image to the by the bounding box
nan_mat = isnan(frame);
frame = replace_nan_with_value(frame, 0);
[frame, bbox] = crop_image_by_bounding_box(frame);
nan_mat = nan_mat(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
sz = size(frame);

bound_ratio = ip.Results.BoundaryDistanceRatio; % only within 30% sz to the boundary
bound_sz = round(sz * bound_ratio);

lb = min([bound_sz; ws(1), ws(1) ws(2)]);
ub = max([bound_sz; ws(1), ws(1), ws(2)]);

c_x = 1; c_y = 1; c_z = 1;
for i = 1 : win_num
    % make sure the centers of window are not nan
    start_flag = true;
    while start_flag || nan_mat(c_y, c_x,  c_z)
        c_x = randi([lb(2), ub(2)]);
        c_y = randi([lb(1), ub(1)]);
        c_z = randi([lb(3), ub(3)]);

        % decide which side
        rs = rand(3, 1) > 0.5;
        c_x = c_x .* (1 - rs(1)) + (sz(2) - rs(1)) * rs(1);
        c_y = c_y .* (1 - rs(2)) + (sz(1) - rs(1)) * rs(2);
        c_z = c_z .* (1 - rs(3)) + (sz(3) - rs(1)) * rs(3);
        start_flag = false;
    end
    xa = max(1, c_x - ws(1)) : min(sz(2), c_x + ws(1));
    ya = max(1, c_y - ws(1)) : min(sz(1), c_y + ws(1));
    za = max(1, c_z - ws(2)) : min(sz(3), c_z + ws(2));
    
    I_crop = frame(ya, xa, za);
    nan_mat_crop = nan_mat(ya, xa, za);
    
    if mean(isnan(I_crop), 'all') > 1 - ip.Results.minValidRatio
        backgroundInfo_mat(i, :) = nan;
        continue;
    end
    
    I_crop = I_crop(~nan_mat_crop);
    backgroundInfo_mat(i, :) = [mean(I_crop), median(I_crop), std(I_crop)];
end

backgroundInfo_mat = backgroundInfo_mat(~any(isnan(backgroundInfo_mat), 2), :);
backgroundInfo = median(backgroundInfo_mat);


end


