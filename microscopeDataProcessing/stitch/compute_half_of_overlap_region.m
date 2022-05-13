function [mregion_1, mregion_2] = compute_half_of_overlap_region(cuboid_1, cuboid_2, px, xyz_factors, varargin)
% reduce overlap region size for two images
% 
% 
% Author: Xiongtao Ruan (03/30/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('cuboid_1', @isnumeric);
ip.addRequired('cuboid_2', @isnumeric);
ip.addRequired('px', @isnumeric);
ip.addRequired('xyz_factors', @isnumeric);
ip.addParameter('overlapType', 'half', @ischar);
ip.addParameter('halfOrder', [2, 3, 1], @isnumeric);
ip.addParameter('stitch2D', false, @islogical);

ip.parse(cuboid_1, cuboid_2, px, xyz_factors, varargin{:});

overlapType = ip.Results.overlapType;
halfOrder = ip.Results.halfOrder;
stitch2D = ip.Results.stitch2D;

[is_overlap, cubiod_overlap] = cuboids_overlaps(cuboid_1, cuboid_2, stitch2D);

s1 = round((cubiod_overlap(:, 1) - cuboid_1(:, 1)) ./ (px * xyz_factors)) + 1;
s2 = round((cubiod_overlap(:, 1) - cuboid_2(:, 1)) ./ (px * xyz_factors)) + 1;

t1 = round((cubiod_overlap(:, 2) - cuboid_1(:, 1)) ./ (px * xyz_factors)) + 1;
t2 = round((cubiod_overlap(:, 2) - cuboid_2(:, 1)) ./ (px * xyz_factors)) + 1;


% split for the half of the 
olLen= t1 - s1 + 1;

split_axis = halfOrder(1);
if olLen(split_axis(1)) <= 1
    if olLen(halfOrder(2)) > 1
        split_axis = halfOrder(2);
    elseif olLen(halfOrder(3)) > 1
        split_axis = halfOrder(3);
    end
end

mregion_1 = [s1, t1];
mregion_2 = [s2, t2];

switch overlapType
    case 'zero'
        fc = 0;
    case 'none'
        fc = 0.5;
    case 'half'
        fc = 0.75;
    case 'full'
        fc = 1;
end

rm_len = round(olLen(split_axis)  * fc);
if t1(split_axis) > t2(split_axis)
    mregion_1(split_axis, 1) = mregion_1(split_axis, 1) + rm_len;
    mregion_2(split_axis, 2) = mregion_2(split_axis, 2) - rm_len;
else
    mregion_1(split_axis, 2) = mregion_1(split_axis, 2) - rm_len;
    mregion_2(split_axis, 1) = mregion_2(split_axis, 1) + rm_len;
end                

% convert to absolute positions
% mregion_1 = (mregion_1 - 1 ) .* (px * xyz_factors) + cuboid_1(:, 1);
% mregion_2 = (mregion_2 - 1 ) .* (px * xyz_factors) + cuboid_2(:, 1);


end
