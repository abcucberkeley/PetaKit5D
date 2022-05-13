function [is_overlap, cubiod_overlap] = cuboids_overlaps(cuboid_1, cuboid_2, stitch2D)
% Determine whether two polyhedra cuboids overlap with each other. 
% 
% cuboid_1 and cuboid_2 are 3 X 2 matrix: [min_x, max_x; min_y, max_y; min_z, max_z].
% 
% Author: Xiongtao Ruan


ol_1 = cuboid_1(:, 1) <= cuboid_2(:, 1) & cuboid_2(:, 1) <= cuboid_1(:, 2);
ol_2 = cuboid_1(:, 1) <= cuboid_2(:, 2) & cuboid_2(:, 2) <= cuboid_1(:, 2);
ol_3 = cuboid_2(:, 1) <= cuboid_1(:, 2) & cuboid_1(:, 2) <= cuboid_2(:, 2);
ol_4 = cuboid_2(:, 1) <= cuboid_1(:, 2) & cuboid_1(:, 2) <= cuboid_2(:, 2);

is_overlap = all(ol_1 | ol_2 | ol_3 | ol_4);

cubiod_overlap = [];
if is_overlap && nargout > 1
    cubiod_overlap = sort([cuboid_1, cuboid_2], 2);
    cubiod_overlap = cubiod_overlap(:, [2, 3]);
    % if there is a singularity edge, return false in overlap.
    if (~stitch2D && any(cubiod_overlap(:, 1) == cubiod_overlap(:, 2))) || ...
            (stitch2D && any(cubiod_overlap(1 : 2, 1) == cubiod_overlap(1 : 2, 2)))
        cubiod_overlap = [];
        is_overlap = false;
    end
end

end