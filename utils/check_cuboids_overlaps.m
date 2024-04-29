function [is_overlap_mat, cuboid_overlap_mat] = check_cuboids_overlaps(cuboid_1, cuboid_2, stitch2D)
% Determine whether two sets of polyhedra cuboids overlap with each other. 
% 
% cuboid_1 and cuboid_2 are n X 6 matrix: [min_x, min_y, min_z, max_x, max_y, max_z].
% 
% Author: Xiongtao Ruan
%
% xruan (05/06/2023): change input size and add support for multiple pairs
% of cuboids at once


if size(cuboid_1, 2) ~= 6 || size(cuboid_2, 2) ~= 6
    error('The size of the cuboids must be n X 6!');
end

nP_1 = size(cuboid_1, 1);
nP_2 = size(cuboid_2, 1);

if nP_1 < nP_2
    if nP_1 == 1
        nP = nP_2;
        cuboid_1 = repmat(cuboid_1, nP, 1);
    else
        error('The number of cuboids in the two sets must equal or one set need to have one cuboid!');
    end
elseif nP_1 > nP_2 
    if nP_2 == 1
        nP = nP_1;
        cuboid_2 = repmat(cuboid_2, nP, 1);
    else
        error('The number of cuboids in the two sets must equal or one set need to have one cuboid!');
    end
else
    nP = nP_1;
end

cuboid_11 = cuboid_1(:, 1 : 3);
cuboid_12 = cuboid_1(:, 4 : 6);
cuboid_21 = cuboid_2(:, 1 : 3);
cuboid_22 = cuboid_2(:, 4 : 6);

% ol_1 = cuboid_11 <= cuboid_21 & cuboid_21 <= cuboid_12;
% ol_2 = cuboid_11 <= cuboid_22 & cuboid_22 <= cuboid_12;
% ol_3 = cuboid_21 <= cuboid_11 & cuboid_11 <= cuboid_22;
% ol_4 = cuboid_21 <= cuboid_12 & cuboid_12 <= cuboid_22;
% 
% is_overlap_mat = all(ol_1 | ol_2 | ol_3 | ol_4, 2);

% check overlap from the non-overlap perspective
nol_1 = cuboid_22 < cuboid_11;
nol_2 = cuboid_12 < cuboid_21;
is_overlap_mat = ~any(nol_1 | nol_2, 2);

cuboid_overlap_mat = [];
if nargout > 1 && any(is_overlap_mat)
    cuboid_overlap_mat = zeros(nP, 6);        
    for p = 1 : nP
        if ~is_overlap_mat(p)
            continue;
        end
        cuboid_1p = cuboid_1(p, :);
        cuboid_2p = cuboid_2(p, :);
        
        cubiod_overlap_p = sort([reshape(cuboid_1p, [], 2), reshape(cuboid_2p, [], 2)], 2);
        cubiod_overlap_p = cubiod_overlap_p(:, [2, 3]);
        cuboid_overlap_mat(p, :) = cubiod_overlap_p(:);
        % if there is a singularity edge, return false in overlap.
        if (~stitch2D && any(cubiod_overlap_p(:, 1) == cubiod_overlap_p(:, 2))) || ...
                (stitch2D && any(cubiod_overlap_p(1 : 2, 1) == cubiod_overlap_p(1 : 2, 2)))
            cuboid_overlap_mat(p, :) = 0;
            is_overlap_mat(p) = false;
        end
    end 
end

end
