function [major_cover_mat, uncover_bbox_mat] = check_major_tile_cover(coords_mat, major_inds)
% check if all the major tile covers the region of minor tiles


major_cover_mat = true(numel(major_inds), 1);
uncover_bbox_mat = zeros(numel(major_inds), 6);
 
for i = 1 : numel(major_inds)
    if major_inds(i)
        continue;
    end
    coords_i = coords_mat(i, :);

    coords_i_mat = coords_i;
    cover_i_mat = false(1, 1);

    for j = 1 : numel(major_inds)
        if ~major_inds(j)
            continue;
        end
        if all(cover_i_mat)
            break;
        end
        
        coords_j = coords_mat(j, :);
        
        % check if current tile covers the (sub)cuboids from tile i
        is_cover_j = all(coords_j(1 : 3) <= coords_i_mat(:, 1 : 3) & coords_i_mat(:, 4 : 6) <= coords_j(4 : 6), 2);
        cover_i_mat = cover_i_mat | is_cover_j;
        
        % check if current tile overlaps the (sub)cuboids from tile i
        nol_1 = coords_i_mat(:, 4 : 6) < coords_j(1 : 3);
        nol_2 = coords_j(4 : 6) < coords_i_mat(:, 1 : 3);
        is_overlap_mat = ~any(nol_1 | nol_2, 2);
        
        cubs_diff_cell = cell(numel(is_overlap_mat), 1);
        for k = 1 : numel(is_overlap_mat)
            if is_cover_j(k)  || ~is_overlap_mat(k)
                continue;
            end
            
            coords_k = coords_i_mat(k, :);
            cubs_diff = cuboid_diff(coords_k, coords_j);
            cubs_diff_cell{k} = cubs_diff;
        end
        cubs_diff_mat = cat(1, cubs_diff_cell{:});
        if ~isempty(cubs_diff_mat)
            coords_i_mat = [coords_i_mat(~cover_i_mat & ~is_overlap_mat, :); cubs_diff_mat];
            cover_i_mat = [cover_i_mat(~cover_i_mat & ~is_overlap_mat); false(size(cubs_diff_mat, 1), 1)];
        end
    end
        
    major_cover_mat(i) = all(cover_i_mat);
    if ~major_cover_mat(i)
        uncover_bbox_mat(i, :) = [min(coords_i_mat(:, 1 : 3), [], 1), min(coords_i_mat(:, 4 : 6), [], 1)];
    end
end

end


function [cub_d] = cuboid_diff(cub_A, cub_B)
% A - B return at most four different cuboids, assume A and B have overlap

cub_d = [];
if all(cub_A(1 : 3) >= cub_B(1 : 3)) && all(cub_A(4 : 6) <= cub_B(4 : 6))
    return;
end

% get the overlap region
cub_overlap_p = sort([reshape(cub_A, [], 2), reshape(cub_B, [], 2)], 2);
cub_overlap_p = cub_overlap_p(:, [2, 3]);
cub_overlap_p = cub_overlap_p(:)';

% use the overlap region to get the difference
coords_d = [reshape(cub_overlap_p, [], 2), reshape(cub_A, [], 2)]';
% coords_d = sort(unique([reshape(cub_overlap_p, [], 2), reshape(cub_A, [], 2)]', 'rows'), 2);
ys = sort(unique(coords_d(:, 1)));
xs = sort(unique(coords_d(:, 2)));
zs = sort(unique(coords_d(:, 3)));
if numel(ys) == 1
    ys = [ys; ys];
end
if numel(xs) == 1
    xs = [xs; xs];
end
if numel(zs) == 1
    zs = [zs; zs];
end

[Y, X, Z] = ndgrid(1 : numel(ys)-1, 1 : numel(xs)-1, 1 : numel(zs)-1);
cub_d = [ys(Y(:)), xs(X(:)), zs(Z(:)), ys(Y(:)+1), xs(X(:)+1), zs(Z(:)+1)];

cub_d(all(cub_d == cub_overlap_p, 2), :) = [];

end

