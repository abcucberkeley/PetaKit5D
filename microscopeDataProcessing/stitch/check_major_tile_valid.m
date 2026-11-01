function [is_valid, unvalid_bbox] = check_major_tile_valid(fmat, major_inds)
% check if all major tiles are valid at all xyz voxels, that is, check if
% all elements are zeros in all major tiles

% is_valid_mat = false(numel(major_inds), 1);
% the mex takes uint16 and single tiles and at least one major tile
if (isa(fmat, 'uint16') || isa(fmat, 'single')) && any(major_inds)
    try
        % is_valid_mat(major_inds) = any_4th_dim_mex(fmat, major_inds);
        [is_valid, unvalid_bbox] = any_4th_dim_mex(fmat, major_inds);
        if ~is_valid
            unvalid_bbox = double(unvalid_bbox);
        else
            unvalid_bbox = [];
        end
        return;
    catch ME
        mex_fallback('any_4th_dim_mex', ME);
    end
end

% is_valid_mat(major_inds) = squeeze(all(fmat(:, :, :, major_inds), 1:3));
is_valid = squeeze(all(any(fmat(:, :, :, major_inds), 4), 1:3));
unvalid_bbox = [];
if ~is_valid
    im_i = ~any(fmat(:, :, :, major_inds), 4);
    inds = find(im_i);
    [y, x, z] = ind2sub(size(im_i), inds);
    unvalid_bbox = [min(y), min(x), min(z), max(y), max(x), max(z)];
end

end
