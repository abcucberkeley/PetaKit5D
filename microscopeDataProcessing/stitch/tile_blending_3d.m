function [nv_block] = tile_blending_3d(tim_f_block, mboxes, BlendMethod)
% blend the tiles of a stitching batch for the blend methods 'none', 'mean',
% 'max' and 'median'
%
% tim_f_block: the tiles of the batch in the 4th dimension (0 means no data)
% mboxes: [tile, ymin, xmin, zmin, ymax, xmax, zmax] per row, regions where a
% tile counts as 0 (the overlap regions that other tiles fill)
%
% A voxel gets the blend of the tiles outside these regions, or, where that is
% 0, the blend of all tiles. The mex computes it in one pass (uint8, uint16
% and single tiles), with the result in the class of the tiles for 'none' and
% 'max'; otherwise MATLAB computes it in single.


if ~any(strcmp(BlendMethod, {'none', 'mean', 'max', 'median'}))
    error('Unsupported blending method: %s.', BlendMethod);
end

if (isa(tim_f_block, 'uint16') || isa(tim_f_block, 'uint8') || isa(tim_f_block, 'single')) ...
        && exist('tile_blending_3d_mex', 'file') == 3
    try
        nv_block = tile_blending_3d_mex(tim_f_block, mboxes, BlendMethod);
        return;
    catch ME
        disp(ME);
    end
end

tim_block = tim_f_block;
for k = 1 : size(mboxes, 1)
    tim_block(mboxes(k, 2) : mboxes(k, 5), mboxes(k, 3) : mboxes(k, 6), mboxes(k, 4) : mboxes(k, 7), mboxes(k, 1)) = 0;
end

% convert to single for processing
tim_block = single(tim_block);
tim_f_block = single(tim_f_block);

if ~strcmp(BlendMethod, 'none')
    tim_block(tim_block == 0) = nan;
    tim_f_block(tim_f_block == 0) = nan;
end

% merge voxel info for the tiles
switch BlendMethod
    case 'none'
        for j = 1 : size(tim_block, 4)
            if j == 1
                nv_block = tim_block(:, :, :, j);
                nv_f_block = tim_f_block(:, :, :, j);
            else
                nv_block = nv_block + tim_block(:, :, :, j) .* (nv_block == 0);
                nv_f_block = nv_f_block + tim_f_block(:, :, :, j) .* (nv_f_block == 0);
            end
        end
    case 'mean'
        nv_block = nanmean(tim_block, 4);
        nv_f_block = nanmean(tim_f_block, 4);
    case 'max'
        nv_block = max(tim_block, [], 4);
        nv_f_block = max(tim_f_block, [], 4);
    case 'median'
        nv_block = nanmedian(tim_block, 4);
        nv_f_block = nanmedian(tim_f_block, 4);
end

if ~strcmp(BlendMethod, 'none')
    nv_block = replace_nans(nv_block);
end
nv_f_block = replace_nans(nv_f_block);
% fill the voxels the block doesn't cover yet
nv_zero_inds = nv_block == 0 & nv_f_block ~= 0;
nv_block = nv_block .* (1 - nv_zero_inds) + nv_f_block .* nv_zero_inds;

end


function [x] = replace_nans(x)
% NaNs to 0, with the mex where it is built (not for ARM macOS)

if exist('replace_nan_with_value', 'file') == 3
    try
        x = replace_nan_with_value(x, 0);
        return;
    catch ME
        disp(ME);
    end
end
x(isnan(x)) = 0;

end
