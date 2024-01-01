function [mask] = decon_mask_edge_erosion(mask, EdgeErosion)
% erode the given mask based on direct assignment or 2d disk/line erosion


mask = mask> 0;

% for regular rectangler full image, directly set edges as false
if sum(mask, 'all') == numel(mask)
    mask([1 : EdgeErosion, end - EdgeErosion + 1 : end], :, :) = false;
    mask(:, [1 : EdgeErosion, end - EdgeErosion + 1 : end], :) = false;
    if ndims(mask) == 3
        mask(:, :, [1 : EdgeErosion, end - EdgeErosion + 1 : end]) = false;
    end
    return;
end

% for irregular shape image, use imerode to erode edges. 
mask([1, end], :, :) = false;
mask(:, [1, end], :) = false;
switch ndims(mask)
    case 2
        mask = imerode(mask, strel('disk', EdgeErosion - 1));
    case 3
        mask(:, :, [1, end]) = false;
        mask = imerode(mask, strel('disk', EdgeErosion - 1));
        mask = permute(mask, [3, 1, 2]);
        mask = imerode(mask, true(EdgeErosion * 2 - 1, 1));
        mask = permute(mask, [2, 3, 1]);
end

end