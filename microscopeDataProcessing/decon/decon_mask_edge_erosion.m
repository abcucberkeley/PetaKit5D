function [mask] = decon_mask_edge_erosion(mask, EdgeErosion)
% erode the given mask based on direct assignment or 2d disk/line erosion
%
% For irregular masks, decon_mask_edge_erosion_mex (CPU) or
% decon_mask_edge_erosion_gpu_mex (gpuArray masks), if available, apply the
% decomposition of strel('disk', EdgeErosion - 1) and the line along z as
% imerode does, with the same result, in one pass.


if ~islogical(mask)
    mask = mask > 0;
end

% for regular rectangler full image, directly set edges as false
if all(mask, 'all')
    mask([1 : EdgeErosion, end - EdgeErosion + 1 : end], :, :) = false;
    mask(:, [1 : EdgeErosion, end - EdgeErosion + 1 : end], :) = false;
    if ndims(mask) == 3
        mask(:, :, [1 : EdgeErosion, end - EdgeErosion + 1 : end]) = false;
    end
    return;
end

% in the mex, which clears the border and erodes with imerode's decomposition and padding
isGPU = isa(mask, 'gpuArray');
mexNames = {'decon_mask_edge_erosion_mex', 'decon_mask_edge_erosion_gpu_mex'};
if ndims(mask) <= 3 && ~issparse(mask) && isa(EdgeErosion, 'double') && isscalar(EdgeErosion) && EdgeErosion >= 1 ...
        && EdgeErosion == round(EdgeErosion) && exist(mexNames{isGPU + 1}, 'file') == 3
    [nhoods, padUL, padLR] = disk_decomposition(EdgeErosion - 1);
    rz = 0;
    if ndims(mask) == 3
        rz = EdgeErosion - 1;
    end
    if isGPU
        mask = decon_mask_edge_erosion_gpu_mex(mask, nhoods, padUL, padLR, rz, true);
    else
        mask = decon_mask_edge_erosion_mex(mask, nhoods, padUL, padLR, rz, true, maxNumCompThreads);
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


function [nhoods, padUL, padLR] = disk_decomposition(r)
% the neighborhoods of the decomposition of strel('disk', r) that imerode
% applies in turn, and its padding of the image for them (morphop: when the
% decomposition has several non-orthogonal parts, the padding of getpadsize,
% before and after swapped for erosion)

persistent lastR lastNhoods lastPadUL lastPadLR

if isempty(lastR) || lastR ~= r
    seq = decompose(strel('disk', r));
    lastNhoods = arrayfun(@getnhood, seq, 'unif', 0);
    lastPadUL = [0, 0];
    lastPadLR = [0, 0];
    if numel(seq) > 1 && ~isdecompositionorthogonal(seq)
        [padUL, padLR] = getpadsize(seq);
        lastPadUL = padLR(1 : 2);
        lastPadLR = padUL(1 : 2);
    end
    lastR = r;
end
nhoods = lastNhoods;
padUL = lastPadUL;
padLR = lastPadLR;

end
