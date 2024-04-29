function [tim_block, mex_compute] = feather_blending_3d(tim_f_block, tim_d_block, tim_block, bbox)
% 3d feather blending function 
%
% if the out volume (tim_block) and the bounding box to place the overlap
% region is provided, use the mex function combined with indexing


mex_compute = true;
try 
    if nargin == 2
        tim_block = feather_blending_3d_mex(tim_f_block, tim_d_block);
    else
        feather_blending_3d_with_indexing_mex(tim_f_block, tim_d_block, tim_block, bbox);
    end
catch ME
    disp(ME);
    mex_compute = false;
    tim_w_block = tim_d_block .* (tim_f_block ~= 0); 
    if nargin == 2
        tim_block = sum(single(tim_f_block) .* tim_w_block, 4) ./ sum(tim_w_block, 4);
    else
        tim_block_crop = sum(single(tim_f_block) .* tim_w_block, 4) ./ sum(tim_w_block, 4);
        indexing4d(tim_block, tim_block_crop, [bbox(1 : 3), 1, bbox(4 : 6), 1]);
    end
end

end
