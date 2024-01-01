function [frame_out] = project3DImageto2D(frame, method)
% project 3d image to 2d
% supported methods: central xy, central xz, central yz, MIP_xy, MIP_yz,
%   MIP_xz, mean_xy, mean_yz, mean_xz
% 
% Author: Xiongtao Ruan, January 2020


switch lower(method)
    case 'central_xy'
        frame_out = frame(:, :, floor((size(frame, 3) + 1) / 2));
    case 'central_yz'
        frame_out = squeeze(frame(:, floor((size(frame, 2) + 1) / 2), :));
    case 'central_xz'
        frame_out = squeeze(frame(floor((size(frame, 1) + 1) / 2), :, :));
    case 'mip_xy'
        frame_out = max(frame, [], 3);
    case 'mip_yz'
        frame_out = squeeze(max(frame, [], 2));
    case 'mip_xz'
        frame_out = squeeze(max(frame, [], 1));
    case 'mean_xy'
        frame_out = nanmean(frame, 3);
    case 'mean_yz'
        frame_out = squeeze(nanmean(frame, 2));
    case 'mean_xz'
        frame_out = squeeze(nanmean(frame, 1));                
    otherwise
        frame_out = frame;
end


end