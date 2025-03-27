function [xyz, data_order_mat] = stitch_process_coordinates(xyz, imSizes, flippedTile, order_sign_mat, data_order_mat, xp, yp, zp, theta, IOScan, objectiveScan, DS, DSR)
% process coordinations for stitching for different spaces and data orders
% Here we define: y: light sheet propogation direction; x: light sheet
% direction; z: scanning direction.
% 
% author: Xiongtao Ruan (03/25/2025)


% adjust the x coordinates for flipped tiles.
if ~isempty(flippedTile)
    flippedTile = flippedTile > 0;
    if DSR    
        offsize = (imSizes(:, data_order_mat(1)) - 1) * xp - imSizes(:, data_order_mat(3)) * cot(theta) * zp;
        % offsize = - imSizes(:, 3) * cot(theta) * zf;
        % offsize = 0;
        xyz(:, 1) = xyz(:, 1) - offsize .* flippedTile(:) * order_sign_mat(2, 1);
    else
        offsize = (imSizes(:, data_order_mat(3)) - 1) * zp / sin(theta);
        xyz(:, 1) = xyz(:, 1) - offsize .* flippedTile(:) * order_sign_mat(2, 1);        
    end
end

if ~IOScan && ~objectiveScan && ~DS && ~DSR
    % convert coordinates in DSR space to skewned space
   xyz = [xyz(:, 3) / sin(theta), xyz(:, 2), -xyz(:, 1) * sin(theta) + xyz(:, 3) * cos(theta)];
end

if objectiveScan || DS
    % convert coordinates in Objective Scan / DS space
    xyz = [xyz(:, 1) * cos(theta) + xyz(:, 3) * sin(theta), xyz(:, 2), -xyz(:, 1) * sin(theta) + xyz(:, 3) * cos(theta)];
end

% identify pairs between pairs
% first slight shift xyz such that the distance between pairs are in
% integer folds of resolution.
xyz = round((xyz - min(xyz, [], 1)) ./ [xp, yp, zp]) .* [xp, yp, zp];

end
