function [dmat] = feather_distance_map_resize_3d(dmat, bbox, wd)
% feather blending distance map resize 


try
    dmat = feather_distance_map_resize_3d_mex(dmat, bbox(4 : 6) - bbox(1 : 3) + 1, wd);
catch ME
    disp(ME);
    % im_d_j =  im_d_j .^ (1 / wd);
    if bbox(3) == bbox(6)
        dmat = imresize(dmat, bbox(4 : 5) - bbox(1 : 2) + 1, 'bilinear');
    else
        if size(dmat, 3) == 1 
            dmat = repmat(dmat, 1, 1, 2);
        end
        dmat = imresize3(dmat, bbox(4 : 6) - bbox(1 : 3) + 1, 'linear');
    end
    dmat = fastPower(dmat, wd);
end

end
