function [out] = min_bbox_3d(orig, bbox)
% min 3d for a given bbox matlab wrapper

try
    out = min_bbox_3d_mex(orig, bbox);
catch ME
    disp(ME);
    
    out = min(crop3d(orig, bbox));
end

end
