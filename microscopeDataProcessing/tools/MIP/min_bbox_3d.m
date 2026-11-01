function [out] = min_bbox_3d(orig, bbox)
% min 3d for a given bbox matlab wrapper

try
    out = min_bbox_3d_mex(orig, bbox);
catch ME
    mex_fallback('min_bbox_3d_mex', ME, {'mat:dataTypeError'});
    
    out = min(crop3d(orig, bbox));
end

end
