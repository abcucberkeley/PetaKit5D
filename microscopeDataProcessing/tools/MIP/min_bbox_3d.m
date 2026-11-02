function [out] = min_bbox_3d(orig, bbox)
% min 3d for a given bbox matlab wrapper

try
    out = min_bbox_3d_mex(orig, bbox);
catch ME
    mex_fallback('min_bbox_3d_mex', ME, {'mat:dataTypeError'});

    % the minimum of the region, as the mex: starting from the largest value
    % of the class, so NaN (never smaller) is skipped and an empty region
    % gives that value
    region = crop3d(orig, bbox);
    if islogical(region)
        maxValue = true;
    elseif isinteger(region)
        maxValue = intmax(class(region));
    else
        maxValue = realmax(class(region));
    end
    out = min([maxValue; region(:)]);
end

end
