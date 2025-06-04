function [out] = max_pooling_3d(orig, poolsize)
% max pooling 3d matlab wrapper

try
    out = max_pooling_3d_mex(orig, poolsize);
catch ME
    disp(ME);

    sz = size(orig);
    orig = padarray(orig, ceil(sz ./ poolsize) .* poolsize - sz, 0, 'post');
    sz = size(orig);    
    orig = reshape(orig, poolsize(1), sz(1) ./ poolsize(1), poolsize(2), sz(2) ./ poolsize(2), poolsize(3), sz(3) ./ poolsize(3));
    out = squeeze(max(orig, [], [1, 3, 5]));
end

end

