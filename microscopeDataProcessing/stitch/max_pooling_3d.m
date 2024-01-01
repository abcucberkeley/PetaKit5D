function [out] = max_pooling_3d(orig, poolsize)
% max pooling 3d in matlab as a backup method

sz = size(orig);
orig = padarray(orig, ceil(sz ./ poolsize) .* poolsize - sz, 0, 'post');
orig = reshape(orig, poolsize(1), sz(1) ./ poolsize(1), poolsize(2), sz(2) ./ poolsize(2), poolsize(3), sz(3) ./ poolsize(3));
out = squeeze(max(orig, [], [1, 3, 5]));

end