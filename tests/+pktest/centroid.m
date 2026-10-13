function [c] = centroid(vol, center, radius)
% intensity weighted centroid [y, x, z] of the spot near center, using the
% voxels within a cube of half width radius around round(center). The
% window minimum is subtracted so a flat background doesn't bias it.

center = round(center);
lo = max(center - radius, 1);
hi = min(center + radius, size(vol, 1 : 3));
win = double(vol(lo(1) : hi(1), lo(2) : hi(2), lo(3) : hi(3)));
win = win - min(win(:));

[y, x, z] = ndgrid(lo(1) : hi(1), lo(2) : hi(2), lo(3) : hi(3));
c = [sum(y(:) .* win(:)), sum(x(:) .* win(:)), sum(z(:) .* win(:))] / sum(win(:));

end
