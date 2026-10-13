function [c] = centroids(vol, positions, radius)
% intensity weighted centroids (#spots x 3, [y, x, z]) of the spots near
% the given approximate positions, see pktest.centroid

c = zeros(size(positions, 1), 3);
for i = 1 : size(positions, 1)
    c(i, :) = pktest.centroid(vol, positions(i, :), radius);
end

end
