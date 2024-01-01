function [bbox] = getImageBoundingBox(I)
% get image (2d or 3d) bbounding box

nd = ndims(I);

I_y = sum(I, [2, 3]);
I_x = sum(I, [1, 3]);

y1 = find(I_y > 0, 1, 'first');
y2 = find(I_y > 0, 1, 'last');

x1 = find(I_x > 0, 1, 'first');
x2 = find(I_x > 0, 1, 'last');

if nd == 3
    I_z = sum(I, [1, 2]);    
    z1 = find(I_z > 0, 1, 'first');
    z2 = find(I_z > 0, 1, 'last');
    bbox = [y1, x1, z1, y2, x2, z2];
else
    bbox = [y1, x1, y2, x2];
end

end