function rotate3D_wrapper(tiff, angle, aspratio, z_trans, zrange)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ischar(angle)
    angle = str2double(angle);
end

if ischar(aspratio)
    aspratio = str2double(aspratio);
end

if ischar(z_trans)
    z_trans = str2double(z_trans);
end

if ischar(zrange)
    zrange = str2num(zrange);
end

img = loadtiff(tiff);

imgR=rotate3D(img, angle, aspratio, z_trans);
write3Dtiff(imgR(:,:,zrange(1):zrange(2)), strcat('rotated/', tiff))

end

