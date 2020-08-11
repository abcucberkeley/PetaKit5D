function [ output_args ] = rotate3DAll(folder, angle, aspratio, zrange)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

cd(folder);
mkdir('rotated')
tiffs=dir('*.tif');
for i=1:length(tiffs)
    img=loadtiff(tiffs(i).name);
    imgR=rotate3D(img, angle, aspratio);
    write3Dtiff(imgR(:,:,zrange(1):zrange(2)), strcat('rotated/', tiffs(i).name))
end

end

