function [imn] = XR_normalize_z_stack(img_in)
% normalize z slice in 3D volume data
% 
% Author: Xiongtao Ruan (04/09/2020)



bg = 105;
T = bg;
im = img_in;

imT = im - T;

Med = zeros(size(im, 3), 1);
for j = 1:size(im,3)
    ims = imT(:,:,j);
    % Avg(j) = mean(ims(ims>0));
    Med(j) = median(ims(ims>0));
end
% figure,
% plot(Avg)
% hold on
% plot(Med)

m = nanmedian(Med);
spd = 10;
n = Med;
n(isnan(n)) = m;
n(n<m-spd) = m-spd;
n(n>m+spd) = m+spd;
% figure,
% plot(Avg)
% hold on
% plot(n)
nn = n/max(n);
% figure, plot(nn)
imn = zeros(size(im), class(im));
for j = 1:size(im,3)
  imn(:,:,j) = imT(:,:,j)/nn(j);  
end
% Tn = thresholdOtsu(imn(:));
% imn = imn-Tn;


end

