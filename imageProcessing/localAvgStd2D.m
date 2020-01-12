% Computes the local average and standard deviation within a square window of side 'w'.

% Francois Aguet
% Last modified on 09/19/2011

function [avg sigma] = localAvgStd2D(img, w)

if mod(w+1, 2)
    error('The window length w must be an odd integer.');
end

nanMask = isnan(img);

% kernel
h = ones(w,w);

% count of non-NaN elements
n = imfilter(double(~nanMask), h, 'replicate');
img(nanMask) = 0;

E = imfilter(img, h, 'replicate');
E2 = imfilter(img.^2, h, 'replicate');

sigma = E2 - E.^2./n;
sigma(sigma<0) = 0;
sigma = sqrt(sigma./(n - 1));
avg = E./n;

avg(nanMask) = NaN;
sigma(nanMask) = NaN;
