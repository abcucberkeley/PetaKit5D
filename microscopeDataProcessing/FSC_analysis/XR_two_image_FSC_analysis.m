function [fsc, res] = XR_two_image_FSC_analysis(im1, im2, varargin)
% perform FSC analysis for two images


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('im1', @(x) ischar(x) || isnumeric(x));
ip.addRequired('im2', @(x) ischar(x) || isnumeric(x));
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.108, @isnumeric);
ip.addParameter('dr', 1 , @isnumeric);
ip.addParameter('dtheta', pi / 12 , @isnumeric);
ip.addParameter('resThreshMethod', 'one-bit', @ischar);
ip.addParameter('resThresh', 0.2, @isnumeric);
ip.addParameter('N', [501, 501, 501], @isnumeric);
ip.addParameter('debug', false, @islogical);

ip.parse(im1, im2, varargin{:});

pr = ip.Results;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
dr = pr.dr;
dtheta = pr.dtheta;
resThreshMethod = pr.resThreshMethod;
resThresh = pr.resThresh;
N = pr.N;
debug = pr.debug;

if ischar(im1)
    im1 = readtiff(im1);
end
im1 = double(im1);

if ischar(im2)
    im2 = readtiff(im2);
end
im2 = double(im2);

% resample, crop, pad images to isotropic resolution
im1 = resample_crop_pad_image_isotropic(im1, xyPixelSize, dz, N * 2);
im2 = resample_crop_pad_image_isotropic(im2, xyPixelSize, dz, N * 2);

im1 = reshape(im1, 2, N(1), 2, N(2), 2, N(3));
im1 = permute(im1, [2, 4, 6, 1, 3, 5]);

im2 = reshape(im2, 2, N(1), 2, N(2), 2, N(3));
im2 = permute(im2, [2, 4, 6, 1, 3, 5]);

% fsc computing
fsc_cell = cell(8, 1);
res_cell = cell(8, 1);
for i = 1 : 8
    [y, x, z] = ind2sub([2, 2, 2], i);
    [fsc_i, res_i] = XR_FSC(im1(:, :, :, y, x, z), im2(:, :, :, y, x, z), ...
        'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, 'N', N);
    fsc_cell{i} = fsc_i;
    res_cell{i} = res_i;    
end

% resolution computing
fsc_mu = fsc_i;
fsc_vals = cellfun(@(x) x.fsc, fsc_cell, 'unif', 0);
fsc_mu.fsc = mean(cat(3, fsc_vals{:}), 3);

[res_mu] = XR_FSC_resolution(fsc_mu, 'xyPixelSize', xyPixelSize, 'dz', dz, 'resThreshMethod', resThreshMethod, ...
    'resThresh', resThresh, 'debug', debug);

fsc = fsc_mu;
res = res_mu;

end


function [im] = resample_crop_pad_image_isotropic(im, xyPixelSize, dz, N)
% resample, crop, pad images to isotropic resolution for given cubic size


psz = min(xyPixelSize, dz);
if dz ~= xyPixelSize
    im = imresize3(im, round(size(im) .* [xyPixelSize, xyPixelSize, dz] / psz));
end

% crop image if any size is greater than the 2 folder of split images
sz = size(im);
if any(sz > N)
    hsz = floor(sz / 2);
    s = max(1, hsz);
    t = min(sz, s + N - 1);
    
    im = im(s(1) : t(1), s(2) : t(2), s(3) : t(3));
    sz = size(im);    
end

if any(sz < N)
    lsz = floor((N - sz) / 2);
    rsz = N - sz - lsz;
    im = padarray(im, lsz, 0, 'pre');
    im = padarray(im, rsz, 0, 'post');
end


end

