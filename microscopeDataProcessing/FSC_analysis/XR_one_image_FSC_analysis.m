function [fsc_mu, res_mu, fsc_cell, res_cell] = XR_one_image_FSC_analysis(fn, varargin)
% perform FSC analysis for a single image
% 
% xruan (10/19/2021): add support for bounding box for the region to calculate FSC
% xruan (06/09/2022): add support for clipping very bright spots


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fn', @(x) ischar(x) || isnumeric(x));
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.108, @isnumeric);
ip.addParameter('dr', 1 , @isnumeric);
ip.addParameter('dtheta', pi / 12 , @isnumeric);
ip.addParameter('resThreshMethod', 'fixed', @ischar);
ip.addParameter('resThresh', 0.2, @isnumeric);
ip.addParameter('N', [501, 501, 501], @isnumeric);
ip.addParameter('bbox', [], @isnumeric);
ip.addParameter('resAxis', 'xz', @ischar);
ip.addParameter('skipConeRegion', true, @islogical);
ip.addParameter('clipPer', [], @isnumeric); % clip intensity higher than the given percentile
ip.addParameter('debug', false, @islogical);

ip.parse(fn,  varargin{:});

pr = ip.Results;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
dr = pr.dr;
dtheta = pr.dtheta;
resThreshMethod = pr.resThreshMethod;
resThresh = pr.resThresh;
N = pr.N;
bbox = pr.bbox;
resAxis = pr.resAxis;
skipConeRegion = pr.skipConeRegion;
clipPer = pr.clipPer;
debug = pr.debug;

persistent im_orig fn_orig;
if ischar(fn)
    if isempty(fn_orig) || ~strcmp(fn, fn_orig)
        fn_orig = fn;
        try
            im_orig = parallelReadTiff(fn);
        catch
            im_orig = readtiff(fn);
        end
    end
    im = im_orig; 
else
    im = fn;
end
im = double(im);

% add noise to empty space
if false
    im_rand = randn(size(im)) * 4 + 100;
    im = im + im_rand .* (im == 0);
end

if ~isempty(clipPer)
    clip_val = prctile(im(im > 0), clipPer);
    im = im .* (im <= clip_val) + clip_val .* (im > clip_val);
end

if ~isempty(bbox)
    im = im(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
end

% test permutation of x and y axes
switch resAxis 
    case 'xz'
    case 'yz'
        im = permute(im, [2, 1, 3]);        
    case 'xy'
        im = permute(im, [3, 2, 1]);
end

psz = min(xyPixelSize, dz);
if dz ~= xyPixelSize
    im = imresize3(im, round(size(im) .* [xyPixelSize, xyPixelSize, dz] / psz), 'linear');
end

% crop image if any size is greater than the 2 folder of split images
sz = size(im);
if any(sz > N * 2)
    hsz = floor((sz - N * 2) / 2);
    s = max(1, hsz);
    t = min(sz, s + N * 2 - 1);
    
    im = im(s(1) : t(1), s(2) : t(2), s(3) : t(3));
    sz = size(im);    
end

if any(sz < N * 2)
    lsz = floor((N * 2 - sz) / 2);
    rsz = N * 2 - sz - lsz;
    im = padarray(im, lsz, 0, 'pre');
    im = padarray(im, rsz, 0, 'post');
end

im = reshape(im, 2, N(1), 2, N(2), 2, N(3));
im = permute(im, [2, 4, 6, 1, 3, 5]);

% fsc computing
tic
[fsc_1, res_1] = XR_FSC(im(:, :, :, 1, 1, 1), im(:, :, :, 2, 2, 2), ...
    'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, 'N', N, ...
    'skipConeRegion', skipConeRegion);
[fsc_2, res_2] = XR_FSC(im(:, :, :, 2, 1, 1), im(:, :, :, 1, 2, 2), ...
    'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, 'N', N, ...
    'skipConeRegion', skipConeRegion);    
[fsc_3, res_3] = XR_FSC(im(:, :, :, 1, 2, 1), im(:, :, :, 2, 1, 2), ...
    'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, 'N', N, ...
    'skipConeRegion', skipConeRegion);    
[fsc_4, res_4] = XR_FSC(im(:, :, :, 1, 1, 2), im(:, :, :, 2, 2, 1), ...
    'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, 'N', N, ...
    'skipConeRegion', skipConeRegion);    
toc

fsc_cell = {fsc_1, fsc_2, fsc_3, fsc_4};
res_cell = {res_1, res_2, res_3, res_4};


% resolution computing
fsc_mu = fsc_1;
fsc_mu.fsc = (fsc_1.fsc + fsc_2.fsc + fsc_3.fsc + fsc_4.fsc) / 4;

[res_mu] = XR_FSC_resolution(fsc_mu, 'xyPixelSize', xyPixelSize, 'dz', dz, 'resThreshMethod', resThreshMethod, ...
    'resThresh', resThresh, 'debug', debug);

end


