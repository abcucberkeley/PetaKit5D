function [fsc_mu, res_mu, fsc_cell, res_cell] = XR_one_image_FSC_analysis(fn, varargin)
% perform FSC analysis for a single image
%
% The FSC is between interleaved sub-lattices of the one image (every other
% voxel along each axis, averaged over the four pairs of opposite ones),
% whose noise must be independent. Deconvolution and deskew/rotation (its
% interpolation) correlate the noise of neighboring voxels, so on low-SNR
% data processed that way the sub-lattices share noise and the FSC can keep
% improving with more iterations or a lower Wiener alpha, with no optimum.
% There, compare two independent images instead (XR_two_image_FSC_analysis,
% e.g., two acquisitions), or choose the parameters by cross-validation
% (e.g., with two recorrupted copies of one image).
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
ip.addParameter('halfSize', [501, 501, 501], @isnumeric);
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
halfSize = pr.halfSize;
bbox = pr.bbox;
resAxis = pr.resAxis;
skipConeRegion = pr.skipConeRegion;
clipPer = pr.clipPer;
debug = pr.debug;

% the image of the last file read is kept while the file is unchanged (name,
% modification time and size)
persistent im_orig fn_orig info_orig;
if ischar(fn)
    dir_info = dir(fn);
    info = [dir_info.datenum, dir_info.bytes];
    if isempty(fn_orig) || ~strcmp(fn, fn_orig) || ~isequal(info, info_orig)
        fn_orig = '';
        im_orig = [];
        try
            im_orig = parallelReadTiff(fn);
        catch ME
            mex_fallback('parallelReadTiff', ME, {'tiff:dataTypeError'});
            im_orig = readtiff(fn);
        end
        fn_orig = fn;
        info_orig = info;
    end
    im = im_orig; 
else
    im = fn;
end

% add noise to empty space
if false
    im = double(im);
    im_rand = randn(size(im)) * 4 + 100;
    im = im + im_rand .* (im == 0);
end

if ~isempty(clipPer)
    im = double(im);
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

% resize to isotropic voxels and crop to 2 * halfSize if any size is greater:
% only the kept region is resized (from the input region it depends on), in
% double
sz = size(im, 1 : 3);
psz = min(xyPixelSize, dz);
if dz ~= xyPixelSize
    rsz = round(sz .* [xyPixelSize, xyPixelSize, dz] / psz);
else
    rsz = sz;
end
s = ones(1, 3);
t = rsz;
if any(rsz > halfSize * 2)
    hsz = floor((rsz - halfSize * 2) / 2);
    s = max(1, hsz);
    t = min(rsz, s + halfSize * 2 - 1);
end
if dz ~= xyPixelSize
    im = imresize3_fast(im, rsz, 'linear', [s, t], 'double');
else
    im = double(im(s(1) : t(1), s(2) : t(2), s(3) : t(3)));
end
sz = size(im, 1 : 3);

if any(sz < halfSize * 2)
    lsz = floor((halfSize * 2 - sz) / 2);
    rsz = halfSize * 2 - sz - lsz;
    im = padarray(im, lsz, 0, 'pre');
    im = padarray(im, rsz, 0, 'post');
end

im = reshape(im, 2, halfSize(1), 2, halfSize(2), 2, halfSize(3));
im = permute(im, [2, 4, 6, 1, 3, 5]);

% fsc computing
[fsc_1, res_1] = XR_FSC(im(:, :, :, 1, 1, 1), im(:, :, :, 2, 2, 2), ...
    'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, 'halfSize', halfSize, ...
    'skipConeRegion', skipConeRegion);
[fsc_2, res_2] = XR_FSC(im(:, :, :, 2, 1, 1), im(:, :, :, 1, 2, 2), ...
    'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, 'halfSize', halfSize, ...
    'skipConeRegion', skipConeRegion);    
[fsc_3, res_3] = XR_FSC(im(:, :, :, 1, 2, 1), im(:, :, :, 2, 1, 2), ...
    'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, 'halfSize', halfSize, ...
    'skipConeRegion', skipConeRegion);    
[fsc_4, res_4] = XR_FSC(im(:, :, :, 1, 1, 2), im(:, :, :, 2, 2, 1), ...
    'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, 'halfSize', halfSize, ...
    'skipConeRegion', skipConeRegion);    

fsc_cell = {fsc_1, fsc_2, fsc_3, fsc_4};
res_cell = {res_1, res_2, res_3, res_4};

% resolution computing
fsc_mu = fsc_1;
fsc_mu.fsc = (fsc_1.fsc + fsc_2.fsc + fsc_3.fsc + fsc_4.fsc) / 4;

[res_mu] = XR_FSC_resolution(fsc_mu, 'xyPixelSize', xyPixelSize, 'dz', dz, 'resThreshMethod', resThreshMethod, ...
    'resThresh', resThresh, 'debug', debug);

end


