function [] = XR_psf_detection_and_cropping(fn, resultDir, varargin)
% check psf and crop them from the raw calibration image. 
%
% The idea is to detect local maximum in radon image for the skewed angle in xz to
% get the peaks of psfs. Then, remove peaks if they are close to each other
% or to the boarder. Finally, crop the psfs for kept peaks for given crop
% size. 
%
% xruan (06/11/2021): also save peak info and as placehoder
% xruan (07/15/2021): discard saturated psfs

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fn');
ip.addRequired('resultDir'); 
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric);
ip.addParameter('angle', 32.45, @isnumeric);
ip.addParameter('cropSize', [256, 128, 201], @isnumeric);
ip.addParameter('distThresh', [256, 128, 201], @isnumeric);
ip.addParameter('prefix', 'test_', @ischar);
% ip.addParameter('prefix', 'test_', @ischar);
ip.parse(fn, resultDir, varargin{:});

pr = ip.Results;
dz = pr.dz;
pixelSize = pr.xyPixelSize;
angle = pr.angle;
cropSize = pr.cropSize;
distThresh = pr.distThresh;
prefix = pr.prefix;

tic
info_fn = [resultDir, prefix, '.mat'];
if exist(info_fn, 'file')
    return;
end

if ischar(fn)
    fn = {fn};
end
fprintf('Start PSF detection for %s ...\n', fn{1});

% read images (merge parts if the image is in serval partial images). 
ims = cell(numel(fn), 1);
for f = 1 : numel(fn)
    ims{f} = readtiff(fn{f});
end
        
im = cat(3, ims{:});
clear ims
sz = size(im);
im(isnan(im)) = 0;

%% radon transform for the given angle

theta = angle;
theta_1 = atand(tand(theta) * pixelSize / (dz * sind(theta)));
radon_y = cell(sz(1), 1);
for y = 1 : sz(1)
    radon_y{y} = radon(squeeze(im(y, :, :)), theta_1);
end

im_radon = cat(2, radon_y{:})';


%%
% local maximum
BW = imregionalmax(im_radon);

[py, px] = find(BW);

% use local thresholding to keep detected peaks
hbbox = [5, 5];
std_thrsh = 2;
kept_flag = false(numel(py), 1);
im_std = std(im_radon(:));
for i = 1 : numel(py)
    yi = max(1, py(i) - hbbox(1)) : min(size(im_radon, 1), py(i) + hbbox(1));
    xi = max(1, px(i) - hbbox(2)) : min(size(im_radon, 2), px(i) + hbbox(2));
    
    im_radon_i = im_radon(yi, xi);
    
    peak_i = im_radon(py(i), px(i));
    
    bg_i = median(im_radon_i(:));
    std_i = std(im_radon_i(:));
    
    if peak_i - bg_i > std_thrsh * (std_i * 0.95 + im_std * .05)
        kept_flag(i) = true;
    end
end

if false
    figure, imshow(im_radon, [])
    hold on, plot(px(kept_flag), py(kept_flag), 'o')
end

%% decide the center for each detected psf
im_i = double(squeeze(im(1, :, :)));
[rd, Xp] = radon(im_i, theta_1);

pyc = py(kept_flag);
pxc = px(kept_flag);

pcoords = zeros(numel(pyc), 3);
k = tand(theta_1 + 90);
[Z, X] = ndgrid(1 : size(im_i, 1), 1 : size(im_i, 2));
rc = floor((size(im_i)+1)/2);
Z = flip(Z, 1);
Z = Z - rc(1);
X = X - rc(2);

for i = 1 : numel(pyc)
    peak_y = pyc(i);
    im_i = double(squeeze(im(peak_y, :, :)));
    
    c = Xp(pxc(i)) * sind(theta_1) - Xp(pxc(i)) * cosd(theta_1) * k;
    % c = (size(im_i, 1));
    
    im_line = abs(X * k + c - Z) < 10;
    
    [~, peak_ind] = max(reshape(im_line .* im_i, [], 1));
    [peak_x, peak_z] = ind2sub(size(im_i), peak_ind);
    pcoords(i, :) = [peak_y, peak_x, peak_z];
end


%% search to exclude nearby psfs

% dist_thrsh = [256, 150, 201];
% dist_thrsh = [128, 200, 251];
dist_thrsh = distThresh;

include_flag = false(size(pcoords, 1), 1);
for i = 1 : size(pcoords, 1)
    pc_i = pcoords(i, :);
    
    close_i = all(abs(pc_i - pcoords) <= dist_thrsh, 2);
    if sum(close_i) == 1
        include_flag(i) = true;
    end
end
pcoords = pcoords(include_flag, :);


% exclude peaks close to the border
border_thrsh = (cropSize - 1) / 2;

include_flag = all(pcoords >= border_thrsh & size(im) - pcoords >= border_thrsh, 2);

pcoords = pcoords(include_flag, :);


%% crop picked psfs
% cropSize = [256, 128, 251];
% cropSize = [128, 200, 251];
hsize = round((cropSize - 1) / 2);

result_dir = resultDir;
mkdir(result_dir);
mip_dir = [result_dir, 'MIPs/'];
mkdir(mip_dir);
include_flag = true(size(pcoords, 1), 1);
for i = 1 : size(pcoords, 1)
    pc_i = pcoords(i, :);

    bbox = [pc_i(1) - hsize(1), pc_i(2) - hsize(2), pc_i(3) - hsize(3), pc_i(1) + hsize(1), pc_i(2) + hsize(2), pc_i(3) + hsize(3)];
    im_i = im(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
    
    % not include saturated (or close to saturation) psfs
    if max(im_i(:)) > 60000
        include_flag(i) = false;
        continue;
    end

    fn_i = sprintf('%s/%s_crop_%d_%d_%d.tif', result_dir, prefix, pc_i(1), pc_i(2), pc_i(3));
    writetiff(im_i, fn_i);
    
    mip_fn_i = sprintf('%s/%s_crop_%d_%d_%d_MIP_z.tif', mip_dir, prefix, pc_i(1), pc_i(2), pc_i(3));
    writetiff(max(im_i, [], 3), mip_fn_i);
end

pcoords = pcoords(include_flag, :);

save('-v7.3', info_fn, 'pcoords', 'pr');

fprintf('Done!\n');
toc

end

