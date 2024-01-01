function [b_omw, OTF_bp_omw, abs_OTF_c, OTF_mask] = omw_backprojector_generation(psf, alpha, skewed, varargin)
% generate OTF masked wiener back projector
%
% Author: Xiongtao Ruan (11/10/2022)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('psf'); 
ip.addRequired('alpha'); 
ip.addRequired('skewed'); 
ip.addParameter('OTFCumThresh', 0.9, @isnumeric);
ip.addParameter('OTFAreaThresh', 100, @isnumeric);
ip.addParameter('hanWinBounds', [0.8, 1], @(x) isnumeric(x) && numel(x) == 2);
ip.addParameter('minIntThrsh', 2.5e-3, @(x) isnumeric(x));

ip.parse(psf, alpha, skewed, varargin{:});

pr = ip.Results;
otf_thresh = pr.OTFCumThresh;
area_thresh = pr.OTFAreaThresh;
han_bounds = pr.hanWinBounds;
minIntThrsh = pr.minIntThrsh;

if isempty(psf) || all(psf == 0, 'all')
    error('A valid nonzero PSF needs to be provided!')
end

psf = single(psf);
psf = psf ./ sum(psf, 'all');

OTF = decon_psf2otf(psf);
abs_OTF = abs(OTF);
abs_OTF_c = fftshift(abs_OTF);
OTF_vals = sort(double(abs_OTF(:)), 'descend');

% segment out OTF support based on the accumulate OTF intensity
tind = find(cumsum(OTF_vals) > sum(OTF_vals) * otf_thresh, 1, 'first');
dc_thresh = OTF_vals(tind) / abs_OTF(1, 1, 1);
% (11/12/2022): add a lower bound 2.5e-3 in case of too noisy
dc_thresh = max(dc_thresh, minIntThrsh);
fprintf('OTF threshold: %0.2d\n', dc_thresh);
OTF_mask = abs_OTF_c > abs_OTF(1, 1, 1) * dc_thresh;
OTF_mask = imopen(OTF_mask, strel('sphere', 2));
OTF_mask = bwareaopen(OTF_mask, area_thresh);
OTF_mask = imclose(OTF_mask, strel('sphere', 2));
OTF_mask = imopen(OTF_mask, strel('sphere', 2));
OTF_mask = bwareaopen(OTF_mask, area_thresh);

if ~any(OTF_mask, 'all')
    error('The OTF mask is empty, check the PSF and OTF-related parameters!')
end

% for skewed space data, the OTF mask has three main componets, need to
% concatenate them along z. 
% first automatically decide if the OTF is in skewed space
if ~false && isempty(skewed)
    skewed = false;
    OTF_mask_xz = squeeze(sum(OTF_mask, 1)) > 0;
    CC = bwconncomp(OTF_mask_xz);
    % if there are more than one connected components, check if the x projected 
    % line cover the whole z range, and also if peak is at/close to center.
    if CC.NumObjects > 1
        % OTF_mask_xy_line = squeeze(sum(OTF_mask, [1, 2]));
        OTF_mask_xy_line = sum(OTF_mask_xz, 1);
        [~, pind] = max(OTF_mask_xy_line);
        if all(OTF_mask_xy_line >0) && abs(pind - (numel(OTF_mask_xy_line) + 1) / 2) > 1
            skewed = true;
        end
    end
end

if skewed
    L = bwlabeln(OTF_mask);
    if numel(unique(L)) > 4
        CC = bwconncomp(OTF_mask, 26);
        vols = cellfun(@numel, CC.PixelIdxList);
        [~, max_inds] = maxk(vols, 3);
        OTF_mask = false(size(OTF_mask));
        OTF_mask(cat(1, CC.PixelIdxList{max_inds})) = true;
        L = bwlabeln(OTF_mask);
    elseif numel(unique(L)) == 2
        L = (L + 1) .* OTF_mask;
    end
    OTF_mask_c = L == 2;
    OTF_mask_l = L == 1;
    OTF_mask_r = L == 3;
    OTF_mask = cat(3, OTF_mask_r, OTF_mask_c, OTF_mask_l);
else
    CC = bwconncomp(OTF_mask, 26);
    [~, max_ind] = max(cellfun(@numel, CC.PixelIdxList));
    OTF_mask = false(size(OTF_mask));
    OTF_mask(CC.PixelIdxList{max_ind}) = true;
end
psz = size(OTF_mask, 1 : 3);

% get convex hull for the OTF mask
[zinds] = find(sum(OTF_mask, [1, 2]))';
[xinds] = find(sum(OTF_mask, [1, 3]));
[yinds] = find(sum(OTF_mask, [2, 3]))';
for z = zinds
    OTF_mask(:, :, z) = bwconvhull(OTF_mask(:, :, z), "union");
end
for y = yinds
    OTF_mask(y, :, :) = bwconvhull(squeeze(OTF_mask(y, :, :)), "union");
end
for x = xinds
    OTF_mask(:, x, :) = bwconvhull(squeeze(OTF_mask(:, x, :)), "union");
end

% direct compute 3D convex hull image, it is slower, so disable it for now.
if false
    [c, r, p] = ndgrid(1 : psz(1), 1 : psz(2), 1 : psz(3));
    stats = regionprops3(OTF_mask, 'ConvexHull');
    hull = stats.ConvexHull{1};
    dt = delaunayTriangulation(hull(:, [2, 1, 3]));
    idx = pointLocation(dt, c(:), r(:), p(:));
    OTF_mask = reshape(~isnan(idx), psz);    
end

% use relative star distance to the edge to define the distance
OTF_mask_edge = OTF_mask & (~imerode(OTF_mask, strel('sphere', 1)));
c = (psz + 1) / 2;
[Y, X, Z] = ndgrid(1 : psz(1), 1 : psz(2), 1 : psz(3));
Y = Y - c(1);
X = X - c(2);
Z = Z - c(3);

[thetas, phis, rs] = cart2sph(X, Y, Z);

% normalized distances 
[Xn, Yn, Zn] = sph2cart(thetas, phis, ones(psz));

% use knn to find the nearest border points, and take average for r (k = 2 is good enough). 
[inds, d_mat] = knnsearch([Yn(OTF_mask_edge), Xn(OTF_mask_edge), Zn(OTF_mask_edge)], [Yn(:), Xn(:), Zn(:)], 'K', 2);
Rn = rs(OTF_mask_edge);
d_1 = d_mat(inds(:, 1), 1);
d_2 = d_mat(inds(:, 2), 2);
rn = (Rn(inds(:, 1)) .* d_2 + Rn(inds(:, 2)) .* d_1) ./ (d_1 + d_2); 
rn = reshape(rn, psz);
bw_dist = rs ./ rn;
clear X Y Z thetas phis rs Xn Yn Zn Rn rn;

% apodization function: using cosine square function
l = han_bounds(1);
u = han_bounds(2);
win_func = @(x) (x < l) + cos(pi * (x - l)/ 2 / (u - l)) .^ 2 .* (x >= l & x < u);
mask = double(win_func(bw_dist));

if skewed
    mask = sum(reshape(mask, psz(1), psz(2), psz(3) / 3, 3), 4);
    OTF_mask = sum(reshape(OTF_mask, psz(1), psz(2), psz(3) / 3, 3), 4);
end

OTF_bp_om = ifftshift(mask);

% generate wiener back projector
OTF_bp_w = conj(OTF) ./(abs_OTF.^2 + alpha); % Wiener filter

% OTF mask and wiener
OTF_bp_omw = OTF_bp_om .* OTF_bp_w;

b_omw = fftshift(real(ifftn(OTF_bp_omw)));


end

