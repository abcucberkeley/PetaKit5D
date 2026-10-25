function [b_omw, OTF_bp_omw, abs_OTF_c, OTF_mask] = omw_backprojector_generation(psf, alpha, skewed, varargin)
% generate OTF masked wiener back projector
%
% Author: Xiongtao Ruan (11/10/2022)
%
% With omw_backprojector_mex, the sorted OTF values, the connected components,
% the convex hulls of the OTF mask slices and the window are computed with the
% same results as the MATLAB code (the fallback without the mex): the hulls
% with bwconvhull's algorithm in parallel, the window from the nearest edge
% points only where bounds of the edge radii around a voxel's direction don't
% decide it.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('psf'); 
ip.addRequired('alpha'); 
ip.addRequired('skewed'); 
ip.addParameter('OTFCumThresh', 0.9, @isnumeric);
ip.addParameter('OTFAreaThresh', 100, @isnumeric);
ip.addParameter('hannWinBounds', [0.8, 1], @(x) isnumeric(x) && numel(x) == 2);
ip.addParameter('minIntThrsh', 1e-3, @(x) isnumeric(x));

ip.parse(psf, alpha, skewed, varargin{:});

pr = ip.Results;
otf_thresh = pr.OTFCumThresh;
area_thresh = pr.OTFAreaThresh;
han_bounds = pr.hannWinBounds;
minIntThrsh = pr.minIntThrsh;

if isempty(psf) || all(psf == 0, 'all')
    error('A valid nonzero PSF needs to be provided!')
end

psf = single(psf);
% summed in double, so that the back projector doesn't depend on the number
% of threads (a single precision sum is split among them)
psf = psf ./ sum(psf, 'all', 'double');

OTF = decon_psf2otf(psf);
abs_OTF = abs(OTF);
abs_OTF_c = fftshift(abs_OTF);
OTF_vals = sort_descend(abs_OTF);

% segment out OTF support based on the accumulate OTF intensity
tind = find(cumsum(OTF_vals) > sum(OTF_vals) * otf_thresh, 1, 'first');
dc_thresh = OTF_vals(tind) / abs_OTF(1, 1, 1);
% (11/12/2022): add a lower bound 2.5e-3 in case of too noisy
dc_thresh = max(dc_thresh, minIntThrsh);
fprintf('OTF threshold: %0.2d\n', dc_thresh);
OTF_mask = abs_OTF_c > abs_OTF(1, 1, 1) * dc_thresh;
OTF_mask = imopen(OTF_mask, strel('sphere', 2));
OTF_mask = area_open(OTF_mask, area_thresh);
OTF_mask = imclose(OTF_mask, strel('sphere', 2));
OTF_mask = imopen(OTF_mask, strel('sphere', 2));
OTF_mask = area_open(OTF_mask, area_thresh);

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
    % the number of distinct labels (bwlabeln numbers the components 1 to max)
    nL = max(L(:)) + any(L(:) == 0);
    if nL > 4
        OTF_mask = largest_components(OTF_mask, 3);
        L = bwlabeln(OTF_mask);
    elseif nL == 2
        L = (L + 1) .* OTF_mask;
    end
    OTF_mask_c = L == 2;
    OTF_mask_l = L == 1;
    OTF_mask_r = L == 3;
    OTF_mask = cat(3, OTF_mask_r, OTF_mask_c, OTF_mask_l);
else
    OTF_mask = largest_components(OTF_mask, 1);
end
psz = size(OTF_mask, 1 : 3);

% get convex hull for the OTF mask
[zinds] = find(sum(OTF_mask, [1, 2]))';
[xinds] = find(sum(OTF_mask, [1, 3]));
[yinds] = find(sum(OTF_mask, [2, 3]))';
if exist('omw_backprojector_mex', 'file') == 3 && all(psz > 1)
    % the same hulls as bwconvhull, the slices of each loop in parallel
    OTF_mask = omw_backprojector_mex('convhull', OTF_mask, zinds, yinds, xinds, maxNumCompThreads);
else
    for z = zinds
        OTF_mask(:, :, z) = bwconvhull(OTF_mask(:, :, z), "union");
    end
    for y = yinds
        OTF_mask(y, :, :) = bwconvhull(squeeze(OTF_mask(y, :, :)), "union");
    end
    for x = xinds
        OTF_mask(:, x, :) = bwconvhull(squeeze(OTF_mask(:, x, :)), "union");
    end
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

% apodization function: using cosine square function
l = han_bounds(1);
u = han_bounds(2);
win_func = @(x) (x < l) + cos(pi * (x - l)/ 2 / (u - l)) .^ 2 .* (x >= l & x < u);

mask = omw_window_with_mex(OTF_mask_edge, psz, han_bounds, win_func);
if isempty(mask)
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
    % (d_mat is indexed with edge point indices, so the weights are the
    % distances of the first voxels; kept as it is)
    d_1 = d_mat(inds(:, 1), 1);
    d_2 = d_mat(inds(:, 2), 2);
    rn = (Rn(inds(:, 1)) .* d_2 + Rn(inds(:, 2)) .* d_1) ./ (d_1 + d_2); 
    rn = reshape(rn, psz);
    bw_dist = rs ./ rn;
    clear X Y Z thetas phis rs Xn Yn Zn Rn rn;

    mask = double(win_func(bw_dist));
end

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


function [vals] = sort_descend(a)
% sort(double(a(:)), 'descend'), in the mex for single values that are not NaN
% or negative

vals = [];
if exist('omw_backprojector_mex', 'file') == 3 && isa(a, 'single') && isreal(a)
    vals = omw_backprojector_mex('sortdesc', a, maxNumCompThreads);
end
if isempty(vals)
    vals = sort(double(a(:)), 'descend');
end

end


function [bw] = area_open(bw, p)
% bwareaopen(bw, p), in the mex for 3D masks

if exist('omw_backprojector_mex', 'file') == 3 && islogical(bw) && ndims(bw) == 3 && isa(p, 'double') ...
        && isscalar(p) && isreal(p) && p >= 0 && p == round(p) && p <= 9e15
    bw = omw_backprojector_mex('components', bw, 'minsize', p, maxNumCompThreads);
else
    bw = bwareaopen(bw, p);
end

end


function [bw] = largest_components(bw, k)
% the k largest 26-connected components (ties to the first one in
% bwconncomp's order, as max and maxk choose), in the mex for 3D masks

if exist('omw_backprojector_mex', 'file') == 3 && islogical(bw) && ndims(bw) == 3
    bw = omw_backprojector_mex('components', bw, 'largest', k, maxNumCompThreads);
    return;
end
CC = bwconncomp(bw, 26);
vols = cellfun(@numel, CC.PixelIdxList);
if k == 1
    [~, max_inds] = max(vols);
else
    [~, max_inds] = maxk(vols, k);
end
bw = false(size(bw));
bw(cat(1, CC.PixelIdxList{max_inds})) = true;

end


function [mask] = omw_window_with_mex(OTF_mask_edge, psz, han_bounds, win_func)
% the window of the code above with omw_backprojector_mex: it is 1 or 0 where
% the radius of a voxel is below or above the bounds of the edge radii around
% its direction, and computed as above for the other voxels, with the
% distances of the first voxels that weigh the radii. Empty if the mex can't
% be used, or the edge points or distances could make the window NaN.

mask = [];
if exist('omw_backprojector_mex', 'file') ~= 3 || ~isa(han_bounds, 'double') || ~isreal(han_bounds) ...
        || ~(han_bounds(1) > 0 && han_bounds(1) < han_bounds(2) && han_bounds(2) < Inf)
    return;
end
c = (psz + 1) / 2;
eInds = find(OTF_mask_edge);
nE = numel(eInds);
if nE < 2
    return;
end
[ey, ex, ez] = ind2sub(psz, eInds);
[thetas, phis, Rn] = cart2sph(ex - c(2), ey - c(1), ez - c(3));
[Xn, Yn, Zn] = sph2cart(thetas, phis, ones(nE, 1));
E = [Yn, Xn, Zn];
[mask, uInds] = omw_backprojector_mex('window', psz, E, Rn, han_bounds, maxNumCompThreads);

% nearest edge points of the first nE voxels and the voxels to compute
qInds = [(1 : nE)'; uInds];
[qy, qx, qz] = ind2sub(psz, qInds);
[thetas, phis, rs] = cart2sph(qx - c(2), qy - c(1), qz - c(3));
[Xn, Yn, Zn] = sph2cart(thetas, phis, ones(numel(qInds), 1));
[inds, d_mat] = knnsearch(E, [Yn, Xn, Zn], 'K', 2);
d_mat = d_mat(1 : nE, :);
% the bounds assume positive weights and radii
if ~(min(d_mat(:, 2)) > 0 && min(Rn) > 0)
    mask = [];
    return;
end
inds = inds(nE + 1 : end, :);
d_1 = d_mat(inds(:, 1), 1);
d_2 = d_mat(inds(:, 2), 2);
rn = (Rn(inds(:, 1)) .* d_2 + Rn(inds(:, 2)) .* d_1) ./ (d_1 + d_2);
mask(uInds) = win_func(rs(nE + 1 : end) ./ rn);

end

