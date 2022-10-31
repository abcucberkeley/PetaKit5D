function [im_int] = skewed_space_interp_defined_stepsize(im, xstep, stepsize, varargin)
% skewed space interpolation along z axis for given image, such that the result
% has defined evenly spaced number of z slices instead of interpolate from the n-th even points
% between two slices. stepsize <= 1
%
% Author: Xiongtao Ruan (10/28/2022)



ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('im');
ip.addRequired('xstep'); % relative steps in real space between two z slice
ip.addRequired('stepsize'); % sampling z step size 
ip.addParameter('Reverse', true, @islogical);
% ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.parse(im, xstep, stepsize, varargin{:});


pr = ip.Results;
Reverse = pr.Reverse;

sz = size(im);
if stepsize == 1
    im_int = im; 
    return;
end

xa = ceil(xstep);

% start index for s and t
zout = 1 : stepsize : sz(3);
nint = numel(zout);
zs_mat = floor(zout);
im_int = zeros(sz(1), sz(2), nint, 'single');

if Reverse
    sz_mat = xa - (zout - zs_mat) * xstep + 1;
    tz_mat = sz_mat - (xa - xstep);    
else
    sz_mat = (zout - zs_mat) * xstep + 1;
    tz_mat = sz_mat + (xa - xstep);    
end

s_mat = floor(sz_mat);
t_mat = floor(tz_mat);
sw_mat = sz_mat - s_mat;
tw_mat = tz_mat - t_mat;
zw_mat = zout - zs_mat;

% change to two layers of loops with outer loop as original z
% tic
for z = 1 : sz(3)
    % tic
    if z == sz(3)
        im_int(:, :, end) = im(:, :, z);
        break;
    end

    % pad images
    im_s = zeros(sz(1), sz(2) + xa, 'single');
    im_t = zeros(sz(1), sz(2) + xa, 'single');
    
    if Reverse
        im_s(:, xa + 1 : xa + sz(2)) = im(:, :, z);
        im_t(:, 1 : sz(2)) = im(:, :, z + 1);
    else
        im_s(:, 1 : sz(2)) = im(:, :, z);
        im_t(:, xa + 1 : xa + sz(2)) = im(:, :, z + 1);
    end

    zs_z = find(zs_mat == z);

    for i = 1 : numel(zs_z)
        zint = zs_z(i);

        s = s_mat(zint);
        t = t_mat(zint);
        sw = sw_mat(zint);
        tw = tw_mat(zint);
        zw = zw_mat(zint);

        % assign the slice from the original image
        if sw == 0
            im_int(:, :, zint) = im(:, :, z);
            continue;
        end
        
        % interpolation
        im_sq = im_s(:, s : s + sz(2) - 1) .* (1 - sw) + im_s(:, s + 1 : s + sz(2)) .* sw;

        im_tq = im_t(:, t : t + sz(2) - 1) .* (1 - tw) + im_t(:, t + 1 : t + sz(2)) .* tw;

        % im_q = (im_sq + im_tq) / 2;
        im_q = im_sq * (1 - zw) + im_tq * zw;
                
        im_int(:, :, zint) = im_q;
    end
    % toc
end
% toc

end

