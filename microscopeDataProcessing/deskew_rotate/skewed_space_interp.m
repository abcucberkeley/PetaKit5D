function [im_int] = skewed_space_interp(im, xstep, nint, varargin)
% skewed space interpolation along z axis for given image, with xstep and
% number of slices interpolated (evenly between two z slices). 
%
% Author: Xiongtao Ruan (05/28/2022)



ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('im');
ip.addRequired('xstep'); % relative steps in real space between two z slice
ip.addRequired('nint'); % number of step to interpolation
ip.addParameter('Reverse', true, @islogical);
% ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.parse(im, xstep, nint, varargin{:});


pr = ip.Results;
Reverse = pr.Reverse;

sz = size(im);

im_int = zeros(sz(1), sz(2), (sz(3) - 1) * nint + 1, 'single');

xa = ceil(xstep);

% start index for s and t
s_mat = zeros(nint - 1, 1);
t_mat = zeros(nint - 1, 1);
sw_mat = zeros(nint - 1, 1);
tw_mat = zeros(nint - 1, 1);
zw_mat = zeros(nint - 1, 1);

for i = 2 : nint
    zw = (i - 1) / nint;
    if Reverse
        s = xa - xstep .* zw + 1;
        t = xstep .* (1 - zw) + 1;
    else
        s = xstep .* (1 - zw) + 1;
        t = xa - xstep .* zw + 1;
    end
    
    % distance to start
    sw = s - floor(s);
    s = floor(s);
    tw = t - floor(t);
    t = floor(t);   

    s_mat(i - 1) = s;
    t_mat(i - 1) = t;
    sw_mat(i - 1) = sw;
    tw_mat(i - 1) = tw;
    zw_mat(i - 1) = zw;
end

for z = 1 : sz(3)
    % tic
    if z == sz(3)
        im_int(:, :, (z - 1) * nint + 1) = im(:, :, z);
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

    for i = 1 : nint
        zint = (z - 1) * nint + i;

        % assign the slice from the original image
        if i == 1
            im_int(:, :, zint) = im(:, :, z);
            continue;
        end
        
        s = s_mat(i - 1);
        t = t_mat(i - 1);
        sw = sw_mat(i - 1);
        tw = tw_mat(i - 1);
        zw = zw_mat(i - 1);

        % interpolation
        im_sq = im_s(:, s : s + sz(2) - 1) .* (1 - sw) + im_s(:, s + 1 : s + sz(2)) .* sw;

        im_tq = im_t(:, t : t + sz(2) - 1) .* (1 - tw) + im_t(:, t + 1 : t + sz(2)) .* tw;

        im_q = im_sq * (1 - zw) + im_tq * zw;
                
        im_int(:, :, zint) = im_q;
    end
    % toc
end

end

