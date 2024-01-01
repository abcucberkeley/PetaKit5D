function [fsc, res] = XR_FSC(im1, im2, varargin)
% calculate FSC for different angles in xz 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('im1', @isnumeric);
ip.addRequired('im2', @isnumeric);
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.108, @isnumeric);
ip.addParameter('dr', 1 , @isnumeric);
ip.addParameter('dtheta', pi / 12 , @isnumeric);
ip.addParameter('resThreshMethod', 'fixed', @ischar);
ip.addParameter('resThresh', 0.2, @isnumeric);
ip.addParameter('skipConeRegion', true, @islogical);
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
skipConeRegion = pr.skipConeRegion;
N = pr.N;
debug = pr.debug;

hfN = (N - 1) / 2;

% don't resample image 
% psz = min(xyPixelSize, dz);
% if dz ~= xyPixelSize
%     im1 = imresize3(im1, round(size(im1) .* [xyPixelSize, xyPixelSize, dz] / psz));
%     im2 = imresize3(im2, round(size(im2) .* [xyPixelSize, xyPixelSize, dz] / psz));
% end

if any(size(im1) ~= N)
    sz1 = size(im1);
    im1 = padarray(im1, max(0, floor((N - sz1) / 2)), 0, 'pre');
    im1 = padarray(im1, max(0, ceil((N - sz1) / 2)), 0, 'post');
    
    c_1 = round((size(im1) + 1) / 2);
    im1 = im1(c_1(1) - hfN(1) : c_1(1) + hfN(1), c_1(2) - hfN(2) : c_1(2) + hfN(2), c_1(3) - hfN(3) : c_1(3) + hfN(3));
end

if any(size(im2) ~= N)
    sz2 = size(im2);    
    im2 = padarray(im2, max(0, floor((N - sz2) / 2)), 0, 'pre');
    im2 = padarray(im2, max(0, ceil((N - sz2) / 2)), 0, 'post');
    
    c_2 = round((size(im2) + 1) / 2);
    im2 = im2(c_2(1) - hfN(1) : c_2(1) + hfN(1), c_2(2) - hfN(2) : c_2(2) + hfN(2), c_2(3) - hfN(3) : c_2(3) + hfN(3));
end

% apply hamming windowing
if ~false 
    hwd = hamming(N(1));
    window = (hwd .* reshape(hwd, 1, numel(hwd)) .* reshape(hwd, 1, 1, numel(hwd))) .^ (1 / 3);
    
    im1 = im1 .* window;
    im2 = im2 .* window;
end


% don't caclculate shells every time
persistent azimuth r N_orig

if isempty(N_orig) || N_orig ~= N(1)
    N_orig = N(1);
    [X, Y, Z] = meshgrid(-(N(2) - 1) / 2 : (N(2) - 1) / 2, -(N(1) - 1) / 2 : (N(1) - 1) / 2, -(N(3) - 1) / 2 : (N(3) - 1) / 2);

    [azimuth,elevation,r] = cart2sph(X, Y, Z);
    % use y axis as the rotation axis
    azimuth = permute(azimuth, [3, 1, 2]);
    azimuth = pi * (azimuth < 0) + azimuth;
end

% generate shells for 
shells_r = round((r + dr / 2) / dr);
% shells_e = round((elevation + pi / 2 + pi / 24) / (pi / 12));
shells_a = round((azimuth) / dtheta) + 1;
% merge disconnected regions
shells_a(shells_a == max(shells_a(:))) = min(shells_a(:));
% uniq_re = unique([shells_r(:), shells_e(:)], 'row');
% uniq_re(uniq_re(:, 1) > (N(1) - 1) / 2 + 1, :) = []; 

shells_r((shells_r - 1) * dr > (N(1) - 1) / 2 + 1) = round(((N(1) - 1) / 2 + 1) / dr) + 1;
shells_a(shells_r == max(shells_r(:))) = max(shells_a(:)) + 1;

max_r = max(shells_r(:));
max_e = max(shells_a(:));
shells = sub2ind([max_r, max_e], shells_r, shells_a);

% remove a section along z-axis
if ~false && skipConeRegion
    shells(azimuth < (pi / 36) | azimuth > pi - (pi / 36) ) = max_r * max_e;
    % shells(abs(azimuth - pi / 2) < pi / 96)  = max_r * max_e;
end

% remove a section along x-axis
if ~false && skipConeRegion
    % shells(azimuth < (pi / 36) | azimuth > pi - (pi / 36) ) = max_r * max_e;
    shells(abs(azimuth - pi / 2) < pi / 48)  = max_r * max_e;
end


num_xi = max_r * max_e;

% Compute Fourier transforms
im1_fft = fftshift(fftn(im1));
im2_fft = fftshift(fftn(im2));

% a = im1_fft .* conj(im2_fft);
% b = abs(im1_fft) .^ 2;
% c = abs(im2_fft) .^ 2;
frc = real(accumarray(shells(:), im1_fft(:) .* conj(im2_fft(:)), [num_xi, 1]) ) ...
        ./ sqrt(accumarray(shells(:), abs(im1_fft(:)).^2, [num_xi, 1]) .* accumarray(shells(:), abs(im2_fft(:)).^2, [num_xi, 1]));

N_points = accumarray(shells(:), ones(numel(shells), 1), [num_xi, 1]);

xi = [1 : max_r * max_e]';
xi(isnan(frc)) = [];
xi(end) = [];
N_points(isnan(frc)) = [];
N_points(end) = [];
frc(isnan(frc)) = [];
frc(end) = [];

[xri, xei] = ind2sub([max_r, max_e], xi);

% frc for each angle
unique_e = unique(xei);
frc_cell = cell(numel(unique_e), 1);
r_cell = cell(numel(unique_e), 1);
npoints_mat = zeros(max_r, max_e);
frc_mat = zeros(max_r, max_e);
for i = 1 : numel(unique_e)
    inds = xei == unique_e(i);
    frc_i = frc(inds);
    r_i = xri(inds);    
    
    npoints_mat(r_i, i) = N_points(inds);
    
    frc_cell{i} = frc_i;
    frc_mat(r_i, i) = frc_i;
    r_cell{i} = r_i;
end
npoints_mat(end, :) = [];
npoints_mat(:, end) = [];

frc_mat(end , :) = [];
frc_mat(: , end) = [];
% set missing value for r = 0 for some angles
frc_mat(1, frc_mat(1, :) == 0) = 1;

% Number of Fourier frequencies in each cell
% n_xi = accumarray(shells(:), ones([numel(shells),1]));

fsc.r = ((1 : (max_r-1))' - 1) * dr / ((N(1) - 1) / 2);
fsc.thetas = (unique_e - 1) * dtheta;
fsc.fsc = frc_mat;
fsc.npoints = npoints_mat;


% calculate resulotion
[res] = XR_FSC_resolution(fsc, 'xyPixelSize', xyPixelSize, 'dz', dz, 'resThreshMethod', resThreshMethod, ...
    'resThresh', resThresh, 'debug', debug);


end





