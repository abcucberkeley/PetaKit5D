%[A_est, c_est] = estGaussianAmplitude3D(vol, sigma, varargin) calculates the 
% amplitude and background coefficient for Gaussians centered on all voxels of
% the input volume.
%
% INPUTS
%     vol : input volume
%   sigma : standard deviation of the Gaussian PSF
%           If the PSF is anisotropic, 'sigma' should be a two-element 
%           vector: [sigma_xy sigma_z]
%
% OPTIONS (as 'specifier'-value pairs):
%    'WindowSize' : Window size for the fit, in pixels. Default: 2*sigma,
%                   i.e., [-2*sigma ... 2*sigma] in each dimension
%
% OUTPUTS
%   A_est : volume of estimated amplitudes
%   c_est : volume of estimated background coefficients

% Francois Aguet, 08/2013 (last modified 10/27/2013)
% Xiongtao Ruan, July 2019 use matlab functions for convolution, which
% turns out to be much faster. 

function [A_est, c_est] = estGaussianAmplitude3D(vol, sigma, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isnumeric);
ip.addRequired('sigma', @isnumeric);
ip.addParamValue('WindowSize', []);
ip.parse(vol, sigma, varargin{:});

if numel(sigma)==1
    sigma = [sigma sigma];
end

if isempty(ip.Results.WindowSize)
    wx = ceil(2*sigma(1));
    wz = ceil(2*sigma(2));
else
    wx = ip.Results.WindowSize(1);
    if numel(ip.Results.WindowSize)==1
        wz = wx;
    else
        wz = ip.Results.WindowSize(2);
    end
end

% gx = exp(-(0:wx).^2/(2*sigma(1)^2));
% gz = exp(-(0:wz).^2/(2*sigma(2)^2));
% fg = conv3fast(vol, gx, gx, gz);
% fu = conv3fast(vol, ones(1,wx+1), ones(1,wx+1), ones(1,wz+1));
% xruan 07/26/2019 use Matlab convolution functions
gauss_kernel = fspecial3('gauss', [wx, wx, wz] * 2 + 1, [sigma(1), sigma(1), sigma(2)]);
fg = imfilter(vol, gauss_kernel, 'conv', 'same', 'symmetric');

sum_kernel = ones(wx * 2 + 1, wx * 2 + 1, wz * 2 + 1);
fu =  imfilter(vol, sum_kernel, 'conv', 'same', 'symmetric');

% Gaussian kernel (spatial)
[x,y,z] = meshgrid(-wx:wx,-wx:wx,-wz:wz);
g = exp(-(x.^2+y.^2)/(2*sigma(1)^2)) .* exp(-z.^2/(2*sigma(2)^2));
n = numel(g);
gsum = sum(g(:));
g2sum = sum(g(:).^2);

% xruan 07/25/2019 adapt to the old method, because of no normalization of
% the convolution
fg = fg * gsum;

% solution to linear system
A_est = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
c_est = (fu - A_est*gsum)/n;
