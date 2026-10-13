function [psf] = deconPSF(sz, sigma, amplitude)
% single precision Gaussian PSF with its peak at the FFT center
% floor(sz / 2) + 1. The first and last 5 z slices are zero, so psf_gen_new
% skips its background removal and the PSF is used as is.
%
%          sz : 1x3 vector (default: [31, 31, 31]). Odd sizes keep the peak at the exact center.
%       sigma : 1x2 vector [sigmaYX, sigmaZ] (default: [1.5, 2]) in voxels.
%   amplitude : a number (default: 10000). Peak value.

if nargin < 1 || isempty(sz)
    sz = [31, 31, 31];
end
if nargin < 2 || isempty(sigma)
    sigma = [1.5, 2];
end
if nargin < 3 || isempty(amplitude)
    amplitude = 10000;
end

psf = pktest.beads(sz, floor(sz / 2) + 1, sigma, amplitude, 0);
psf(:, :, [1 : 5, end - 4 : end]) = 0;

end
