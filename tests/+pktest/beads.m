function [vol] = beads(sz, positions, sigma, amplitude, background)
% generate a single precision volume of Gaussian spots
%
%          sz : 1x3 vector. Volume size [ny, nx, nz].
%   positions : #spots x 3. Spot centers [y, x, z], 1-based, can be fractional.
%       sigma : scalar or 1x2 vector [sigmaYX, sigmaZ] in voxels.
%   amplitude : scalar or #spots x 1 (default: 1000). Peak heights.
%  background : scalar (default: 0). Constant offset.

if nargin < 4 || isempty(amplitude)
    amplitude = 1000;
end
if nargin < 5 || isempty(background)
    background = 0;
end
if isscalar(sigma)
    sigma = [sigma, sigma];
end
amplitude = amplitude(:) .* ones(size(positions, 1), 1);

[y, x, z] = ndgrid(1 : sz(1), 1 : sz(2), 1 : sz(3));
vol = zeros(sz, 'single') + background;
for i = 1 : size(positions, 1)
    p = positions(i, :);
    vol = vol + single(amplitude(i) * exp(-((y - p(1)) .^ 2 + (x - p(2)) .^ 2) / (2 * sigma(1) ^ 2) ...
        - (z - p(3)) .^ 2 / (2 * sigma(2) ^ 2)));
end

end
