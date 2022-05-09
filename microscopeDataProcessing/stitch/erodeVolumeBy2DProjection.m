function [volout] = erodeVolumeBy2DProjection(vol, esize, varargin)
% Edge erosion by erosion of the max projection in xz, instead of erode the
% whole 3D volume. This is more efficient for DS and DSR-like volume
% 
% Author: Xiongtao Ruan (10/11/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isnumeric);
ip.addRequired('esize', @isnumeric); % erode size scalar.
ip.parse(vol, esize, varargin{:});

sz = size(vol);

if ismatrix(vol)
    mask = vol ~= 0;
    volout = vol .* cast(mask, class(vol)); 
    return;
end

% caculate xz projection and erode the projection
MIP = squeeze(max(vol, [], 1) > 0);
MIP_pad = false(sz([2, 3]) + 2);
MIP_pad(2 : end - 1, 2 : end - 1) = MIP;
MIP_pad = imerode(MIP_pad, strel('square', esize * 2 + 1));
MIP = MIP_pad(2 : end - 1, 2 : end - 1);

% create 3D mask and also erode in y-axis
mask = false(sz);
mask(esize + 1 : end - esize, :, :) = repmat(permute(MIP, [3, 1, 2]), sz(1) - 2 * esize, 1, 1);
volout = vol .* cast(mask, class(vol)); 

end