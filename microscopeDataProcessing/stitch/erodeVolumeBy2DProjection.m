function [vol] = erodeVolumeBy2DProjection(vol, esize, varargin)
% Edge erosion by erosion of the max projection in xz, instead of erode the
% whole 3D volume. This is more efficient for DS and DSR-like volume
% 
% Author: Xiongtao Ruan (10/11/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isnumeric);
ip.addRequired('esize', @isnumeric); % erode size scalar.
ip.parse(vol, esize, varargin{:});

if esize == 0
    return;
end
sz = size(vol);

if ismatrix(vol)
    mask = vol ~= 0;
    vol = vol .* cast(mask, class(vol)); 
    return;
end

% caculate xz projection and erode the projection
MIP = squeeze(max(vol, [], 1) > 0);
if all(MIP > 0, 'all')
    MIP_erode = false(size(MIP));
    MIP_erode(esize + 1 : end - esize, esize + 1 : end - esize) = MIP(esize + 1 : end - esize, esize + 1 : end - esize);
else
    MIP_pad = false(sz([2, 3]) + 2);
    MIP_pad(2 : end - 1, 2 : end - 1) = MIP;
    MIP_pad = imerode(MIP_pad, strel('square', esize * 2 + 1));
    MIP_erode = MIP_pad(2 : end - 1, 2 : end - 1);
end

% create 3D mask and also erode in y-axis
% mask = false(sz);
% mask(esize + 1 : end - esize, :, :) = repmat(permute(MIP, [3, 1, 2]), sz(1) - 2 * esize, 1, 1);
% volout = vol .* cast(mask, class(vol)); 
vol = vol .* cast(permute(MIP_erode, [3, 1, 2]), class(vol));
try 
    indexing3d_mex(vol, [1, 1, 1, esize, sz(2), sz(3)], zeros(esize, sz(2), sz(3), class(vol)));
    indexing3d_mex(vol, [sz(1) - esize + 1, 1, 1, sz(1), sz(2), sz(3)], zeros(esize, sz(2), sz(3), class(vol)));
catch ME
    disp(ME);
    vol([1 : esize, sz(1) - esize : sz(1)], :, :) = 0;
end

end
