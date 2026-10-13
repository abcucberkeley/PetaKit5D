function [vol] = texture(sz, sigma, intRange)
% generate a smooth random texture (single), for tests that need unique
% content everywhere, e.g., registration and stitching. Uses the global
% random stream, which pktest.TestCase seeds for every test.
%
%        sz : 1x2 or 1x3 vector. Image size.
%     sigma : a number (default: 2). Gaussian smoothing sigma in voxels.
%  intRange : 1x2 vector (default: [100, 2000]). Output intensity range.

if nargin < 2 || isempty(sigma)
    sigma = 2;
end
if nargin < 3 || isempty(intRange)
    intRange = [100, 2000];
end

vol = randn(sz, 'single');
if numel(sz) == 2 || sz(3) == 1
    vol = imgaussfilt(vol, sigma);
else
    vol = imgaussfilt3(vol, sigma);
end
vol = (vol - min(vol(:))) / (max(vol(:)) - min(vol(:)));
vol = single(intRange(1) + vol * (intRange(2) - intRange(1)));

end
