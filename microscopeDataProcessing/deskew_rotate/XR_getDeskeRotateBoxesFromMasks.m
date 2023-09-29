function [inBbox, outBbox, dsr] = XR_getDeskeRotateBoxesFromMasks(maskFns, BorderSize, xyPixelSize, dz, SkewAngle, Reverse, resample, inputBbox)
% estimate input and output bbox for deskew and rotation by providing the
% MIP masks. 
%
% Author: Xiongtao Ruan (08/02/2023)


if nargin < 7
    resample = [];
end

if nargin < 8
    inputBbox = [];
end

bboxes = zeros(3, 4);
sz = zeros(1, 3);

for i = 1 : 3
    fn = maskFns{i};
    
    [~, ~, ext] = fileparts(fn);
    switch ext
        case {'.tif', '.tiff'}
            im = readtiff(fn);
        case '.zarr'
            im = readzarr(fn);
    end
    bbox = getImageBoundingBox(im);

    bboxes(i, :) = bbox;
    if i == 1
        sz(1 : 2) = size(im, [1, 2]);
    elseif i == 2
        sz(3) = size(im, 2);
        im_y = im;
    end
end

ys = min(bboxes(1, 1), bboxes(3, 1));
xs = min(bboxes(1, 2), bboxes(2, 1));
zs = min(bboxes(2, 2), bboxes(3, 2));

yt = max(bboxes(1, 3), bboxes(3, 3));
xt = max(bboxes(1, 4), bboxes(2, 3));
zt = max(bboxes(2, 4), bboxes(3, 4));

bbox = [ys, xs, zs, yt, xt, zt];

% load MIP mask y, and decide output size
if isempty(inputBbox)
    inputBbox = [1, 1, 1, sz];
end
inBbox = [max(inputBbox(1 : 3), bbox(1 : 3) - BorderSize), min(inputBbox(4 : 6), bbox(4 : 6) + BorderSize)];

% deskew and rotate MIP y
im_y = im_y(inBbox(2) : inBbox(5), inBbox(3) : inBbox(6));
frame = single(permute(im_y, [3, 1, 2]));
ObjectiveScan = false;

rs = [1, 1, 1];
if ~isempty(resample)
    rs(1) = resample(2);
    rs(3) = resample(3);
end
Interp = 'linear';
xStepThresh = 2;

dsr = deskewRotateFrame3D(frame, SkewAngle, dz, xyPixelSize, 'reverse', Reverse, ...
    'bbox', [], 'ObjectiveScan', ObjectiveScan, 'resample', rs, 'Interp', Interp, ...
    'xStepThresh', xStepThresh);

dsr = squeeze(dsr);

bbox = getImageBoundingBox(dsr);
outBbox = [1, max(1, bbox(1 : 2) - BorderSize(2 : 3)), round((inBbox(4) - inBbox(1)  + 1) ./ resample(1)), ...
    min(size(dsr), bbox(3 : 4) + BorderSize(2 : 3))];

end

