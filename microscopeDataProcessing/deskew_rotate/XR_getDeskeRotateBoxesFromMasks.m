function [inBbox, outBbox, dsr] = XR_getDeskeRotateBoxesFromMasks(maskFullpaths, BorderSize, xyPixelSize, dz, SkewAngle, Reverse, resample, inputBbox, flipZstack)
% estimate input and output bbox for deskew and rotation by providing the
% MIP masks. 
%
% resample: 1, 2 or 3 values in xyz order, as deskewRotateFrame3D takes them
% flipZstack: the data is flipped along z after the crop by inBbox
%
% Author: Xiongtao Ruan (08/02/2023)


if nargin < 7
    resample = [];
end

if nargin < 8
    inputBbox = [];
end

if nargin < 9
    flipZstack = false;
end

bboxes = zeros(3, 4);
sz = zeros(1, 3);

for i = 1 : 3
    fn = maskFullpaths{i};
    
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
if flipZstack
    im_y = flip(im_y, 2);
end
frame = uint16(permute(im_y, [3, 1, 2]));
objectiveScan = false;

% the MIP along y has one row: resampled along x and z here, the output rows
% are the input rows resampled along y
rs = [1, 1, 1];
rsy = 1;
if ~isempty(resample)
    rs_in = resample(:)';
    rs_in = [ones(1, 4 - numel(rs_in)) * rs_in(1), rs_in(2 : end)];
    rs(1) = rs_in(1);
    rs(3) = rs_in(3);
    rsy = rs_in(2);
end
interpMethod = 'linear';
xStepThresh = 2;
save16bit = all(rs == 1);

dsr = deskewRotateFrame3D(frame, SkewAngle, dz, xyPixelSize, 'reverse', Reverse, ...
    'bbox', [], 'objectiveScan', objectiveScan, 'resample', rs, 'interpMethod', interpMethod, ...
    'xStepThresh', xStepThresh, 'save16bit', save16bit);

dsr = squeeze(dsr);

bbox = getImageBoundingBox(dsr);
outBbox = [1, max(1, bbox(1 : 2) - BorderSize(2 : 3)), round((inBbox(4) - inBbox(1)  + 1) ./ rsy), ...
    min(size(dsr), bbox(3 : 4) + BorderSize(2 : 3))];

end

