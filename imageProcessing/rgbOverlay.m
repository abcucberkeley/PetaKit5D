% Overlays a color mask on a grayscale input image.
%
% Inputs:
%         img      : Grayscale input image.
%         masks    : Overlay masks (cell array). Nonzero elements are colored.
%         colors   : Overlay colors (cell array of 3-element RGB vectors).
%         {iRange} : dynamic range of 'img'.

% Francois Aguet, April 2011.

function imgRGB = rgbOverlay(img, masks, colors, iRange)

if nargin<4
    iRange = [];
end

if isnumeric(masks)
    masks = {masks};
end
if isnumeric(colors)
    colors = {colors};
end

nm = length(masks);
img = scaleContrast(img, iRange);

combIdx = sum(cat(3, masks{:}),3) ~= 0;
imgRGB = img;
imgRGB(combIdx) = 0;

imgRGB = repmat(imgRGB, [1 1 3]);
for c = 1:3
    cMask = arrayfun(@(k) colors{k}(c)*masks{k}, 1:nm, 'UniformOutput', false);
    cMask = min(sum(cat(3, cMask{:}),3), 1);
    idx = find(cMask==1);
    tmp = imgRGB(:,:,c);
    tmp(idx) = img(idx);
    imgRGB(:,:,c) = tmp;
end
imgRGB = uint8(imgRGB);
