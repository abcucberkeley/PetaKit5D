% Overlays a color mask on a grayscale input image.
%
% Inputs:
%         img      : Grayscale input image.
%         masks    : Overlay masks (cell array). Nonzero elements are colored.
%         colors   : Overlay colors (cell array of 3-element RGB vectors).
%         {iRange} : dynamic range of 'img'.
%
% Copyright (C) 2017, Danuser Lab - UTSouthwestern 
%
% This file is part of CMEAnalysis_Package.
% 
% CMEAnalysis_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CMEAnalysis_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with CMEAnalysis_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

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
