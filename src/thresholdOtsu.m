function [level, varargout] = thresholdOtsu(imageIn,varargin)
% Select a thresholding level using Otsu's method
%
% level = thresholdOtsu(imageIn)
% level = thresholdOtsu(imageIn,showPlots)
% 
% This function selects a threshold for the input fluorescence image by
% analyzing the image's intensity distribution. This requires good signal-
% to-noise, and a significant amount of background in the image. The
% threshold is selected using Otsu's method
% 
% Input:
% 
%   imageIn - The N-Dimensional image to be thresholded.
% 
% 
%   showPlots - If true, a plot of the histogram and an overlay of the mask
%   on the image will be shown. The overlay plot only works if the image is
%   2D.
% 
% 
% Output:
% 
% 
%   level - The intensity value selected for thresholding.
%   mask (optional) -- the thresholded foreground mask is returned as 
%                      an optional output argument
%                      (added by Deepak on Mar, 2013) 
%
% Revamped from OtsuSeg
%
% Sebastien Besson, 5/2011

ip=inputParser;
ip.addRequired('imageIn',@isnumeric);
ip.addOptional('showPlots',0,@isnumeric)
ip.parse(imageIn,varargin{:});
showPlots=ip.Results.showPlots;

%Convert to double if necessary
imageIn = double(imageIn);

%find nonzero values (due to masking)
nzInd = find(imageIn);

%get minumum and maximum pixel values in image
minSignal = min(imageIn(nzInd));
maxSignal = max(imageIn(nzInd));

%normalize nonzero value between 0 and 1
imageInNorm = zeros(size(imageIn));
imageInNorm(nzInd) = (imageIn(nzInd)- minSignal) / (maxSignal - minSignal);

level = graythresh(imageInNorm);

level = level*(maxSignal - minSignal)+minSignal;

if showPlots 
    imageMask = imageIn >= level;
    figure;
    imagesc(imageIn);
    hold on
    contour(imageMask,'w')
    colormap hot
end

if nargout > 1
    varargout{1} = double(imageIn >= level);
end