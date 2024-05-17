function frame = processFFCorrectionFrame(frame, LSImage, BackgroundImage, constOffset, LowerLimit, removeFFImBackground, LSRescale, castDataType)
% flat field correction for a frame in memory

if nargin < 8
    castDataType =  true;
end
if nargin < 7
    LSRescale =  true;
end
if nargin < 6
    removeFFImBackground = true;
end

if nargin < 5
    LowerLimit = 0.4;
end

% ff correction
LSIm = readtiff(LSImage);
BKIm = readtiff(BackgroundImage);
frame = XR_LSFlatFieldCorrection(frame,LSIm,BKIm,'LowerLimit', LowerLimit, ...
    'constOffset', constOffset, 'removeFFImBackground', removeFFImBackground, ...
    'LSRescale', LSRescale, 'castDataType', castDataType);

end

