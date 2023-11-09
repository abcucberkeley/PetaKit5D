function frame = processFFCorrectionFrame(frame, LSImage, BackgroundImage, constOffset, LowerLimit, LSBackground, LSRescale, castDataType)
% flat field correction for a frame in memory

if nargin < 8
    castDataType =  false;
end
if nargin < 7
    LSRescale =  true;
end
if nargin < 6
    LSBackground = true;
end

if nargin < 5
    LowerLimit = 0.4;
end
% constOffset = [];

dtype = class(frame);
frame = single(frame);

% if flipZstack < 0
%     frame = flip(frame, 3);
% end
% 
% ff correction
LSIm = readtiff(LSImage);
BKIm = readtiff(BackgroundImage);
frame = XR_LSFlatFieldCorrection(frame,LSIm,BKIm,'LowerLimit', LowerLimit, ...
    'constOffset', constOffset, 'LSBackground', LSBackground, 'LSRescale', LSRescale);

if castDataType
    frame = cast(frame, dtype);
end

end

