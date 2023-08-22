function frame = processFFCorrectionFrame(frame, LSImage, BackgroundImage, constOffset, LowerLimit, LSBackground)
% flat field correction for a frame in memory

if nargin < 6
    LSBackground = true;
end

if nargin < 5
    LowerLimit = 0.4;
end
% constOffset = [];

frame = single(frame);

% if flipZstack < 0
%     frame = flip(frame, 3);
% end
% 
% ff correction
LSIm = readtiff(LSImage);
BKIm = readtiff(BackgroundImage);
frame = XR_LSFlatFieldCorrection(frame,LSIm,BKIm,'LowerLimit', LowerLimit, ...
    'constOffset', constOffset, 'LSBackground', LSBackground);

end

