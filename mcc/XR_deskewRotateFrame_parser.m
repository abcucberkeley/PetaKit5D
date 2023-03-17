function [] = XR_deskewRotateFrame_parser(framePath, xyPixelSize, dz, varargin)

%#function XR_deskewRotateFrame

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('framePath'); 
ip.addRequired('xyPixelSize'); 
ip.addRequired('dz'); 
ip.addParameter('ObjectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ZstageScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Crop', false, @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('Reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Rotate', false, @(x) islogical(x) || ischar(x));
ip.addParameter('InputBbox', [], @(x) isnumeric(x) || ischar(x)); % bounding box apply to input
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
% sCMOS camera flip
ip.addParameter('sCMOSCameraFlip', false, @(x) islogical(x) || ischar(x));
% LLSM Flat-FieldCorrection
ip.addParameter('LLFFCorrection', false, @(x) islogical(x) || ischar(x));
ip.addParameter('LowerLimit', 0.4, @(x) isnumeric(x) || ischar(x)); % this value is the lowest
ip.addParameter('LSImage', '' , @ischar);
ip.addParameter('BackgroundImage', '' , @ischar);
ip.addParameter('constOffset', [], @(x) isnumeric(x) || ischar(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('BKRemoval', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Save16bit', false , @(x) islogical(x) || ischar(x)); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('RescaleRotate', false , @(x) islogical(x) || ischar(x)); % Rescale rotated data to [0 65535]
ip.addParameter('save3DStack', true , @(x) islogical(x) || ischar(x)); % option to save 3D stack or not
ip.addParameter('SaveMIP', true , @(x) islogical(x) || ischar(x)); % save MIP-z for ds and dsr. 
ip.addParameter('saveZarr', false , @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('blockSize', [500, 500, 500] , @(x) isnumeric(x) || ischar(x)); % save as zarr
ip.addParameter('xStepThresh', 2.0, @(x) isnumeric(x) || ischar(x)); % 2.344 for ds=0.3, 2.735 for ds=0.35
ip.addParameter('aname', '', @ischar); % XR allow user-defined result path
ip.addParameter('ZoffsetCorrection', false, @(x) islogical(x) || ischar(x)); % xruan: add option for correction of z offset
ip.addParameter('DSRCombined', true, @(x) islogical(x) || ischar(x)); % combined processing 
ip.addParameter('resample', [], @(x) isnumeric(x) || ischar(x)); % resampling after rotation 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})) && ischar(x));
ip.addParameter('surffix', '', @ischar); % suffix for the folder
ip.addParameter('uuid', '', @ischar);

ip.parse(framePath, xyPixelSize, dz, varargin{:});

pr = ip.Results;
ObjectiveScan = pr.ObjectiveScan;
ZstageScan = pr.ZstageScan;
Overwrite = pr.Overwrite;
Crop = pr.Crop;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
Rotate = pr.Rotate;
InputBbox = pr.InputBbox;
flipZstack = pr.flipZstack;
sCMOSCameraFlip = pr.sCMOSCameraFlip;
LLFFCorrection = pr.LLFFCorrection;
LowerLimit = pr.LowerLimit;
LSImage = pr.LSImage;
BackgroundImage = pr.BackgroundImage;
constOffset = pr.constOffset;
BKRemoval = pr.BKRemoval;
Save16bit = pr.Save16bit;
RescaleRotate = pr.RescaleRotate;
save3DStack = pr.save3DStack;
SaveMIP = pr.SaveMIP;
ZoffsetCorrection = pr.ZoffsetCorrection;
DSRCombined = pr.DSRCombined;
resample = pr.resample;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
xStepThresh = pr.xStepThresh;
aname = pr.aname;
Interp = pr.Interp;
surffix = pr.surffix;
uuid = pr.uuid;

if ischar(framePath) && strcmp(framePath(1), '{')
    framePath = eval(framePath);
end
if ischar(xyPixelSize)
    xyPixelSize = str2double(xyPixelSize);
end
if ischar(dz)
    dz = str2double(dz);
end
if ischar(ObjectiveScan)
    ObjectiveScan = strcmp(ObjectiveScan, 'true');
end
if ischar(ZstageScan)
    ZstageScan = strcmp(ZstageScan, 'true');
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
end
if ischar(Crop)
    Crop = strcmp(Crop, 'true');
end
if ischar(SkewAngle)
    SkewAngle = str2double(SkewAngle);
end
if ischar(Reverse)
    Reverse = strcmp(Reverse, 'true');
end
if ischar(Rotate)
    Rotate = strcmp(Rotate, 'true');
end
if ischar(InputBbox)
    InputBbox = str2num(InputBbox);
end
if ischar(flipZstack)
    flipZstack = strcmp(flipZstack, 'true');
end
if ischar(sCMOSCameraFlip)
    sCMOSCameraFlip = strcmp(sCMOSCameraFlip, 'true');
end
if ischar(LLFFCorrection)
    LLFFCorrection = strcmp(LLFFCorrection, 'true');
end
if ischar(LowerLimit)
    LowerLimit = str2double(LowerLimit);
end
if ischar(constOffset)
    constOffset = str2num(constOffset);
end
if ischar(BKRemoval)
    BKRemoval = strcmp(BKRemoval, 'true');
end
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit, 'true');
end
if ischar(RescaleRotate)
    RescaleRotate = strcmp(RescaleRotate, 'true');
end
if ischar(save3DStack)
    save3DStack = strcmp(save3DStack, 'true');
end
if ischar(SaveMIP)
    SaveMIP = strcmp(SaveMIP, 'true');
end
if ischar(saveZarr)
    saveZarr = strcmp(saveZarr, 'true');
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(xStepThresh)
    xStepThresh = str2double(xStepThresh);
end
if ischar(ZoffsetCorrection)
    ZoffsetCorrection = strcmp(ZoffsetCorrection, 'true');
end
if ischar(DSRCombined)
    DSRCombined = strcmp(DSRCombined, 'true');
end
if ischar(resample)
    resample = str2num(resample);
end

XR_deskewRotateFrame(framePath, xyPixelSize, dz, ObjectiveScan=ObjectiveScan, ...
    ZstageScan=ZstageScan, Overwrite=Overwrite, Crop=Crop, SkewAngle=SkewAngle, ...
    Reverse=Reverse, Rotate=Rotate, InputBbox=InputBbox, flipZstack=flipZstack, ...
    sCMOSCameraFlip=sCMOSCameraFlip, LLFFCorrection=LLFFCorrection, LowerLimit=LowerLimit, ...
    LSImage=LSImage, BackgroundImage=BackgroundImage, constOffset=constOffset, ...
    BKRemoval=BKRemoval, Save16bit=Save16bit, RescaleRotate=RescaleRotate, ...
    save3DStack=save3DStack, SaveMIP=SaveMIP, saveZarr=saveZarr, blockSize=blockSize, ...
    xStepThresh=xStepThresh, aname=aname, ZoffsetCorrection=ZoffsetCorrection, ...
    DSRCombined=DSRCombined, resample=resample, Interp=Interp, surffix=surffix, uuid=uuid);

end

