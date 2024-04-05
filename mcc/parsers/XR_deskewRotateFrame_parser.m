function [ds, dsr] = XR_deskewRotateFrame_parser(framePath, xyPixelSize, dz, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('framePath', @(x) ischar(x) || iscell(x)); 
ip.addRequired('xyPixelSize', @(x) isscalar(x) || ischar(x)); 
ip.addRequired('dz', @(x) isscalar(x) || ischar(x)); 
ip.addParameter('ObjectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ZstageScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Crop', false, @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('Reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Rotate', false, @(x) islogical(x) || ischar(x));
ip.addParameter('InputBbox', [], @(x) isnumeric(x) || ischar(x)); % bounding box apply to input
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('sCMOSCameraFlip', false, @(x) islogical(x) || ischar(x));
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
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})) || ischar(x));
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
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
xStepThresh = pr.xStepThresh;
aname = pr.aname;
ZoffsetCorrection = pr.ZoffsetCorrection;
DSRCombined = pr.DSRCombined;
resample = pr.resample;
Interp = pr.Interp;
surffix = pr.surffix;
uuid = pr.uuid;

if ischar(framePath) && strcmp(framePath(1), '{')
    framePath = eval(framePath);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(ObjectiveScan)
    ObjectiveScan = str2num(ObjectiveScan);
end
if ischar(ZstageScan)
    ZstageScan = str2num(ZstageScan);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(Crop)
    Crop = str2num(Crop);
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(Reverse)
    Reverse = str2num(Reverse);
end
if ischar(Rotate)
    Rotate = str2num(Rotate);
end
if ischar(InputBbox)
    InputBbox = str2num(InputBbox);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(sCMOSCameraFlip)
    sCMOSCameraFlip = str2num(sCMOSCameraFlip);
end
if ischar(LLFFCorrection)
    LLFFCorrection = str2num(LLFFCorrection);
end
if ischar(LowerLimit)
    LowerLimit = str2num(LowerLimit);
end
if ischar(constOffset)
    constOffset = str2num(constOffset);
end
if ischar(BKRemoval)
    BKRemoval = str2num(BKRemoval);
end
if ischar(Save16bit)
    Save16bit = str2num(Save16bit);
end
if ischar(RescaleRotate)
    RescaleRotate = str2num(RescaleRotate);
end
if ischar(save3DStack)
    save3DStack = str2num(save3DStack);
end
if ischar(SaveMIP)
    SaveMIP = str2num(SaveMIP);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(xStepThresh)
    xStepThresh = str2num(xStepThresh);
end
if ischar(ZoffsetCorrection)
    ZoffsetCorrection = str2num(ZoffsetCorrection);
end
if ischar(DSRCombined)
    DSRCombined = str2num(DSRCombined);
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
    DSRCombined=DSRCombined, resample=resample, Interp=Interp, surffix=surffix, ...
    uuid=uuid);

end

