function [ds, dsr] = XR_deskewRotateFrame_parser(framePath, xyPixelSize, dz, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('framePath', @(x) ischar(x) || iscell(x)); 
ip.addRequired('xyPixelSize', @(x) isscalar(x) || ischar(x)); 
ip.addRequired('dz', @(x) isscalar(x) || ischar(x)); 
ip.addParameter('DSDirName', 'DS/', @ischar);
ip.addParameter('DSRDirName', 'DSR/', @ischar);
ip.addParameter('objectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('zStageScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('crop', false, @(x) islogical(x) || ischar(x));
ip.addParameter('skewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addParameter('reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('rotate', false, @(x) islogical(x) || ischar(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x) || ischar(x)); % bounding box apply to input
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('sCMOSCameraFlip', false, @(x) islogical(x) || ischar(x));
ip.addParameter('FFCorrection', false, @(x) islogical(x) || ischar(x));
ip.addParameter('lowerLimit', 0.4, @(x) isnumeric(x) || ischar(x)); % this value is the lowest
ip.addParameter('FFImage', '' , @ischar);
ip.addParameter('backgroundImage', '' , @ischar);
ip.addParameter('constOffset', [], @(x) isnumeric(x) || ischar(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('BKRemoval', false, @(x) islogical(x) || ischar(x));
ip.addParameter('save16bit', true , @(x) islogical(x) || ischar(x)); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('rescaleRotate', false , @(x) islogical(x) || ischar(x)); % Rescale rotated data to [0 65535]
ip.addParameter('save3DStack', true , @(x) islogical(x) || ischar(x)); % option to save 3D stack or not
ip.addParameter('saveMIP', true , @(x) islogical(x) || ischar(x)); % save MIP-z for ds and dsr. 
ip.addParameter('saveZarr', false , @(x) islogical(x) || ischar(x)); % save as zarr
ip.addParameter('blockSize', [500, 500, 500] , @(x) isnumeric(x) || ischar(x)); % save as zarr
ip.addParameter('xStepThresh', 2.0, @(x) isnumeric(x) || ischar(x)); % 2.344 for ds=0.3, 2.735 for ds=0.35
ip.addParameter('zOffsetCorrection', false, @(x) islogical(x) || ischar(x)); % xruan: add option for correction of z offset
ip.addParameter('DSRCombined', true, @(x) islogical(x) || ischar(x)); % combined processing 
ip.addParameter('resampleFactor', [], @(x) isnumeric(x) || ischar(x)); % resampling after rotation 
ip.addParameter('interpMethod', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})) || ischar(x));
ip.addParameter('surffix', '', @ischar); % suffix for the folder
ip.addParameter('uuid', '', @ischar);

ip.parse(framePath, xyPixelSize, dz, varargin{:});

pr = ip.Results;
DSDirName = pr.DSDirName;
DSRDirName = pr.DSRDirName;
objectiveScan = pr.objectiveScan;
zStageScan = pr.zStageScan;
overwrite = pr.overwrite;
crop = pr.crop;
skewAngle = pr.skewAngle;
reverse = pr.reverse;
rotate = pr.rotate;
inputBbox = pr.inputBbox;
flipZstack = pr.flipZstack;
sCMOSCameraFlip = pr.sCMOSCameraFlip;
FFCorrection = pr.FFCorrection;
lowerLimit = pr.lowerLimit;
FFImage = pr.FFImage;
backgroundImage = pr.backgroundImage;
constOffset = pr.constOffset;
BKRemoval = pr.BKRemoval;
save16bit = pr.save16bit;
rescaleRotate = pr.rescaleRotate;
save3DStack = pr.save3DStack;
saveMIP = pr.saveMIP;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
xStepThresh = pr.xStepThresh;
zOffsetCorrection = pr.zOffsetCorrection;
DSRCombined = pr.DSRCombined;
resampleFactor = pr.resampleFactor;
interpMethod = pr.interpMethod;
surffix = pr.surffix;
uuid = pr.uuid;

if ischar(framePath) && ~isempty(framePath) && strcmp(framePath(1), '{')
    framePath = eval(framePath);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(objectiveScan)
    objectiveScan = str2num(objectiveScan);
end
if ischar(zStageScan)
    zStageScan = str2num(zStageScan);
end
if ischar(overwrite)
    overwrite = str2num(overwrite);
end
if ischar(crop)
    crop = str2num(crop);
end
if ischar(skewAngle)
    skewAngle = str2num(skewAngle);
end
if ischar(reverse)
    reverse = str2num(reverse);
end
if ischar(rotate)
    rotate = str2num(rotate);
end
if ischar(inputBbox)
    inputBbox = str2num(inputBbox);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(sCMOSCameraFlip)
    sCMOSCameraFlip = str2num(sCMOSCameraFlip);
end
if ischar(FFCorrection)
    FFCorrection = str2num(FFCorrection);
end
if ischar(lowerLimit)
    lowerLimit = str2num(lowerLimit);
end
if ischar(constOffset)
    constOffset = str2num(constOffset);
end
if ischar(BKRemoval)
    BKRemoval = str2num(BKRemoval);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(rescaleRotate)
    rescaleRotate = str2num(rescaleRotate);
end
if ischar(save3DStack)
    save3DStack = str2num(save3DStack);
end
if ischar(saveMIP)
    saveMIP = str2num(saveMIP);
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
if ischar(zOffsetCorrection)
    zOffsetCorrection = str2num(zOffsetCorrection);
end
if ischar(DSRCombined)
    DSRCombined = str2num(DSRCombined);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end

XR_deskewRotateFrame(framePath, xyPixelSize, dz, DSDirName=DSDirName, DSRDirName=DSRDirName, ...
    objectiveScan=objectiveScan, zStageScan=zStageScan, overwrite=overwrite, ...
    crop=crop, skewAngle=skewAngle, reverse=reverse, rotate=rotate, inputBbox=inputBbox, ...
    flipZstack=flipZstack, sCMOSCameraFlip=sCMOSCameraFlip, FFCorrection=FFCorrection, ...
    lowerLimit=lowerLimit, FFImage=FFImage, backgroundImage=backgroundImage, ...
    constOffset=constOffset, BKRemoval=BKRemoval, save16bit=save16bit, rescaleRotate=rescaleRotate, ...
    save3DStack=save3DStack, saveMIP=saveMIP, saveZarr=saveZarr, blockSize=blockSize, ...
    xStepThresh=xStepThresh, zOffsetCorrection=zOffsetCorrection, DSRCombined=DSRCombined, ...
    resampleFactor=resampleFactor, interpMethod=interpMethod, surffix=surffix, ...
    uuid=uuid);

end

