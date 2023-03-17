function [] = XR_RotateFrame3D_parser(framePaths, xyPixelSize, dz, varargin)


%#function XR_RotateFrame3D

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('framePaths'); 
ip.addRequired('xyPixelSize'); 
ip.addRequired('dz'); 
ip.addParameter('ObjectiveScan', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Crop', true, @(x) islogical(x) || ischar(x));
ip.addParameter('bbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('resample', [], @(x) isnumeric(x) || ischar(x)); % resampling after rotation 
ip.addParameter('SkewAngle', 31.5, @(x) isscalar(x) || ischar(x));
ip.addParameter('Reverse', false, @(x) islogical(x) || ischar(x));
ip.addParameter('sCMOSCameraFlip', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Save16bit', false , @(x) islogical(x) || ischar(x)); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('uuid', '', @ischar);

ip.parse(framePaths, xyPixelSize, dz, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
Crop = pr.Crop;
bbox = pr.bbox;
resample = pr.resample;
Reverse = pr.Reverse;
SkewAngle = pr.SkewAngle;
ObjectiveScan = pr.ObjectiveScan;
sCMOSCameraFlip = pr.sCMOSCameraFlip;
Save16bit = pr.Save16bit;
uuid = ip.Results.uuid;

if ischar(framePaths) && strcmp(framePaths(1), '{')
    framePaths = eval(framePaths);
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
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
end
if ischar(Crop)
    Crop = strcmp(Crop, 'true');
end
if ischar(bbox)
    bbox = str2num(bbox);
end
if ischar(SkewAngle)
    SkewAngle = str2double(SkewAngle);
end
if ischar(Reverse)
    Reverse = strcmp(Reverse, 'true');
end
if ischar(sCMOSCameraFlip)
    sCMOSCameraFlip = strcmp(sCMOSCameraFlip, 'true');
end
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit, 'true');
end

XR_RotateFrame3D(framePaths, xyPixelSize, dz, ObjectiveScan=ObjectiveScan, ...
    Overwrite=Overwrite, Crop=Crop, bbox=bbox, resample=resample, SkewAngle=SkewAngle, ...
    Reverse=Reverse, sCMOSCameraFlip=sCMOSCameraFlip, Save16bit=Save16bit, uuid=uuid);

end

