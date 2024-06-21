function [] = XR_RotateFrame3D(framePaths, xyPixelSize, dz, varargin)
% Rotate data for a single file or a batch of files with same parameters
% 
% Based on deskewData.m
%
% Author: Xiongtao Ruan (03/18/2020)
% 
% xruan (11/30/2021): add support for bounding box crop and resample

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('framePaths', @(x) ischar(x) || iscell(x)); 
ip.addRequired('xyPixelSize', @isscalar); 
ip.addRequired('dz', @isscalar); 
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('Crop', true, @islogical);
ip.addParameter('bbox', [], @(x) isnumeric(x));
ip.addParameter('resample', [], @(x) isnumeric(x)); % resampling after rotation 
ip.addParameter('SkewAngle', 31.5, @isscalar);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('sCMOSCameraFlip', false, @islogical);
ip.addParameter('save16bit', true , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('uuid', '', @ischar);

ip.parse(framePaths, xyPixelSize, dz, varargin{:});

warning('off', 'MATLAB:MKDIR:DirectoryExists');

if ischar(framePaths)
    framePaths = {framePaths};
end

pr = ip.Results;
Overwrite = pr.Overwrite;
Crop = pr.Crop;
bbox = pr.bbox;
resample = pr.resample;
Reverse = pr.Reverse;
SkewAngle = pr.SkewAngle;
objectiveScan = pr.objectiveScan;
save16bit = pr.save16bit;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
    
% decide zAniso
if objectiveScan
    zAniso = dz / xyPixelSize;
else
    theta = SkewAngle * pi / 180;
    dz0 = sin(theta) * dz; % ~0.25 for dz0 = 0.45
    zAniso = dz0 / xyPixelSize;
end

for f = 1 : numel(framePaths)
    framePath = framePaths{f};
    
    % first check if the file exists
    if ~exist(framePath, 'file')
        warning('%s does not exist!', framePath);
        continue;
    end

    % check if the result exists
    [rt, fsname] = fileparts(framePath);
    fname = [fsname, '.tif'];
    
    rtPath = [rt, '/', 'Rotated/'];
    mkdir(rtPath);
    rtFullname = [rtPath, fname];

    if ~exist(rtFullname, 'file') || Overwrite
        fprintf('Rotate frame %s...\n', framePath);
        if ~exist('ds', 'var')
            im = double(readtiff(framePath));
        end
        im_rt = rotateFrame3D(im, SkewAngle, zAniso, Reverse, 'Crop', Crop, ...
            'resample', resample, 'objectiveScan', objectiveScan);
        
        if ~isempty(bbox)
            im_rt = im_rt(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));            
        end
        
        rtTempName = sprintf('%s%s_%s.tif', rtPath, fsname, uuid);
        if save16bit
            writetiff(uint16(im_rt), rtTempName);
        else
            writetiff(single(im_rt), rtTempName);
        end

        movefile(rtTempName, rtFullname);
    end
end

end
