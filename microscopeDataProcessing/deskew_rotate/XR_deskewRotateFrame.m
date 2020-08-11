function [] = XR_deskewRotateFrame(framePath, xyPixelSize, dz, varargin)
% Deskew and/or rotate data for a single file
% 
%
% Author: Xiongtao Ruan (02/25/2020)
%
% Based on deskewData.m
% xruan (07/13/2020): set not rescale data to true for rotated data
% xruan (07/15/2020): add flatfield correction
% xruan (07/19/2020): add save of MIPs by default


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('framePath'); 
ip.addRequired('xyPixelSize'); 
ip.addRequired('dz'); 
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('MovieSelector', 'cell', @ischar);
ip.addParameter('Crop', false, @islogical);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('Rotate', false, @islogical);
ip.addParameter('CheckFrameMismatch', false, @islogical);
ip.addParameter('LoadSettings', false, @islogical);
% sCMOS camera flip
ip.addParameter('sCMOSCameraFlip', false, @islogical);
% LLSM Flat-FieldCorrection
ip.addParameter('LLFFCorrection', false, @islogical);
ip.addParameter('LowerLimit', 0.4, @isnumeric); % this value is the lowest
ip.addParameter('LSImage', '' , @isstr);
ip.addParameter('BackgroundImage', '' , @isstr);
ip.addParameter('constOffset', [], @(x) isempty(x) || isnumeric(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('Save16bit', false , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('RescaleRotate', false , @islogical); % Rescale rotated data to [0 65535]
ip.addParameter('SaveMIP', true , @islogical); % save MIP-z for ds and dsr. 
ip.addParameter('aname', '', @isstr); % XR allow user-defined result path
ip.addParameter('ZoffsetCorrection', false, @islogical); % xruan: add option for correction of z offset
ip.addParameter('uuid', '', @isstr);

ip.parse(framePath, xyPixelSize, dz, varargin{:});

pr = ip.Results;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
ObjectiveScan = pr.ObjectiveScan;
LLFFCorrection = pr.LLFFCorrection;
LSImage = pr.LSImage;
BackgroundImage = pr.BackgroundImage;
SaveMIP = pr.SaveMIP;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end

if Reverse
    SkewAngle = -abs(SkewAngle);
end

% decide zAniso
if ObjectiveScan
    zAniso = dz / xyPixelSize;
else
    theta = SkewAngle * pi / 180;
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
end

%% deskew frame

fprintf('Deskew frame %s...\n', framePath);

% first check if the file exists
if ~exist(framePath, 'file')
    warning('%s does not exist!', framePath);
    return;
end

% check if the result exists
[rt, fname] = fileparts(framePath);
dsPath = sprintf('%s/DS/', rt);
mkdir(dsPath);

dsFullname = [dsPath, fname, '.tif'];

if ~exist(dsFullname, 'file') || ip.Results.Overwrite
    frame = double(readtiff(framePath));
    % flat field correction
    if LLFFCorrection
        LSIm = readtiff(LSImage);
        BKIm = readtiff(BackgroundImage);
        frame = GU_LSFlatFieldCorrection(frame,LSIm,BKIm,'LowerLimit', ip.Results.LowerLimit, ...
            'constOffset', ip.Results.constOffset);
    end

    ds = deskewFrame3D(frame, ip.Results.SkewAngle, dz, xyPixelSize, ip.Results.Reverse, ...
        'Crop', ip.Results.Crop); 
    
    % save MIP
    if SaveMIP
        dsMIPPath = sprintf('%s/MIPs/', dsPath);
        if ~exist(dsMIPPath, 'dir')
            mkdir(dsMIPPath);
            fileattrib(dsMIPPath, '+w', 'g');            
        end
        
        dsMIPname = sprintf('%s%s_MIP_z.tif', dsMIPPath, fname);
        writetiff(uint16(max(ds, [], 3)), dsMIPname);
    end
    
    dsTempname = sprintf('%s%s_%s.tif', dsPath, fname, uuid);
    if ip.Results.Save16bit
        writetiff(uint16(ds), dsTempname);
    else
        writetiff(single(ds), dsTempname);
    end
    movefile(dsTempname, dsFullname);
end

%% rotate frame

if ip.Results.Rotate 
    dsrPath = sprintf('%s/DSR/', rt);
    mkdir(dsrPath);

    dsrFullname = [dsrPath, fname, '.tif'];
    
    if ~exist(dsrFullname, 'file')
        fprintf('Rotate frame %s...\n', framePath);
        if ~exist('ds', 'var')
            ds = double(readtiff(dsFullname));
        end
        dsr = rotateFrame3D(ds, ip.Results.SkewAngle, zAniso, ip.Results.Reverse,...
            'Crop', true, 'ObjectiveScan', ObjectiveScan);
        
        % save MIP
        if SaveMIP
            dsrMIPPath = sprintf('%s/MIPs/', dsrPath);
            if ~exist(dsrMIPPath, 'dir')
                mkdir(dsrMIPPath);
                fileattrib(dsrMIPPath, '+w', 'g');
            end

            dsrMIPname = sprintf('%s%s_MIP_z.tif', dsrMIPPath, fname);
            writetiff(uint16(max(dsr, [], 3)), dsrMIPname);
        end

        dsrTempName = sprintf('%s%s_%s.tif', dsrPath, fname, uuid);
%         if ip.Results.Save16bit
%             if ip.Results.RescaleRotate
%                 iRange = [min(ds(:)), max(ds(:))];
%                 writetiff(uint16(scaleContrast(dsr, iRange, [0 65535])), dsrTempName);
%             else
%                 writetiff(uint16(dsr), dsrTempName);
%             end
%         else
%             writetiff(single(dsr), dsrTempName);
%         end

        if ip.Results.RescaleRotate
            iRange = [min(ds(:)), max(ds(:))];
            writetiff(uint16(scaleContrast(dsr, iRange, [0 65535])), dsrTempName);
        else
            writetiff(uint16(dsr), dsrTempName);
        end

        movefile(dsrTempName, dsrFullname);
    end
end


end
