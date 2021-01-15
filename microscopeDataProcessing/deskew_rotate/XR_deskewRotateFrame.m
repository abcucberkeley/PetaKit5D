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
% xruan (10/05/2020): convert to single rather than double to save space.
%                     Add image block based processing for tall images
% xruan (10/06/2020): add support of a list of files and combine them first
%                     before deskew. The files should be orderred from top to bottom
% xruan (10/11/2020): add the combined deskew, rotate and resampling function. 
% xruan (12/05/2020): add option to remove background 
%                     add option to flip z stack in raw data (for negative X interval)
% xruan (12/18/2020): add support to not save 3D stack


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
ip.addParameter('flipZstack', false, @islogical);
% sCMOS camera flip
ip.addParameter('sCMOSCameraFlip', false, @islogical);
% LLSM Flat-FieldCorrection
ip.addParameter('LLFFCorrection', false, @islogical);
ip.addParameter('LowerLimit', 0.4, @isnumeric); % this value is the lowest
ip.addParameter('LSImage', '' , @isstr);
ip.addParameter('BackgroundImage', '' , @isstr);
ip.addParameter('constOffset', [], @(x) isempty(x) || isnumeric(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('BKRemoval', false, @islogical);
ip.addParameter('Save16bit', false , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('RescaleRotate', false , @islogical); % Rescale rotated data to [0 65535]
ip.addParameter('save3DStack', true , @islogical); % option to save 3D stack or not
ip.addParameter('SaveMIP', true , @islogical); % save MIP-z for ds and dsr. 
ip.addParameter('aname', '', @isstr); % XR allow user-defined result path
ip.addParameter('ZoffsetCorrection', false, @islogical); % xruan: add option for correction of z offset
ip.addParameter('DSRCombined', true, @islogical); % combined processing 
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x)); % resampling after rotation 
ip.addParameter('saveZarr', false, @islogical); % save as zarr
ip.addParameter('blockSize', [500, 500, 500], @isnumeric); % save as zarr
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.addParameter('surffix', '', @isstr); % suffix for the folder
ip.addParameter('uuid', '', @isstr);

ip.parse(framePath, xyPixelSize, dz, varargin{:});

pr = ip.Results;
Crop = pr.Crop;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
ObjectiveScan = pr.ObjectiveScan;
flipZstack = pr.flipZstack;
LLFFCorrection = pr.LLFFCorrection;
BKRemoval = pr.BKRemoval;
LSImage = pr.LSImage;
BackgroundImage = pr.BackgroundImage;
save3DStack = pr.save3DStack;
SaveMIP = pr.SaveMIP;
DSRCombined = pr.DSRCombined;
resample = pr.resample;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
Interp = pr.Interp;
surffix = pr.surffix;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end

% not convert to negative it for the wrapper
% if Reverse
%     SkewAngle = -abs(SkewAngle);
% end

% decide zAniso
if ObjectiveScan
    zAniso = dz / xyPixelSize;
else
    theta = SkewAngle * pi / 180;
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
end

%% deskew frame
% check if needs to combine parts of volume
if ischar(framePath)
    framePath = {framePath};
end
combinedFrame = false;
if numel(framePath) > 1
    % use the first one as the name of the combined file.
    % assume the combined frame can be fitted to memory. (In future may
    % implement Image Block method)
    combinedFrame = true;
end

% first check if the file exists
for i = 1 : numel(framePath)
    if ~exist(framePath{i}, 'file')
        warning('%s does not exist!', framePath{i});
        return;
    end
end

% Create DS result dire
[rt, fsname] = fileparts(framePath{1});
if ~DSRCombined    
    dsPath = sprintf('%s/DS%s/', rt, surffix);
    mkdir(dsPath);
    fileattrib(dsPath, '+w', 'g');
    dsFullname = [dsPath, fsname, '.tif'];
end

if (~DSRCombined && (~exist(dsFullname, 'file') || ip.Results.Overwrite)) || DSRCombined
    % frame = double(readtiff(framePath));
    if combinedFrame
        frame_cell = cell(numel(framePath), 1);
        for i = 1 : numel(framePath)
            frame_cell{i} = readtiff(framePath{i});
        end
        frame = single(cat(3, frame_cell{:}));
        clear frame_cell;
    else
        frame = single(readtiff(framePath{1}));
    end
    
    if flipZstack
        frame = flip(frame, 3);
    end
        
    % flat field correction
    if LLFFCorrection
        LSIm = readtiff(LSImage);
        BKIm = readtiff(BackgroundImage);
        frame = GU_LSFlatFieldCorrection(frame,LSIm,BKIm,'LowerLimit', ip.Results.LowerLimit, ...
            'constOffset', ip.Results.constOffset);
    end
    % remove camera background
    if BKRemoval
        BKIm = readtiff(BackgroundImage);
        frame = XR_CameraBackgroundRemoval(frame, BKIm, 'constOffset', ip.Results.constOffset);
    end

    if ~DSRCombined
        fprintf('Deskew frame %s...\n', framePath{1});
        try 
            ds = deskewFrame3D(frame, SkewAngle, dz, xyPixelSize, Reverse, ...
                'Crop', Crop, 'Interp', Interp); 
            clear frame;
        catch ME
            disp(ME);
            sz = size(frame);
            if sz(3) > 200
                fprintf('Use image block method for deskew...\n');
                % only split in y-axis, it is not right when splitting from
                % other axes. 
                blockSize = sz;
                blockSize(1) = min(ceil(blockSize(1) / 8), 150);
                bim = blockedImage(frame, 'BlockSize', blockSize);
                OutputLocation = sprintf('%s/%s_%s', dsPath, fsname, uuid);
                BorderSize = [5, 0, 0];
                TrimBorder = true;

                bo = apply(bim, @(bs) deskewFrame3D(single(bs.Data), SkewAngle, dz, ...
                    xyPixelSize, Reverse, 'crop', Crop, 'Interp', Interp), 'blockSize', bim.BlockSize, ...
                    'OutputLocation', OutputLocation, 'BorderSize', BorderSize, 'TrimBorder', TrimBorder, ...
                    'useParallel', false);
                clear frame;
                ds = gather(bo);
                rmdir(OutputLocation, 's');
                clear bim bo;
            end
        end

        % save MIP
        if SaveMIP
            dsMIPPath = sprintf('%s/MIPs/', dsPath);
            if ~exist(dsMIPPath, 'dir')
                mkdir(dsMIPPath);
                fileattrib(dsMIPPath, '+w', 'g');            
            end

            dsMIPname = sprintf('%s%s_MIP_z.tif', dsMIPPath, fsname);
            if ip.Results.Save16bit
                writetiff(uint16(max(ds, [], 3)), dsMIPname);
            else
                writetiff(single(max(ds, [], 3)), dsMIPname);
            end            
        end

        dsTempname = sprintf('%s%s_%s.tif', dsPath, fsname, uuid);
        if save3DStack
            if ip.Results.Save16bit
                writetiff(uint16(ds), dsTempname);
            else
                writetiff(single(ds), dsTempname);
            end
            movefile(dsTempname, dsFullname);
        end
    end
end

%% rotate frame

if ip.Results.Rotate || DSRCombined
    dsrPath = sprintf('%s/DSR%s/', rt, surffix);
    mkdir(dsrPath);
    
    if saveZarr
        dsrFullname = [dsrPath, fsname, '.zarr'];        
    else
        dsrFullname = [dsrPath, fsname, '.tif'];
    end
    
    if (~saveZarr && ~exist(dsrFullname, 'file')) || (saveZarr && ~exist(dsrFullname, 'dir'))
        if ~DSRCombined
            fprintf('Rotate frame %s...\n', framePath{1});
            if ~exist('ds', 'var')
                ds = single(readtiff(dsFullname));
            end
            dsr = rotateFrame3D(ds, ip.Results.SkewAngle, zAniso, ip.Results.Reverse,...
                'Crop', true, 'ObjectiveScan', ObjectiveScan, 'Interp', Interp);
            clear ds;
        else
            fprintf('Deskew, Rotate and resample for frame %s...\n', framePath{1});            
            dsr = deskewRotateFrame3D(frame, ip.Results.SkewAngle, dz, xyPixelSize, ...
                'reverse', ip.Results.Reverse, 'Crop', true, 'ObjectiveScan', ObjectiveScan, ...
                'resample', resample, 'Interp', Interp);
        end
        
        % save MIP
        if SaveMIP
            dsrMIPPath = sprintf('%s/MIPs/', dsrPath);
            if ~exist(dsrMIPPath, 'dir')
                mkdir(dsrMIPPath);
                fileattrib(dsrMIPPath, '+w', 'g');
            end

            dsrMIPname = sprintf('%s%s_MIP_z.tif', dsrMIPPath, fsname);
            if ip.Results.Save16bit
                writetiff(uint16(max(dsr, [], 3)), dsrMIPname);
            else
                writetiff(single(max(dsr, [], 3)), dsrMIPname);
            end
        end
        
        if saveZarr
            dsrTempName = sprintf('%s%s_%s.zarr', dsrPath, fsname, uuid);
        else
            dsrTempName = sprintf('%s%s_%s.tif', dsrPath, fsname, uuid);            
        end
        
        if ip.Results.RescaleRotate
            iRange = [min(dsr(:)), max(dsr(:))];
            dsr = scaleContrast(dsr, iRange, [0 65535]);
        end
        
        if ip.Results.Save16bit
            dsr = uint16(dsr);
        else
            dsr = single(dsr);
        end
        
        if save3DStack
            if saveZarr
                bim = blockedImage(dsr);
                blockSize = min(size(dsr), blockSize);
                write(bim, dsrTempName, 'Adapter', ZarrAdapter, 'BlockSize', blockSize);
            else
                writetiff(dsr, dsrTempName);
            end
            movefile(dsrTempName, dsrFullname);
        end
    end
end


end
