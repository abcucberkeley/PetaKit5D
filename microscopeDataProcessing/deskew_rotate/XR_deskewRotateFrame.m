function [ds, dsr] = XR_deskewRotateFrame(framePath, xyPixelSize, dz, varargin)
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
% xruan (07/27/2021): add support for z-stage scan
% xruan (08/17/2021): add support to not load the full image of DS to the memory if it is too large. 
% xruan (10/21/2021): add support for zarr file as input
% xruan (01/25/2022): add support for bbox crop before processing
% xruan (07/07/2022): add support to resample skewed data for combined DSR
%   for big image (in case of OOM issue).
% xruan (12/06/2022): add predefined outSize for rotate in separate deskew/rotate

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('framePath', @(x) ischar(x) || iscell(x)); 
ip.addRequired('xyPixelSize', @isscalar); 
ip.addRequired('dz', @isscalar); 
ip.addParameter('DSDirName', 'DS/', @ischar);
ip.addParameter('DSRDirName', 'DSR/', @ischar);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('zStageScan', false, @islogical);
ip.addParameter('overwrite', false, @islogical);
ip.addParameter('crop', false, @islogical);
ip.addParameter('skewAngle', 32.45, @isscalar);
ip.addParameter('reverse', false, @islogical);
ip.addParameter('rotate', false, @islogical);
ip.addParameter('inputBbox', [], @isnumeric); % bounding box apply to input
ip.addParameter('flipZstack', false, @islogical);
% sCMOS camera flip
ip.addParameter('sCMOSCameraFlip', false, @islogical);
% LLSM Flat-FieldCorrection
ip.addParameter('FFCorrection', false, @islogical);
ip.addParameter('lowerLimit', 0.4, @isnumeric); % this value is the lowest
ip.addParameter('FFImage', '' , @ischar);
ip.addParameter('backgroundImage', '' , @ischar);
ip.addParameter('constOffset', [], @(x) isnumeric(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('BKRemoval', false, @islogical);
ip.addParameter('save16bit', true , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('rescaleRotate', false , @islogical); % Rescale rotated data to [0 65535]
ip.addParameter('save3DStack', true , @islogical); % option to save 3D stack or not
ip.addParameter('saveMIP', true , @islogical); % save MIP-z for ds and dsr. 
ip.addParameter('saveZarr', false , @islogical); % save as zarr
ip.addParameter('blockSize', [500, 500, 500] , @isnumeric); % save as zarr
ip.addParameter('xStepThresh', 2.0, @isnumeric); % 2.344 for ds=0.3, 2.735 for ds=0.35
ip.addParameter('zOffsetCorrection', false, @islogical); % xruan: add option for correction of z offset
ip.addParameter('DSRCombined', true, @islogical); % combined processing 
ip.addParameter('resampleFactor', [], @(x) isnumeric(x)); % resampling after rotation 
ip.addParameter('interpMethod', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.addParameter('surffix', '', @ischar); % suffix for the folder
ip.addParameter('uuid', '', @ischar);

ip.parse(framePath, xyPixelSize, dz, varargin{:});

pr = ip.Results;
DSDirName = pr.DSDirName;
DSRDirName = pr.DSRDirName;
crop = pr.crop;
skewAngle = pr.skewAngle;
reverse = pr.reverse;
objectiveScan = pr.objectiveScan;
zStageScan = pr.zStageScan;
inputBbox = pr.inputBbox;
flipZstack = pr.flipZstack;
FFCorrection = pr.FFCorrection;
lowerLimit = pr.lowerLimit;
BKRemoval = pr.BKRemoval;
FFImage = pr.FFImage;
backgroundImage = pr.backgroundImage;
constOffset = pr.constOffset;
save16bit = pr.save16bit;
save3DStack = pr.save3DStack;
saveMIP = pr.saveMIP;
DSRCombined = pr.DSRCombined;
resampleFactor = pr.resampleFactor;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
xStepThresh = pr.xStepThresh;
interpMethod = pr.interpMethod;
surffix = pr.surffix;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end

% not convert to negative it for the wrapper
% if reverse
%     skewAngle = -abs(skewAngle);
% end

% decide zAniso
if objectiveScan
    zAniso = dz / xyPixelSize;
elseif zStageScan
    theta = skewAngle * pi / 180;
    zAniso = cos(abs(theta)) * dz / xyPixelSize;    
else
    theta = skewAngle * pi / 180;
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
end

if zStageScan
    skewAngle_1 = 90 - skewAngle;
    reverse =  ~reverse;
else
    skewAngle_1 = skewAngle;
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
    dsPath = sprintf('%s/%s%s/', rt, DSDirName, surffix);
    if ~exist(dsPath, 'dir')
        mkdir(dsPath);
        if ~ispc
            fileattrib(dsPath, '+w', 'g');
        end
    end
    if saveZarr
        dsFullname = [dsPath, fsname, '.zarr'];
    else
        dsFullname = [dsPath, fsname, '.tif'];
    end
end

if (~DSRCombined && (~exist(dsFullname, 'file') || ip.Results.overwrite)) || DSRCombined
    t_read = tic;
    % frame = double(readtiff(framePath));
    if combinedFrame
        frame_cell = cell(numel(framePath), 1);
        for i = 1 : numel(framePath)
            frame_cell{i} = readtiff(framePath{i});
        end
        frame = cat(3, frame_cell{:});
        clear frame_cell;
    else
        [~, ~, ext] = fileparts(framePath{1});
        switch ext
            case {'.tif', '.tiff'}
                frame = readtiff(framePath{1});
            case {'.zarr'}
                frame = readzarr(framePath{1});
        end                
    end
    t_read = toc(t_read);
    fprintf('Image read time: %0.6f s\n', t_read);

    if ~isempty(inputBbox)
        try 
            frame = crop3d_mex(frame, inputBbox);
        catch ME
            disp(ME);
            frame = frame(inputBbox(1) : inputBbox(4), inputBbox(2) : inputBbox(5), inputBbox(3) : inputBbox(6));
        end
    end

    if flipZstack
        frame = flip(frame, 3);
    end

    % flat field correction
    if FFCorrection
        fprintf(['Flat-field correction for frame %s...\n', ...
            '  Flat-field image: %s\n  Background image: %s\n'], ...
            framePath{1}, FFImage, backgroundImage);
        LSIm = readtiff(FFImage);
        BKIm = readtiff(backgroundImage);            
        frame = XR_LSFlatFieldCorrection(frame,LSIm,BKIm,'lowerLimit', lowerLimit, ...
            'constOffset', constOffset);
    end
    % remove camera background
    if BKRemoval
        BKIm = readtiff(backgroundImage);
        frame = XR_CameraBackgroundRemoval(frame, BKIm, 'constOffset', constOffset);
    end

    if ~DSRCombined
        fprintf('Deskew frame %s...\n', framePath{1});
        splitCompute = false;
        sz = size(frame);
        % 03/23/2021, for image with more than 1000 slices, directly split to blocks for the processing
        if sz(3) > 1000
            splitCompute = true;
        end

        frame = single(frame);
        if ~splitCompute
            try 
                ds = deskewFrame3D(frame, skewAngle_1, dz, xyPixelSize, reverse, ...
                    'crop', crop, 'Interp', interpMethod); 
                clear frame;
            catch ME
                disp(ME);
                sz = size(frame);
                if sz(3) > 200
                    splitCompute = true;
                end
            end
        end

        if splitCompute
            fprintf('Use image block method for deskew...\n');
            % only split in y-axis, it is not right when splitting from
            % other axes. 
            blockSize = sz;
            blockSize(1) = min(ceil(blockSize(1) / 8), 150);
            bim = blockedImage(frame, 'blockSize', blockSize);
            OutputLocation = sprintf('%s/%s_%s', dsPath, fsname, uuid);
            BorderSize = [5, 0, 0];
            % TrimBorder = true;
            if save16bit 
                bo = apply(bim, @(bs) uint16(trimBorder(deskewFrame3D(single(bs.Data), skewAngle_1, dz, ...
                    xyPixelSize, reverse, 'crop', crop, 'Interp', interpMethod), BorderSize, 'both')), 'blockSize', bim.BlockSize, ...
                    'OutputLocation', OutputLocation, 'BorderSize', BorderSize, 'useParallel', false);
            else
                bo = apply(bim, @(bs) trimBorder(deskewFrame3D(single(bs.Data), skewAngle_1, dz, ...
                    xyPixelSize, reverse, 'crop', crop, 'Interp', interpMethod), BorderSize, 'both'), 'blockSize', bim.BlockSize, ...
                    'OutputLocation', OutputLocation, 'BorderSize', BorderSize, 'useParallel', false);
            end
            clear frame bim;
            % ds = gather(bo);
            % rmdir(OutputLocation, 's');
        end            

        % save MIP
        if saveMIP
            dsMIPPath = sprintf('%s/MIPs/', dsPath);
            if ~exist(dsMIPPath, 'dir')
                mkdir(dsMIPPath);
                if ~ispc                
                    fileattrib(dsMIPPath, '+w', 'g');
                end
            end

            if splitCompute
                bmip = apply(bo, @(bs) max(bs.Data, [], 3), 'blockSize', [bo.blockSize(1:2), bo.Size(3)], 'useParallel', false);
                mip = gather(bmip);
            else
                mip = max(ds, [], 3);
            end
            if save16bit 
                mip = uint16(mip);
            end
            if saveZarr
                dsMIPname = sprintf('%s%s_MIP_z.zarr', dsMIPPath, fsname);
                writezarr(mip, dsMIPname);
            else
                dsMIPname = sprintf('%s%s_MIP_z.tif', dsMIPPath, fsname);
                writetiff(mip, dsMIPname);
            end
        end

        dsTempname = sprintf('%s%s_%s.tif', dsPath, fsname, uuid);
        if save3DStack
            if splitCompute
                if saveZarr
                    write(bo, dsTempname, 'blockSize', min(bo.Size, blockSize), 'Adapter', ZarrAdapter);
                else
                    write(bo, dsTempname, 'blockSize', [bo.Size(1), bo.Size(2), min(bo.Size(3), 100)], 'Adapter', MPageTiffAdapter);
                end
                rmdir(OutputLocation, 's');
                clear bo
            else                
                if save16bit
                    if saveZarr
                        writezarr(uint16(ds), dsTempname, 'blockSize', blockSize);
                    else
                        writetiff(uint16(ds), dsTempname);
                    end
                else
                    if saveZarr
                        writezarr(single(ds), dsTempname, 'blockSize', blockSize);
                    else
                        writetiff(single(ds), dsTempname);
                    end
                end
            end
            movefile(dsTempname, dsFullname);
        else
            fclose(fopen(dsTempname, 'w'));
            movefile(dsTempname, dsFullname);            
        end
    end
end

%% rotate frame

if ip.Results.rotate || DSRCombined
    dsrPath = sprintf('%s/%s%s/', rt, DSRDirName, surffix);
    if ~exist(dsrPath, 'dir')
        mkdir(dsrPath);
        if ~ispc
            fileattrib(dsrPath, '+w', 'g');
        end
    end

    if saveZarr
        dsrFullname = [dsrPath, fsname, '.zarr'];        
    else
        dsrFullname = [dsrPath, fsname, '.tif'];
    end

    if ~isempty(resampleFactor)
        rs = resampleFactor(:)';
        % complete rs to 3d in case it is not
        resampleFactor = [ones(1, 4 - numel(rs)) * rs(1), rs(2:end)];    
    end

    if (~saveZarr && ~exist(dsrFullname, 'file')) || (saveZarr && ~exist(dsrFullname, 'dir'))
        t_dsr = tic;
        if ~DSRCombined
            fprintf('rotate frame %s...\n', framePath{1});
            if ~exist('ds', 'var')
                if saveZarr
                    ds = single(readzarr(dsFullname));
                else
                    ds = single(readtiff(dsFullname));
                end
                sz = getImageSize(framePath{1});
                for i = 2 : numel(framePath)
                    sz_i = getImageSize(framePath{i});
                    sz(3) = sz(3) + sz_i(3);
                end
            end

            ny = sz(1);
            nx = sz(2);
            nz = sz(3);
            if ~objectiveScan
                % calculate height; first & last 2 frames have interpMethodolation artifacts
                outSize = round([ny, (nx-1)*cos(theta)+(nz-1)*zAniso/sin(abs(theta)), (nx-1)*sin(abs(theta))-4]);
            else
                % exact proportions of rotated box
                outSize = round([ny, nx*cos(theta)+nz*zAniso*sin(abs(theta)), nz*zAniso*cos(theta)+nx*sin(abs(theta))]);
            end

            dsr = rotateFrame3D(ds, skewAngle_1, zAniso, reverse, 'crop', true, ...
                'resample', resampleFactor, 'objectiveScan', objectiveScan, 'outSize', outSize, 'Interp', interpMethod);

            t_dsr = toc(t_dsr);
            fprintf('Deskew/rotation processing time: %0.6f s\n', t_dsr);

            if nargout == 0
                clear ds;
            end
        else
            % add support for resample before dsr for big data
            if ~isempty(resampleFactor) && any(resampleFactor ~= 1) && size(frame, 3) > 1000
                outPixelSize = resampleFactor * xyPixelSize;
                pre_rs = min(outPixelSize) ./ [xyPixelSize, xyPixelSize, dz .* sind(skewAngle_1)];
                pre_rs(3) = max(1, round(pre_rs(3)));
                
                % resample data pre-DSR
                % binning in z and resize in xy
                if any(pre_rs ~= 1)
                    frame = frame(:, :, round(pre_rs(3) / 2) : pre_rs(3) : end);
                    frame = imresize(frame, round(size(frame, [1, 2]) ./ pre_rs(1 : 2)));

                    % define new xyPixelSize, dz, resample after resample
                    xyPixelSize = xyPixelSize * pre_rs(1);
                    dz = dz * pre_rs(3);
                    resampleFactor = outPixelSize ./ min(outPixelSize);
                end
            end
            if ~save16bit
                frame = single(frame);
            end

            fprintf('Deskew, rotate and resample for frame %s...\n', framePath{1});
            dsr = deskewRotateFrame3D(frame, skewAngle_1, dz, xyPixelSize, ...
                'reverse', reverse, 'crop', true, 'objectiveScan', objectiveScan, ...
                'resampleFactor', resampleFactor, 'interpMethod', interpMethod, ...
                'xStepThresh', xStepThresh, 'save16bit', save16bit);

            t_dsr = toc(t_dsr);
            fprintf('Combined deskew/rotation processing time: %0.6f s\n', t_dsr);

            clear frame;
        end
        
        % save MIP
        if saveMIP
            dsrMIPPath = sprintf('%s/MIPs/', dsrPath);
            if ~exist(dsrMIPPath, 'dir')
                mkdir(dsrMIPPath);
                if ~ispc                
                    fileattrib(dsrMIPPath, '+w', 'g');
                end
            end

            mip = max(dsr, [], 3);            
            if save16bit
                mip = uint16(mip);
            end
            if saveZarr
                dsrMIPname = sprintf('%s%s_MIP_z.zarr', dsrMIPPath, fsname);
                dsrMIPTmpname = sprintf('%s%s_MIP_z.zarr_%s', dsrMIPPath, fsname, uuid);
                writezarr(mip, dsrMIPTmpname);
            else
                dsrMIPname = sprintf('%s%s_MIP_z.tif', dsrMIPPath, fsname);
                dsrMIPTmpname = sprintf('%s%s_MIP_z.tiff_%s', dsrMIPPath, fsname, uuid);
                writetiff(mip, dsrMIPTmpname);
            end
            movefile(dsrMIPTmpname, dsrMIPname);
        end
        
        t_write = tic;
        if saveZarr
            dsrTempName = sprintf('%s%s_%s.zarr', dsrPath, fsname, uuid);
        else
            dsrTempName = sprintf('%s%s_%s.tif', dsrPath, fsname, uuid);            
        end
        
        if ip.Results.rescaleRotate
            iRange = [min(dsr(:)), max(dsr(:))];
            dsr = scaleContrast(dsr, iRange, [0 65535]);
        end
        
        if save16bit
            dsr = uint16(dsr);
        else
            dsr = single(dsr);
        end
        
        if save3DStack
            if saveZarr
                writezarr(dsr, dsrTempName, 'blockSize', blockSize);
            else
                writetiff(dsr, dsrTempName);
            end
            movefile(dsrTempName, dsrFullname);
        else
            fclose(fopen(dsrTempName, 'w'));
            movefile(dsrTempName, dsrFullname);            
        end
        t_write = toc(t_write);
        fprintf('Deskew/rotation result write time: %0.6f s\n', t_write);
    end
end

end

