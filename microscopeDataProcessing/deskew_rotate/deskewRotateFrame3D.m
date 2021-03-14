function volout = deskewRotateFrame3D(vol, angle, dz, xyPixelSize, varargin)
% Applies a shear transform to convert raw light sheet microscopy data
% into a volume with real-world coordinates. After deskew, rotate the view of
% the volume, followed by resampling (optional). 
% 
% Based on deskewFrame3D.m and rotateFrame3D.m, and also add resampling step.
% 
% Author: Xiongtao Ruan (10/08/2020)

% xruan (02/10/2021): add option to directly apply combined processing when x step size 
% is small, and use separate processing when x step size is large (also split to
% parts in the processing when the image is tall. 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol');
ip.addRequired('angle'); % typical value: 32.8
ip.addRequired('dz'); % typical value: 0.2-0.5
ip.addRequired('xyPixelSize'); % typical value: 0.1
ip.addOptional('reverse', false, @islogical);
ip.addParameter('Crop', true, @islogical);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('xStepThresh', 3.2, @isnumeric); % 3.125 for ds=0.4
ip.addParameter('resample', [], @isnumeric); % resample factor in xyz order. 
ip.addParameter('gpuProcess', false, @islogical); % use gpu for the processing. 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.parse(vol, angle, dz, xyPixelSize, varargin{:});

[ny,nx,nz] = size(vol);

theta = ip.Results.angle * pi/180;
dx = cos(theta)*dz/xyPixelSize; % pixels shifted slice to slice in x

if ip.Results.ObjectiveScan
    zAniso = dz / xyPixelSize;
else
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
end

%% deskew
if ~ip.Results.reverse
    xshift = -dx;
    xstep = dx;
else
    xshift = dx + ceil((nz-1)*dx);
    xstep = -dx;
end
nxDs = ceil((nz-1)*dx) + nx; % width of output volume as if there is DS.

% shear transform matrix
if ip.Results.ObjectiveScan
    nxDs = nx;
    ds_S = eye(4);
else
    ds_S = [1 0 0 0;
            0 1 0 0;
            xstep 0 1 0;
            xshift 0 0 1];
end

%% rotate
% nxDs = nxOut;
if ip.Results.reverse
    theta = -theta;
end

center = ([ny nxDs nz]+1)/2;
T1 = [1 0 0 0
      0 1 0 0
      0 0 1 0
      -center([2 1 3]) 1];

S = [1 0 0 0
     0 1 0 0
     0 0 zAniso 0
     0 0 0 1];

% Rotate x,z
R = [cos(theta) 0 -sin(theta) 0; % order for imwarp is x,y,z
     0 1 0 0;
     sin(theta) 0 cos(theta) 0;
     0 0 0 1];

if ~ip.Results.ObjectiveScan
    % outSize = round([ny nxDs/cos(theta) h]);
    % calculate height; first & last 2 frames have interpolation artifacts
    outSize = round([ny, (nx-1)*cos(theta)+(nz-1)*zAniso/sin(abs(theta)), (nx-1)*sin(abs(theta))-4]);
else
    % exact proportions of rotated box
    outSize = round([ny, nx*cos(theta)+nz*zAniso*sin(abs(theta)), nz*zAniso*cos(theta)+nx*sin(abs(theta))]);
end

T2 = [1 0 0 0
      0 1 0 0
      0 0 1 0
      (outSize([2 1 3])+1)/2 1];

%% resampling after deskew and rotate
rs = ip.Results.resample;
if ~isempty(rs)
    RT1 = [1 0 0 0
           0 1 0 0
           0 0 1 0
           -(outSize([2,1,3])+1)/2 1];
    RS =[1/rs(1) 0 0 0
         0 1/rs(2) 0 0
         0 0 1/rs(3) 0
         0 0 0 1];
    outSize = round(outSize ./ rs([2,1,3]));
    RT2 = [1 0 0 0
           0 1 0 0
           0 0 1 0
           (outSize([2,1,3])+1)/2 1];     
else
    RT1 = eye(4);
    RS = eye(4);
    RT2 = eye(4);
end

%% summarized transform
RA = imref3d(outSize, 1, 1, 1);
if ip.Results.gpuProcess
    vol = gpuArray(vol);
end

% for sample scan, check x step size to decide use separate or combined
% processing
combinedProcess = true;
if ~ip.Results.ObjectiveScan && abs(dx) > ip.Results.xStepThresh
    combinedProcess = false;
end

if combinedProcess
    [volout] = imwarp(vol, affine3d(ds_S*(T1*S*R*T2)*(RT1*RS*RT2)), ip.Results.Interp, 'FillValues', 0, 'OutputView', RA);
else
    try 
        [vol] = imwarp(vol, affine3d(ds_S), ip.Results.Interp, 'FillValues', 0);
        [volout] = imwarp(vol, affine3d((T1*S*R*T2)*(RT1*RS*RT2)), ip.Results.Interp, 'FillValues', 0, 'OutputView', RA);
    catch ME
        % for very tall images
        % parameters in y axis
        padSize = 5; 
        batchSize = 100;
        
        nBatch = ceil(ny / batchSize);
        volout = zeros(outSize, class(vol));
        for n = 1 : nBatch
            inds_n = (n - 1) * batchSize + 1 : min(n * batchSize, ny);
            s = max(1, inds_n(1) - padSize); 
            t = min(ny, inds_n(end) + padSize); 
            vol_n = vol(s : t, :, :);
            
            outSize_n = [t - s + 1, outSize(2), outSize(3)];
            RA_n = imref3d(outSize_n, 1, 1, 1);
            [vol_n] = imwarp(vol_n, affine3d(ds_S), ip.Results.Interp, 'FillValues', 0);
            [vol_n] = imwarp(vol_n, affine3d((T1*S*R*T2)*(RT1*RS*RT2)), ip.Results.Interp, 'FillValues', 0, 'OutputView', RA_n);
            
            volout(inds_n, :, :) = vol_n(inds_n(1) - s + 1 : inds_n(end) - s + 1, :, :);
        end
    end
end

% test separate affine transform for DS and dSR
% if false
%     [vol] = imwarp(vol, affine3d(ds_S), ip.Results.Interp, 'FillValues', 0);
%     [volout] = imwarp(vol, affine3d((T1*S*R*T2)*(RT1*RS*RT2)), ip.Results.Interp, 'FillValues', 0, 'OutputView', RA);
% end
if ip.Results.gpuProcess
    volout = gather(volout);
end


end
