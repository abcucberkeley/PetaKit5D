function [vol] = deskewRotateFrame3D(vol, angle, dz, xyPixelSize, varargin)
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
% xruan (03/16/2021): change default xStepThresh to 2.35 (ds=0.3). 
% xruan (01/27/2022): change default xStepThresh to 2.74 (ds=0.35). 
% xruan (05/26/2022): change default xStepThresh to 2.42 (ds=0.31). 
% xruan (05/30/2022): add skewed space interpolation based dsr for large step size
% xruan (06/02/2022): change default xStepThresh to 1.96 (ds=0.25). 
% xruan (07/17/2022): change default xStepThresh to 2.00 (ds=0.255). 
% xruan (10/29/2022): refactor code to first decide whether interpolate in skewed space 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol');
ip.addRequired('angle'); % typical value: 32.8
ip.addRequired('dz'); % typical value: 0.2-0.5
ip.addRequired('xyPixelSize'); % typical value: 0.1
ip.addOptional('reverse', false, @islogical);
ip.addParameter('Crop', true, @islogical);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('xStepThresh', 2.0, @isnumeric); % 2.344 for ds=0.3, 2.735 for ds=0.35
ip.addParameter('resample', [], @isnumeric); % resample factor in xyz order. 
ip.addParameter('gpuProcess', false, @islogical); % use gpu for the processing. 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.parse(vol, angle, dz, xyPixelSize, varargin{:});

pr = ip.Results;
Reverse = pr.reverse;
Crop = pr.Crop;
ObjectiveScan = pr.ObjectiveScan;
xStepThresh = pr.xStepThresh;
resample = pr.resample;
gpuProcess = pr.gpuProcess;
Interp = pr.Interp;

[ny,nx,nz] = size(vol);


theta = angle * pi/180;
dx = cos(theta)*dz/xyPixelSize; % pixels shifted slice to slice in x

if ip.Results.ObjectiveScan
    zAniso = dz / xyPixelSize;
else
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
end

% use original dz to decide outSize
if ~ObjectiveScan
    % outSize = round([ny nxDs/cos(theta) h]);
    % calculate height; first & last 2 frames have interpolation artifacts
    outSize = round([ny, (nx-1)*cos(theta)+(nz-1)*zAniso/sin(abs(theta)), (nx-1)*sin(abs(theta))-4]);
else
    % exact proportions of rotated box
    outSize = round([ny, nx*cos(theta)+nz*zAniso*sin(abs(theta)), nz*zAniso*cos(theta)+nx*sin(abs(theta))]);
end

%% skew space interpolation
if ~ObjectiveScan && abs(dx) > xStepThresh
    % skewed space interplation combined dsr
    fprintf('The step size is greater than the threshold, use skewed space interpolation for combined deskew rotate...\n');
    % for dx only slightly larger than xStepThresh, we interpolate to
    % isotropic (approximate), so it is not oversampled a lot to save memory. 
    if abs(dx) / xStepThresh < 1.5
        dzout = xyPixelSize / sin(theta);
        % round it by the significant digits of dz
        sf = 10^floor(log10(dz));
        dzout = round(dzout / sf ) * sf;
    else
        ndiv = ceil(abs(dx) / xStepThresh);
        dzout = dz / ndiv;
    end
    int_stepsize = dzout / dz;
    fprintf('Input dz: %f , interpolated dz: %f\n', dz, dzout);

    % add the mex version skewed space interpolation as default
    try 
        vol = skewed_space_interp_defined_stepsize_mex(vol, abs(dx), int_stepsize, Reverse);
    catch ME
        disp(ME);
        vol = skewed_space_interp_defined_stepsize(vol, abs(dx), int_stepsize, 'Reverse', Reverse);
    end
    
    % update parameters after interpolation
    [ny,nx,nz] = size(vol);
    dz = dzout;
    dx = cos(theta)*dz/xyPixelSize; % pixels shifted slice to slice in x

    if ip.Results.ObjectiveScan
        zAniso = dz / xyPixelSize;
    else
        zAniso = sin(abs(theta)) * dz / xyPixelSize;
    end
end

%% deskew
if ~Reverse
    xshift = -dx;
    xstep = dx;
else
    xshift = dx + ceil((nz-1)*dx);
    xstep = -dx;
end
nxDs = ceil((nz-1)*dx) + nx; % width of output volume as if there is DS.

% shear transform matrix
if ObjectiveScan
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
if Reverse
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

T2 = [1 0 0 0
      0 1 0 0
      0 0 1 0
      (outSize([2 1 3])+1)/2 1];

%% resampling after deskew and rotate
rs = resample;
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
if gpuProcess
    vol = gpuArray(vol);
end

[vol] = imwarp(vol, affine3d(ds_S*(T1*S*R*T2)*(RT1*RS*RT2)), Interp, 'FillValues', 0, 'OutputView', RA);
if gpuProcess
    vol = gather(vol);
end    


end





