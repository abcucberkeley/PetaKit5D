function [volout] = deskewRotateFrame3D(vol, angle, dz, xyPixelSize, varargin)
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
% xruan (08/02/2023): add support for direct bounding box crop of output


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol');
ip.addRequired('angle'); % typical value: 32.8
ip.addRequired('dz'); % typical value: 0.2-0.5
ip.addRequired('xyPixelSize'); % typical value: 0.1
ip.addParameter('reverse', false, @islogical);
ip.addParameter('Crop', true, @islogical);
ip.addParameter('bbox', [], @isnumeric);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('xStepThresh', 2.0, @isnumeric); % 2.344 for ds=0.3, 2.735 for ds=0.35
ip.addParameter('resampleFactor', [], @isnumeric); % resample factor after rotation, [x y z] (RS scales imwarp x first); [1 1 sin(angle)*dz/xyPixelSize] keeps the raw sampling density
ip.addParameter('gpuProcess', false, @islogical); % use gpu for the processing. 
ip.addParameter('save16bit', false, @islogical); % direct output results as 16bit for mex functions
ip.addParameter('interpMethod', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.parse(vol, angle, dz, xyPixelSize, varargin{:});

pr = ip.Results;
Reverse = pr.reverse;
bbox = pr.bbox;
objectiveScan = pr.objectiveScan;
xStepThresh = pr.xStepThresh;
resampleFactor = pr.resampleFactor;
gpuProcess = pr.gpuProcess;
save16bit = pr.save16bit;
interpMethod = pr.interpMethod;

[ny,nx,nz] = size(vol);

theta = angle * pi/180;
dx = cos(theta)*dz/xyPixelSize; % pixels shifted slice to slice in x

if ip.Results.objectiveScan
    zAniso = dz / xyPixelSize;
else
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
end

% use original dz to decide outSize
if ~objectiveScan
    % outSize = round([ny nxDs/cos(theta) h]);
    % calculate height; first & last 2 frames have interpolation artifacts
    outSize = round([ny, (nx-1)*cos(theta)+(nz-1)*zAniso/sin(abs(theta)), (nx-1)*sin(abs(theta))-4]);
else
    % exact proportions of rotated box
    outSize = round([ny, nx*cos(theta)+nz*zAniso*sin(abs(theta)), nz*zAniso*cos(theta)+nx*sin(abs(theta))]);
end

do_interp = ~objectiveScan && abs(dx) > xStepThresh;
rs = resampleFactor;
if ~isempty(rs)
    validateattributes(rs, {'numeric'}, {'positive', 'finite'}, mfilename, 'resampleFactor');
end
nz_in = nz;
dz_in = dz;
%% skew space interpolation parameters (the volume is only materialized on the imwarp path)
if do_interp
    % skewed space interplation combined dsr
    % for dx only slightly larger than xStepThresh, we interpolate to
    % a step size lower than the threshold, and the ratio between dz /
    % dzout ceils to the faction of 1/n with 10% to the threshold.
    if abs(dx) / xStepThresh < 1.5
        % dzout = xyPixelSize / sin(theta);
        dzout_thresh = xyPixelSize * xStepThresh / cos(theta);
        dzout = dz / (ceil((dz / dzout_thresh) / (1 / 1)) *(1 / 1));
        counter = 1;
        while dzout / dzout_thresh < 0.4
            dzout = dz / (ceil((dz / dzout_thresh) / (1 / (counter + 1))) *(1 / (counter + 1)));
            counter = counter + 1;
        end
        % round it by the significant digits of dz
        % sf = 10^floor(log10(dz));
        % dzout = round(dzout / sf ) * sf;
        % dzout = dz / 2;
    else
        ndiv = ceil(abs(dx) / xStepThresh);
        dzout = dz / ndiv;
    end
    int_stepsize = dzout / dz;

    % geometry after interpolation; the one definition of the plane count (see interpolatedPlanes)
    nz = interpolatedPlanes(nz, int_stepsize);
    dx_orig = dx;
    dz = dzout;
    dx = cos(theta)*dz/xyPixelSize; % pixels shifted slice to slice in x

    if ip.Results.objectiveScan
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
if objectiveScan
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
% M is the forward transform for imwarp (1-based pixel centres). The mex takes the same
% map as a backward matrix in 0-based indices and its own axis order: one conjugation
% by P converts the whole map, so the shear needs no index-dependent adjustment.
M = ds_S*(T1*S*R*T2)*(RT1*RS*RT2);
P = eye(4);
P(4, 1 : 3) = 1;
tmat = eye(4) / (P*M/P)';
tmat = tmat([2, 1, 3, 4], [2, 1, 3, 4]);

% the warp mex implements one affine class (dsr_mex_accepts) for single/uint16 input;
% everything else goes through imwarp with M
in_class = dsr_mex_accepts(tmat);
dtype_ok = isa(vol, 'single') || isa(vol, 'uint16');
linear = strcmpi(interpMethod, 'linear');
use_fast_method = in_class && dtype_ok && linear && ~gpuProcess;

% one line that reproduces the decision from a log: geometry in, geometry out, and the gates
if do_interp
    dz_str = sprintf('%.4g -> %.4g (skew-space interpolation x%d)', dz_in, dz, round(1 / int_stepsize));
else
    dz_str = sprintf('%.4g', dz_in);
end
if isempty(rs)
    rs_str = '[]';
else
    rs_str = mat2str(rs, 4);
end
fprintf('Deskew/rotate: input %dx%dx%d %s, dz %s, xy %.4g, skew %.4g, reverse %d, objectiveScan %d, resampleFactor %s, output %s\n', ...
    ny, nx, nz_in, class(vol), dz_str, xyPixelSize, angle, Reverse, objectiveScan, rs_str, mat2str(outSize));
if use_fast_method
    if do_interp
        fprintf('Deskew/rotate: combined mex path (fused skew-space interpolation + warp).\n');
    else
        fprintf('Deskew/rotate: combined mex path (warp).\n');
    end
else
    fprintf('Deskew/rotate: imwarp path (in mex class: %d, dtype ok: %d, linear: %d, gpu: %d). Backward map (mex axes, 0-based):\n', ...
        in_class, dtype_ok, linear, gpuProcess);
    disp(round(tmat, 4));
end

RA = imref3d(outSize, 1, 1, 1);
if ~isempty(bbox)
    RA = imref3d(bbox(4 : 6) - bbox(1 : 3) + 1, [bbox(2) - 0.5, bbox(5) + 0.5], [bbox(1) - 0.5, bbox(4) + 0.5], [bbox(3) - 0.5, bbox(6) + 0.5]);
end

if ~use_fast_method
    vol_1 = vol;
    if do_interp
        t0 = tic;
        % the interpolation mex handles single/uint16; other types use the MATLAB version
        if isa(vol, 'single') || isa(vol, 'uint16')
            try
                vol_1 = skewed_space_interp_defined_stepsize_mex(vol, abs(dx_orig), int_stepsize, Reverse, save16bit);
            catch ME
                disp(ME);
                vol_1 = skewed_space_interp_defined_stepsize(vol, abs(dx_orig), int_stepsize, 'Reverse', Reverse);
            end
        else
            vol_1 = skewed_space_interp_defined_stepsize(vol, abs(dx_orig), int_stepsize, 'Reverse', Reverse);
        end
        fprintf('Skewed space interpolation time: %f s\n', toc(t0));
        % M was built for interpolatedPlanes(nz_in); a materializer that disagrees would shift the output
        assert(size(vol_1, 3) == nz, 'deskewRotateFrame3D:planeCount', ...
            'skewed-space interpolation produced %d planes, transform expects %d (nz_in %d, step %.6g)', ...
            size(vol_1, 3), nz, nz_in, int_stepsize);
    end
    if gpuProcess
        vol_1 = gpuArray(vol_1);
    end
    [volout] = imwarp(vol_1, affine3d(M), interpMethod, 'FillValues', 0, 'OutputView', RA);
else
    if ~isempty(bbox)
        bbox_in = bbox;
    else
        bbox_in = [1, 1, 1, outSize];
    end

    if do_interp
        try
            volout = skewed_space_interp_volume_deskew_rotate_warp_mex(vol, abs(dx_orig), int_stepsize, Reverse, tmat, bbox_in, save16bit);
        catch ME
            fprintf('Deskew/rotate: fused mex failed, fall back to skew-space interpolation + warp mex.\n');
            disp(ME);
            try
                vol_1 = skewed_space_interp_defined_stepsize_mex(vol, abs(dx_orig), int_stepsize, Reverse, save16bit);
            catch ME
                fprintf('Deskew/rotate: interpolation mex failed, fall back to the MATLAB interpolation.\n');
                disp(ME);
                vol_1 = skewed_space_interp_defined_stepsize(vol, abs(dx_orig), int_stepsize, 'Reverse', Reverse);
            end
            try
                volout = volume_deskew_rotate_warp_mex(vol_1, tmat, bbox_in, save16bit);
            catch ME
                fprintf('Deskew/rotate: warp mex failed, fall back to imwarp.\n');
                disp(ME);
                [volout] = imwarp(vol_1, affine3d(M), interpMethod, 'FillValues', 0, 'OutputView', RA);
            end
        end
    else
        try
            volout = volume_deskew_rotate_warp_mex(vol, tmat, bbox_in, save16bit);
        catch ME
            fprintf('Deskew/rotate: warp mex failed, fall back to imwarp.\n');
            disp(ME);
            [volout] = imwarp(vol, affine3d(M), interpMethod, 'FillValues', 0, 'OutputView', RA);
        end
    end
end
if gpuProcess
    volout = gather(volout);
end

end


function n = interpolatedPlanes(nz, int_stepsize)
% Plane count after skewed-space interpolation to step int_stepsize (fraction of dz). This is
% the arithmetic the interpolation mex uses (rounded at 1e-5 against float slop) and the one
% definition the transform is built from; the imwarp path asserts the materialized volume matches.
n = floor(round((nz - 1) / int_stepsize * 100000) / 100000) + 1;
end


