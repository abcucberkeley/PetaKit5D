%[volout] = rotateFrame3D(vol, angle, zxRatio, varargin) rotates de-skewed light sheet microscope frames for visualization
%
% Inputs: 
%           vol : deskewed volume (or raw volume for objective scan data)
%         angle : scanning angle of the system (in degrees)
%       zxRatio : ratio between true z-step and pixel size
%                 This is dz*sin(angle)/psize for sample scan data,
%                 and dz/psize for objective scan data, where dz is the step size
%
% Options:
%       reverse : true|{false}, selects scanning direction
%                 true: top right to bottom left
%                 false: top left to bottom right
%
% Parameters:
%       'Crop' : Automatically crop empty slices at top/bottom of rotated volume.
%                Default: true
%     'Interp' : 'cubic'|{'linear'} interpolation mode.
%
% See also imwarp, deskewFrame3D

% Author: Francois Aguet
% xruan (11/03/2021): decide crop size by the distance of two boundary
% lines for rotation
% xruan (11/30/2021): add support for resampling
% xruan (12/06/2022): add support for user input outSize

function [volout] = rotateFrame3D(vol, angle, zxRatio, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol');
ip.addRequired('angle'); % typical value: 32.8
ip.addRequired('zxRatio'); % typical value: ~ 2 (sample scan) - 4 (obj. scan)
ip.addOptional('reverse', false, @islogical);
ip.addParameter('Crop', true, @islogical);
ip.addParameter('resample', [], @isnumeric); % resample factor in xyz order. 
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('outSize', [], @isnumeric);
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.parse(vol, angle, zxRatio, varargin{:});

[ny,nx,nz] = size(vol);

theta = ip.Results.angle * pi/180;
if ip.Results.reverse
    theta = -theta;
end

doCrop = false;
if ~ip.Results.ObjectiveScan
    % calculate height of rotated volume
    % project and find left/right edges of deskewed rhomboid, take median width
    proj = squeeze(max(vol,[],1))'~=0;
    proj = diff(proj,1,2);
        
    % xruan: use distance between boundary lines before rotation to decide
    % the height after rotation 
    
    % boarder point coordinates
    [z, x] = find(proj ~= 0);
    % crop size in z
    a = sin(-theta);
    b = cos(-theta);
    c = a * x + b * z * zxRatio; 
    h = round(max(c) - min(c)) - 4;
    
    % use the bounding box to decide whether to crop or not
    doCrop = h < max(abs(a * (nx - 1) + b * (nz - 1) * zxRatio), abs(a * (nx - 1) + b * (1 - nz) * zxRatio));
    if isempty(doCrop)
        doCrop = false;
    end
        
    % crop size in x
    c_x = b * x + (-a) * z * zxRatio;
    w = round(max(c_x) - min(c_x)) - 4;
end

center = ([ny nx nz]+1)/2;
T1 = [1 0 0 0
      0 1 0 0
      0 0 1 0
      -center([2 1 3]) 1];

S = [1 0 0 0
     0 1 0 0
     0 0 ip.Results.zxRatio 0
     0 0 0 1];

% Rotate x,z
R = [cos(theta) 0 -sin(theta) 0; % order for imwarp is x,y,z
     0 1 0 0;
     sin(theta) 0 cos(theta) 0;
     0 0 0 1];

if isempty(ip.Results.outSize)
    if ip.Results.Crop && doCrop
        %outSize = round([ny nx/cos(theta) nz*ip.Results.zxRatio]);
        % outSize = round([ny nx/cos(theta) h]);
        outSize = round([ny w h]);
    else
        % exact proportions of rotated box
        outSize = round([ny nx*cos(theta)+nz*ip.Results.zxRatio/sin(abs(theta)) nz*ip.Results.zxRatio*cos(theta)+nx*sin(abs(theta))]);
    end
else
    outSize = ip.Results.outSize;
end

T2 = [1 0 0 0
      0 1 0 0
      0 0 1 0
      (outSize([2 1 3])+1)/2 1];
  
% resampling after deskew and rotate
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

RA = imref3d(outSize, 1, 1, 1);
[volout] = imwarp(vol, affine3d((T1*S*R*T2)*(RT1*RS*RT2)), ip.Results.Interp, 'FillValues', 0, 'OutputView', RA);

end

