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

function [volout] = rotateFrame3D(vol, angle, zxRatio, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol');
ip.addRequired('angle'); % typical value: 32.8
ip.addRequired('zxRatio'); % typical value: ~ 2 (sample scan) - 4 (obj. scan)
ip.addOptional('reverse', false, @islogical);
ip.addParamValue('Crop', true, @islogical);
ip.addParamValue('ObjectiveScan', false, @islogical);
ip.addParamValue('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
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
    
    startIndex = ones(nz,1);
    endIndex = nx*ones(nz,1);
    [srow,scol] = find(proj==1);
    [erow,ecol] = find(proj==-1);
    startIndex(srow) = scol;
    endIndex(erow) = ecol;
    
    doCrop = any(startIndex>1 & endIndex<nx);
    
    % calculate height; first & last 2 frames have interpolation artifacts
    w = median(diff([startIndex endIndex],[],2));
    h = round(abs(w*tan(theta)*cos(theta)))-4;
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

if ip.Results.Crop && doCrop
    %outSize = round([ny nx/cos(theta) nz*ip.Results.zxRatio]);
    outSize = round([ny nx/cos(theta) h]);
else
    % exact proportions of rotated box
    outSize = round([ny nx*cos(theta)+nz*ip.Results.zxRatio*sin(abs(theta)) nz*ip.Results.zxRatio*cos(theta)+nx*sin(abs(theta))]);
end

T2 = [1 0 0 0
      0 1 0 0
      0 0 1 0
      (outSize([2 1 3])+1)/2 1];

RA = imref3d(outSize, 1, 1, 1);
[volout] = imwarp(vol, affine3d(T1*S*R*T2), ip.Results.Interp, 'FillValues', 0, 'OutputView', RA);
