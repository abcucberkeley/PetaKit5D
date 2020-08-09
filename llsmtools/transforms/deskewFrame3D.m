%deskewFrame3D(vol, angle, dz, xyPixelSize, varargin)
% Applies a shear transform to convert raw light sheet microscopy data
% into a volume with real-world coordinates.
%
% Inputs:
%           vol : raw data volume
%         angle : scanning angle of the system (in degrees)
%            dz : distance between acquisition planes (z-step), in µm
%   xyPixelSize : pixel size in object space, in µm
%
% Options:
%       reverse : {true}|false, selects scanning direction
%                 true: top right to bottom left
%                 false: top left to bottom right
%
% Parameters:
%        'Crop' : {true}|false crops the transformed volume at 50% of the
%                 outermost image planes
%      'Interp' : {'cubic'}|'linear' interpolation mode.
%
% See also imwarp

% Author: Francois Aguet (09/2013)

function volout = deskewFrame3D(vol, angle, dz, xyPixelSize, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol');
ip.addRequired('angle'); % typical value: 32.8
ip.addRequired('dz'); % typical value: 0.2-0.5
ip.addRequired('xyPixelSize'); % typical value: 0.1
ip.addOptional('reverse', false, @islogical);
ip.addParamValue('Crop', true, @islogical);
ip.addParamValue('Interp', 'cubic', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.parse(vol, angle, dz, xyPixelSize, varargin{:});

[ny,nx,nz] = size(vol);

theta = ip.Results.angle * pi/180;
dx = cos(theta)*dz/xyPixelSize; % pixels shifted slice to slice in x

if ~ip.Results.reverse
    xshift = -dx;
    xstep = dx;
else
    xshift = dx + ceil((nz-1)*dx);
    xstep = -dx;
end
nxOut = ceil((nz-1)*dx) + nx; % width of output volume
if ip.Results.Crop
    if nxOut>=2*nx
        cropWidth = floor(nx/2);
    else
        cropWidth = floor(ceil((nz-1)*dx)/2);
    end
    xshift = xshift - cropWidth;
    nxOut = nxOut - cropWidth*2;
end

% shear transform matrix
S = [1 0 0 0;
    0 1 0 0;
    xstep 0 1 0;
    xshift 0 0 1];
tform = affine3d(S);

RA = imref3d([ny nxOut nz], 1, 1, 1);
volout = imwarp(vol, tform, ip.Results.Interp, 'FillValues', 0, 'OutputView', RA);
