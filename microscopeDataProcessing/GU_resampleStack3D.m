function ivol = GU_resampleStack3D(vol, xf, yf, zf, varargin)
% interpolates the volumes to resample based on the x,y,z voxel scaling
% factor
% Gokul Upadhyayula, Nov 2017

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol'); % typical value: 32.8
ip.addRequired('xf'); % typical value: 32.8
ip.addRequired('yf'); % typical value: 32.8
ip.addRequired('zf'); % typical value: 32.8
ip.addParameter('Interp', 'nearest', @(x) any(strcmpi(x, {'cubic', 'linear', 'nearest'})));
ip.addParameter('Save16bit', false, @islogical);
ip.parse(vol, xf, yf, zf, varargin{:});
p = ip.Results;

% interpolate
[ny,nx,nz] = size(vol);
[y,x,z] = ndgrid(1:ny,1:nx,1:nz);
[Y,X,Z] = ndgrid(1:yf:ny,1:xf:nx,1:zf:nz);
ivol = interp3(x,y,z,double(vol),X,Y,Z,p.Interp);

if p.Save16bit
    ivol = uint16(ivol);
else
    ivol = single(ivol);
end