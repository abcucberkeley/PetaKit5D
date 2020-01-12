function sph = binarySphere(radius)
%BINARYSPHERE creates a 3D spherical neighborhood/structuring element for morphological operations
% 
% sphereMat = binarySphere(radius)
% 
% This function genereates a 3D logical matrix with values inside a sphere
% of the specified radius being true and those outside being false. This
% sphere can be used to perform morphological operations such is dilation
% and erosion on 3D binary matrices.
% 
% Input:
% 
%   radius - Positive scalar. The radius of the sphere to generate.
%           Note that for: 
%              0       <  radius < 1       a single true voxel is returned
%              1       <= radius < sqrt(2) the standard 6-connected neighborhood is returned.
%              sqrt(2) <= radius < sqrt(3) the standard 18-connected neighborhood is returned
%              sqrt(3) <= radius < 2       the standard 26-connected neighborhood is returned
%               
% 
% Output:
% 
%   sphereMat - The resulting 3D cubic logical matrix (neighborhood), of size
%   ~2*radius+1 (this is only exact if radius is an integer)
% 
% Hunter Elliott
% 2/2010
%

if nargin < 1 || isempty(radius) || numel(radius)>1 || radius <= 0
    error('You must specify a single positive radius!')
end

w = floor(radius);

if round(radius) ~= radius
    %This avoids numerical error in the inequality below, and gives correct
    %output for rad = sqrt(2) and rad = sqrt(3) etc.
    radius = radius + eps(radius);
end

    
%Get x,y,z coordinate matrices for distance-from-origin calculation
[xx,yy,zz] = meshgrid(-w:w,-w:w,-w:w);

%Return all points which are less than radius away from origin of
%neighborhood
sph = xx .^2 + yy .^2 + zz.^2 <= radius^2;
