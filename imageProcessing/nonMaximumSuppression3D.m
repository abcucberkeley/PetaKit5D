function nms = nonMaximumSuppression3D(u,v,w)
%NONMAXIMUMSUPPRESSION performs non-maximum suppression (NMS) on the input 3D vector image
% 
% nms = nonMaximumSuppression3D(u,v,w)
% 
% This function performs non-maximum suppression (NMS) on the input 3D
% vector-valued data. That is, pixels are found where the magnitude of the
% value is a local maximum in the direction of that pixel's vector.
% 
% Input:
% 
%   u,v,w - The 3D matrices specifying the X, Y and Z components of the
%           data to perform NMS on. These are assumed to correspond to
%           dimensions 2,1,3 of the input matrix respectively.
%           
% Output: 
% 
%   nms - 3D matrix with zeros in non-maximum areas and the magnitude
%         of the input data at the maxima.
%
% Hunter Elliott
% 9/2011
%

if nargin < 3 || ndims(u) ~=3 || ndims(v) ~= 3 || ndims(w) ~= 3 || ~isequal(size(u),size(v),size(w))
    error('The inputs u, v and w must all be 3-dimensional matrices of equal size!')
end

[M,N,P] = size(u);

%Calculate the magnitude of the vector field at each point
nms = sqrt(u .^2 + v .^2 + w .^2);
%Use this to normalize the vector data
u = u ./ nms;
v = v ./ nms;
w = w ./ nms;

u(isnan(u)) = 0;
v(isnan(v)) = 0;
w(isnan(w)) = 0;

%Pad this array so we can detect maxima at the edges
nmsPad = padarrayXT(nms,[1 1 1],'symmetric');
[X,Y,Z] = meshgrid(1:N,1:M,1:P);%Coordinates for interpolation of unpadded matrices
[Xp,Yp,Zp] = meshgrid(0:N+1,0:M+1,0:P+1);%Coordinates for interpolation of padded matrices

%Get matrices containing the interpolated values 1 pixel in the direction
%parallel and antiparallel to the local orientation
mag1 = interp3(Xp,Yp,Zp,nmsPad,X+u,Y+v,Z+w);%In direction of vector
mag2 = interp3(Xp,Yp,Zp,nmsPad,X-u,Y-v,Z-w);%In opposite direction of vector

%Remove pixels which are not greater than both these values
nms(nms<=mag1 | nms <= mag2) = 0;
