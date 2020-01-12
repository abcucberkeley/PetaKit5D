% [c] = b3spline2D(img, boundary)
% 
% Computes the cubic B-spline coefficients.
%
% Inputs: 
%           img      : input image
%           boundary : boundary conditions: 'mirror' (default) or 'periodic'

% Francois Aguet, June 2010

function c = b3spline2D(img, boundary)

if nargin<2
    boundary = 'mirror';
end

c = b3spline1D(img, boundary);
c = b3spline1D(c', boundary)';
