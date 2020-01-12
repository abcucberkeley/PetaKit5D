function [dX,dY] = gradientFilterGauss2D(input, sigma, borderCondition)
% gradientfilterGauss2D :	gradient filters an image with a 2-D Gaussian gradient mask
%
%    [dX,dY] = gradientFilterGauss2D(image, sigma, borderCondition);
%
%       Filters the input matrix using partial derivatives of a gaussian,
%       giving a filtered gradient image.
%
%    INPUT: image           : 2-D input array
%           sigma           : standard deviation of the Gaussian to use
%                             derivatives of for filtering
%           borderCondition : input for 'padarrayXT'. Default: 'symmetric'
%                             Options: 'symmetric', 'replicate', 'circular', 'antisymmetric', or a constant value
%
%    OUTPUT: [dX,dY]    : Matrices filtered with partial derivatives of the
%                         gaussian in the X and Y directions espectively,
%                         corresponding to matrix dimensions 2, and 1
%                         respectively.
%
% Hunter Elliott
% 2/2014


if nargin < 3 || isempty(borderCondition)
    borderCondition = 'symmetric';
end

w = ceil(3*sigma); % cutoff radius of the gaussian kernel
x = -w:w;
g = exp(-x.^2/(2*sigma^2));
dg = -x / sigma^2 .* exp(-x.^2/(2*sigma^2));
gSum = sum(g);
g = g/gSum;
dg = dg/gSum;

dX = convn(padarrayXT(input, [w w], borderCondition), dg, 'valid');
dX = convn(dX, g', 'valid');

dY = convn(padarrayXT(input, [w w], borderCondition), g, 'valid');
dY = convn(dY, dg', 'valid');