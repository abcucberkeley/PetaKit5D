function [out, G] = filterGauss2D(image, sigma, borderCondition)
% filterGauss2D :	filters an image with a 2-D Gaussian mask
%
%    [out, G] = filterGauss2D(image, sigma, borderCondition);
%
%    INPUT: image           : 2-D input array
%           sigma           : standard deviation of the Gaussian
%           borderCondition : input for 'padarrayXT'. Default: 'symmetric'
%                             Options: 'symmetric', 'replicate', 'circular', 'antisymmetric', or a constant value
%
%    OUTPUT: out : filtered image
%            G   : Gaussian mask
%
% Francois Aguet, added 01/21/2010

if nargin < 3 || isempty(borderCondition)
    borderCondition = 'symmetric';
end

w = ceil(3*sigma); % cutoff radius of the gaussian kernel
x = -w:w;
g = exp(-x.^2/(2*sigma^2));
g = g/sum(g);

out = conv2(g', g, padarrayXT(image, [w w], borderCondition), 'valid');

if (nargout>1)
    G = g'*g;
end