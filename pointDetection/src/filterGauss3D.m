function [out] = filterGauss3D(input, sigma, borderCondition)
% filterGauss3D :	filters a data volume with a 3-D Gaussian mask
%
%     out = filterGauss3D(image, sigma, borderCondition);
%
%    INPUT: image           : 3-D input array
%           sigma           : standard deviation of the Gaussian
%           borderCondition : input for 'padarrayXT'. Default: 'symmetric'
%                             Options: 'symmetric', 'replicate', 'circular', 'antisymmetric', or a constant value
%
%    OUTPUT: out : filtered volume
%
% Francois Aguet, added 01/21/2010

if nargin < 3 || isempty(borderCondition)
    borderCondition = 'symmetric';
end

w = ceil(3*sigma); % cutoff radius of the gaussian kernel
g = exp(-(-w(1):w(1)).^2/(2*sigma(1)^2));
g = g/sum(g);

if numel(sigma)>1
    gz = exp(-(-w(2):w(2)).^2/(2*sigma(2)^2));
    gz = reshape(gz/sum(gz), [1 1 2*w(2)+1]);
else
    gz = reshape(g, [1 1 2*w+1]);
    w = [w w];
end

out = convn(padarrayXT(input, [w(1) w(1) w(2)], borderCondition), g', 'valid');
out = convn(out, g, 'valid');
out = convn(out, gz, 'valid');