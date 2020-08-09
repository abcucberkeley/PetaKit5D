%CONV3FAST Fast 3D convolution with symmetric kernels
%
%  Usage:
%    [F] = conv3fast(volume, kernel)
%    [F] = conv3fast(volume, xKernel, yKernel, zKernel)
%
%  Outputs:
%     F : filtered volume
%
%  Example: convolution with a Gaussian kernel
%     s = 2;
%     w = ceil(4*s);
%     g = exp(-(0:w).^2/(2*s^2)); % symmetric kernel starts at '0'
%     F = conv3fast(data, g);
%
%  Note: NaNs in input are allowed

% Francois Aguet, 09/19/2013