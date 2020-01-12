function W = awt1D(signal, nBands)
% W = AWT(I) computes the A Trou Wavelet Transform of image I.
% A description of the algorithm can be found in:
% J.-L. Starck, F. Murtagh, A. Bijaoui, "Image Processing and Data
% Analysis: The Multiscale Approach", Cambridge Press, Cambridge, 2000.
%
% W = AWT(I, nBands) computes the A Trou Wavelet decomposition of the
% image I up to nBands scale (inclusive). The default value is nBands =
% ceil(max(log2(N), log2(M))), where [N M] = size(I).
%
% Output:
% W contains the wavelet coefficients, an array of size N x M x nBands+1.
% The coefficients are organized as follows:
% W(:, :, 1:nBands) corresponds to the wavelet coefficients (also called
% detail images) at scale k = 1...nBands
% W(:, :, nBands+1) corresponds to the last approximation image A_K.
%
% You can use awtDisplay(W) to display the wavelet coefficients.
%
% Sylvain Berlemont, 2009
% Francois Aguet, 2010

N = length(signal);
K = ceil(log2(N));

if nargin<2
    nBands = K;
else
    if nBands < 1 || nBands > K
        error('invalid range for nBands parameter.');
    end
end
W = zeros(N, nBands+1);

signal = double(signal(:));
lastA = signal;

for k = 1:nBands
    newA = convolve(lastA, k);
    W(:,k) = lastA - newA;
    lastA = newA;
end
W(:,nBands+1) = lastA;


function F = convolve(signal, k)
N = length(signal);
k1 = 2^(k - 1);
k2 = 2*k1;

tmp = padarrayXT(signal, [k2 0], 'symmetric');

for i = k2+1:k2+N
    signal(i-k2) = 6*tmp(i) + 4*(tmp(i+k1) + tmp(i-k1)) + tmp(i+k2) + tmp(i-k2);
end
F = signal/16;
