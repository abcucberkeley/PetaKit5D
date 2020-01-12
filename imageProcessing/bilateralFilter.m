%bf = bilateralFilter(img, sigmaS, sigmaR) filters an image with a bilateral filter.
% The implementation is based on a trigonometric expansion of the Gaussian range kernel.
%
% Inputs:
%             img : input image
%          sigmaS : spatial sigma (Gaussian blur)
%          sigmaR : range (intensity) sigma
%
% Based on:
% [1] Chaudhury et al., IEEE Trans. Imag. Proc. 2011
%
% Francois Aguet, 07/12/2011

function bf = bilateralFilter(img, sigmaS, sigmaR)

T = max(img(:));
[ny,nx] = size(img);

gamma = pi/(2*T);
rho = gamma*sigmaR;

% Lookup table for optimal number of coefficients (see ref. [1])
sigmaRThreshold = [200 150 100 80 60 50 40]/255*T;
N0 = [1 2 3 4 5 7 9];

if sigmaR >= sigmaRThreshold(end)
    N = N0(find(sigmaR>=sigmaRThreshold, 1, 'first'));
elseif sigmaR > 1/gamma^2
    N = 50; % arbitrary high N
else
    N = ceil(1/(gamma*sigmaR)^2);
end

% determine cutoff for insignificant coefficients
if N > 20
    bounds = norminv([0.025 0.975], N/2, sqrt(N)/2);
    bounds = [floor(bounds(1)) ceil(bounds(2))];
else
    bounds = [0 N];
end

% binomial coefficients/weights
c = binopdf(0:N,N,0.5);

h = zeros(ny,nx);
g = zeros(ny,nx);
for n = bounds(1):bounds(2);
    tmp = 1i*gamma*(2*n-N)*img/(rho*sqrt(N));
    hx = exp(tmp);
    d = c(n+1)*exp(-tmp);
    h = h + d.*filterGauss2D(hx, sigmaS);
    g = g + d.*filterGauss2D(img.*hx, sigmaS);
end
bf = real(g./h);
