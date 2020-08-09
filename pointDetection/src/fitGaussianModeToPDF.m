%[mu, sigma, xi, g] = fitGaussianModeToCDF(samples, varargin) fits a Gaussian to the lower half of the first mode of the sample distribution
% The fit is performed on the PDF of the samples, calculated using kernel density.
%
% Inputs:
%         samples : data points
%
% Outputs:
%              mu : mean of the Gaussian
%           sigma : standard deviation of the Gaussian
%              xi : sample space vector
%               g : Gaussian calculated on xi

% Francois Aguet, 05/24/2012

function [mu, sigma, xi, g] = fitGaussianModeToPDF(samples, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Display', false, @islogical);
ip.addParamValue('FixMode', false, @islogical);
ip.parse(varargin{:});

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

pct = prctile(samples, [0 99.9]);
xi = linspace(pct(1), pct(2), 1000);
[pdf,xi] = ksdensity(samples,xi);

A0 = 0.5;
%HE - This seems to give more reliable initialization.
[mu0,sigma0] = robustMean(samples(:),1,2);

if ~ip.Results.FixMode
    p = lsqnonlin(@cost, [A0 mu0 sigma0], [0 -Inf 0], [1 Inf Inf], opts, xi, pdf);
    A = p(1);
    mu = p(2);
    sigma = p(3);
else
    % estimate only the sigma parameter. A and mu are fixed by the mode.
    p = lsqnonlin(@costSigma, sigma0, 0, Inf, opts, xi, pdf);
    sigma = p(1);
    [A,idx] = max(pdf);
    mu = xi(idx);
end
g = A*exp(-(xi-mu).^2/(2*sigma^2));


if ip.Results.Display
    figure;
    hold on;
    
    pct = prctile(samples, [0 99.9]);
    xh = linspace(pct(1), pct(2), 50);
    [nh,xh] = hist(samples, xh);
    nh = nh/sum(nh)/(xh(2)-xh(1));
    bar(xh, nh, 'BarWidth', 1, 'FaceColor', 'none', 'EdgeColor', 0.6*[1 1 1]);
    
    plot(xi, pdf, 'k-');
    plot(xi, g, 'r');
    
    hl = legend(' Histogram', ' Density', ' Gaussian fit to first mode');
    set(hl, 'EdgeColor', 'w');
end



function v = cost(p, xi, pdf)
A = p(1);
mu = p(2);
sigma = p(3);

g = A*exp(-(xi-mu).^2/(2*sigma^2));
v = g - pdf;
T = find(xi>mu+sigma, 1, 'first');
v(T:end) = 0;


function v = costSigma(p, xi, pdf)
sigma = p(1);
[A,idx] = max(pdf(:));
mu = xi(idx);

g = A*exp(-(xi-mu).^2/(2*sigma^2));
v = g - pdf;
v(xi>xi(idx)) = 0;

