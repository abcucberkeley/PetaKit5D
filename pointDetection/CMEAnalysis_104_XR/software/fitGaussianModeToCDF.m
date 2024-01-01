%[mu, sigma, xi, g] = fitGaussianModeToCDF(samples, varargin) fits a Gaussian to the lower half of the first mode of the sample distribution
% The fit is performed on the EDF of the samples.
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

function [mu, sigma, xi, g] = fitGaussianModeToCDF(samples, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Display', false, @islogical);
ip.parse(varargin{:});

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

[f_edf, x_edf] = ecdf(samples);

A0 = 0.5;
mu0 = mean(samples);
sigma0 = 0.5*std(samples);
p = lsqnonlin(@cost, [A0 mu0 sigma0], [0 0 0], [1 Inf Inf], opts, x_edf, f_edf);
mu = p(2);
sigma = p(3);

if nargout>2  || ip.Results.Display
    pct = prctile(samples, [0 99.9]);
    xi = linspace(pct(1), pct(2), 1000);
    [f,xi] = ksdensity(samples,xi);
    g = exp(-(xi-mu).^2/(2*sigma^2)) / sqrt(2*pi)/sigma;
    
    % adjust A for density plot
    T = find(xi>mu, 1, 'first')-1;
    A = sum(f(1:T).*g(1:T)) / sum(g(1:T).^2);
    g = A*g;
end

if ip.Results.Display
    figure;
    hold on;
    plot(xi, f, 'k-');
    plot(xi, g, 'r');
end



function v = cost(p, x_edf, f_edf)
A = p(1);
mu = p(2);
sigma = p(3);

cdf = 0.5 * (1 + erf((x_edf-mu)/sqrt(2)/sigma));
v = (A*cdf - f_edf)/numel(f_edf);
T = find(cdf>0.5, 1, 'first');
v(T:end) = 0;
