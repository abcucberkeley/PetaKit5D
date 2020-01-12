% [A_est, c_est, pval] = filterGaussianFit2D(img, sigma)
% Performs a 2D Gaussian fit at every pixel of the input image, returning the
% amplitude and background intensity.
%
% Inputs :   
%           img : input image
%         sigma : standard deviation of the Gaussian PSF
%
% Outputs: 
%         A_est : Gaussian amplitude at each pixel
%         c_est : Background intensity at each pixel
%          pval : p-value for t-test between 'A' and background noise.
%                 This indicates the significance of the amplitude at each pixel.
%
% See also pointSourceDetection
%
% Copyright (C) 2017, Danuser Lab - UTSouthwestern 
%
% This file is part of CMEAnalysis_Package.
% 
% CMEAnalysis_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CMEAnalysis_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with CMEAnalysis_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Francois Aguet, 2013

function [A_est, c_est, varargout] = filterGaussianFit2D(img, sigma, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img', @isnumeric);
ip.addRequired('sigma', @isscalar);
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.parse(img, sigma, varargin{:});
alpha = ip.Results.Alpha;

w = ceil(4*sigma);
g = exp(-(-w:w).^2/(2*sigma^2));
u = ones(1,numel(g));
    
g2 = g'*g;
n = numel(g2);
gsum = sum(g2(:));
g2sum = sum(g2(:).^2);

J = [g2(:) ones(n,1)]; % g_dA g_dc
C = inv(J'*J);

% calculate fitted amplitude & background at each pixel
imgXT = double(padarrayXT(img, [w w], 'symmetric'));
fg = conv2(g', g, imgXT, 'valid');
fu = conv2(u', u, imgXT, 'valid');

A_est = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
c_est = (fu - A_est*gsum)/n;

if nargout>2
    fu2 = conv2(u', u, imgXT.^2, 'valid');

    f_c = fu2 - 2*c_est.*fu + n*c_est.^2; % f-c
    RSS = A_est.^2*g2sum - 2*A_est.*(fg - c_est*gsum) + f_c;
    RSS(RSS<0) = 0; % negative numbers may result from machine epsilon/roundoff precision
    sigma_e2 = RSS/(n-3);
    
    sigma_A = sqrt(sigma_e2*C(1,1));
    
    % standard deviation of residuals
    sigma_res = sqrt(RSS/(n-1));
    
    kLevel = norminv(1-alpha/2.0, 0, 1);
    
    SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
    df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
    scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
    T = (A_est - sigma_res*kLevel) ./ scomb;

    varargout{1} = tcdf(-T, df2);
end
