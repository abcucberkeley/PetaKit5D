% ADTEST(X) implements the Anderson-Darling test for normality
%
% Inputs:
%          x : Sample vector (test only work for n>=5)
%    {alpha} : Significance level. Possible values:
%              0.25 0.15 0.10 0.05 [default] 0.025 0.01 0.005 0.0025
%
% Output:
%          H : Result of the hypothesis test
%              0: Do not reject H0 at given significance level
%              1: Reject H0 at given significance level
%       pval : p-value. Only returned for case 3 (mean and variance unknown)   
%         A2 : (Adjusted) test statistic
%       cval : Critical value for the test
%
% NOTE: As of Matlab 2013a, adtest() is a built-in function. The implementation uses
%       different interpolations and is slower, therefore this function remains available.
%       The function prototypes are identical.
%
%
% For the test and its derivation, see
% [1] Anderson & Darling, Ann. Math. Stat. 23, 1952
% [2] Anderson & Darling, J. Am. Stat. Assoc. 49, 1954
%
% Critical values taken from 
% [3] R.B. D'Agostino and M.A. Stephens, Goodness-of-Fit Techniques, ed. M. Dekker, Inc., 1986
%
% Older, less accurate critical values were given in
% [4] M.A. Stephens, J. Am. Stat. Assoc. 69(347), pp. 730-737, 1974
% [5] M.A. Stephens, Ann. Stat. 4(2), pp. 357-369, 1976
%
% See also
% [6] http://en.wikipedia.org/wiki/Anderson-Darling_test
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

% Francois Aguet (last modified 05/23/2012)

function [H, pval, A2, cval] = adtest1(x, varargin)

alphaVec = [0.5 0.25  0.15  0.10  0.05  0.025 0.01  0.005 0.0025];

ip = inputParser;
ip.addRequired('x');
ip.addOptional('alpha', 0.05, @(x) ismember(x, alphaVec));
ip.addParamValue('mu', []);
ip.addParamValue('sigma', []);
ip.addParamValue('Distribution', 'normal', @(x) any(strcmpi(x, {'normal', 'exponential'})));
ip.parse(x, varargin{:});

n = numel(x);
if n<5
    error('The Anderson-Darling test requires at least 5 samples.');
end
x = reshape(x, [1 n]);

estidx = [0 0]; % [sigma mu]

mu = ip.Results.mu;
if isempty(mu)
    mu = mean(x);
    estidx(2) = 1;
end

sigma = ip.Results.sigma;
if isempty(sigma)
    sigma = std(x);
    estidx(1) = 1;
end


% sort samples in ascending order
x = sort(x);

% case
switch ip.Results.Distribution
    case 'normal'
        c = bin2dec(num2str(estidx))+1;
        z = normcdf(x, mu, sigma);

    case 'exponential'
        c = 5;
        z = expcdf(x, mu);
end

% Look-up table for critical values
%  Case 0: mean and variance known
%  Case 1: mean unknown, variance known
%  Case 2: mean known, variance unknown
%  Case 3: mean and variance unknown
%  Case 4: exponential distribution, scale unknown
ctable = zeros(5,9);
% alpha:      [0.5 0.25  0.15  0.10  0.05  0.025 0.01  0.005 0.0025]
%-----------------------------
% Normal distribution
%-----------------------------
% Case 0: Table 4.2 from, [3], p. 105, valid for n>=5
ctable(1,:) = [NaN 1.248 1.610 1.933 2.492 3.070 3.857 4.500 NaN]; % 6.000 for alpha = 0.001

% Cases 1 and 2: Table 4.6 from [3], p. 122
ctable(2,:) = [NaN 0.644 0.782 0.894 1.087 1.285 1.551 1.756 1.964];
ctable(3,:) = [NaN 1.072 1.430 1.743 2.308 2.898 3.702 4.324 4.954];

% Case 3: Table 4.7 from [3], p. 123
ctable(4,:) = [0.341 0.470 0.561 0.631 0.752 0.873 1.035 1.159 NaN];

%-----------------------------
% Exponential distribution
%-----------------------------
% Case 4: Table 4.11 from [3]
ctable(5,:) = [NaN 0.736 0.916 1.062 1.321 1.591 1.959 2.244 2.534];


% Select critical value
cval = ctable(c,alphaVec==ip.Results.alpha);
      
% Test statistic
i = 1:n;
A2 = -n - 1/n*sum( (2*i-1).*(log(z) + log(1-z(n+1-i))) );


% Note: A2 is not modified for cases 0 (Table 4.2) and 1 & 2 (Table 4.6)
pval = NaN;
if c==4
    % Modified statistic for case 3
    A2 = A2 * (1.0 + 0.75/n + 2.25/n^2); % [1] Table 4.7, replaces the less accurate correction factor (1+4/n-25/n^2)
    
    % p-values for case 3 ([1] p.127, Table 4.9)
    if (0.600<A2)
        pval = exp(1.2937 - 5.709*A2 + 0.0186*A2*A2);
    end
    if (0.340<A2 && A2<=0.600)
        pval = exp(0.9177 - 4.279*A2 - 1.38*A2*A2);
    end
    if (0.200<A2 && A2<=0.340)
        pval = 1 - exp(-8.318 + 42.796*A2 - 59.938*A2*A2);
    end
    if (A2<=0.200)
        pval = 1 - exp(-13.436 + 101.14*A2 - 223.73*A2*A2);
    end
end
if c==5
    A2 = A2 * (1.0 + 0.6/n);
end
    
H = A2 > cval;
