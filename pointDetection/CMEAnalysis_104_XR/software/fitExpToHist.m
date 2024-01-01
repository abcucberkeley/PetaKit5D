function [mu, mu_std, A, pdf] = fitExpToHist(ti, ni)


opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

mu0 = nansum(ni.*ti)/nansum(ni);

[p,resnorm,~,~,~,~,J] = lsqnonlin(@costFct, [mu0 1], [], [], opts, ti, ni);
mu = p(1);
A = p(2);

J = full(J);
mu_std = sqrt( resnorm/(numel(ni)-2) * inv(J'*J) ); %#ok<MINV>

pdf = A * 1/mu * exp(-1/mu*ti);
pdf(isnan(ni)) = NaN;


function v = costFct(p, ti, ni)
mu = p(1);
A = p(2);

pdf = A * 1/mu * exp(-1/mu*ti);
v = pdf-ni;
v(isnan(v)) = [];
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
