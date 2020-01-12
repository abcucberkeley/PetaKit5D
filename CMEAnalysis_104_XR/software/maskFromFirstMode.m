% Francois Aguet, 02/17/2012
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

function [mask, T] = maskFromFirstMode(img, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addParamValue('Connect', true, @islogical);
ip.addParamValue('Display', false, @islogical);
ip.addParamValue('ModeRatio', 0.6, @isscalar);
ip.parse(img, varargin{:});

[ny,nx] = size(img);

g = filterGauss2D(img, 5);

v = g(:);
v(isnan(v)) = [];

pct = prctile(v, [1 99]);
v(v<pct(1) | pct(2)<v) = [];

[f,xi] = ksdensity(v, 'npoints', 100);

% local max/min
lmax = locmax1d(f, 3);
lmin = locmin1d(f, 3);

dxi = xi(2)-xi(1);

T = [];

% identify min after first mode
if ~isempty(lmin)
    idx = find(lmin>lmax(1), 1, 'first');
    if ~isempty(idx) && sum(f(1:lmin(idx(1))))*dxi < ip.Results.ModeRatio
        min0 = lmin(idx);
        T = xi(min0);
        mask = g>T;
        
        % retain largest connected component
        if ip.Results.Connect
            CC = bwconncomp(mask, 8);
            compsize = cellfun(@(i) numel(i), CC.PixelIdxList);
            [~,idx] = sort(compsize, 'descend');
            mask = zeros([ny,nx]);
            mask(vertcat(CC.PixelIdxList{compsize>=compsize(idx(1))})) = 1;
        end
    else
        mask = ones(ny,nx);
    end
else
    mask = ones(ny,nx);
end

if ip.Results.Display
    dx = xi(2)-xi(1);
    ni = hist(v, xi);
    ni = ni/sum(ni)/dx;
    
    figure;
    hold on;
    plot(xi, ni, 'k.-', 'LineWidth', 3, 'MarkerSize', 20);
    plot(xi, f, 'r-', 'LineWidth', 1);
    set(gca, 'LineWidth', 2, 'FontSize', 18);
end