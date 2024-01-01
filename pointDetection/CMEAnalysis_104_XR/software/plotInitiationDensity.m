%plotInitiationDensity(lftRes, xlabels, cv, cf, varargin) compares CS vs. CCP initiation density under different conditions
%
% Inputs:
%
%   lftRes : output structure from runLifetimeAnalysis()
%  xlabels : legend; cell array of strings
%       cv : Nx3 color matrix for edges
%       cv : Nx3 color matrix
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

% Francois Aguet

function ha = plotInitiationDensity(lftRes, xlabels, varargin)

nd = numel(lftRes);
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('lftRes', @iscell);
ip.addRequired('xlabels', @iscell);
ip.addOptional('cv', hsv2rgb([(0:1/nd:1-1/nd)' ones(nd,2)]), @(x) size(x,1)==nd & size(x,2)==3);
ip.addOptional('cf', hsv2rgb([(0:1/nd:1-1/nd)' 0.4*ones(nd,1) ones(nd,1)]), @(x) size(x,1)==nd & size(x,2)==3);
ip.addOptional('SigLink', cell(1,2), @iscell)
ip.addParamValue('YLim', []);
ip.addParamValue('YTick', []);
ip.parse(lftRes, xlabels, varargin{:});
cv = ip.Results.cv;
cf = ip.Results.cf;

fset = loadFigureSettings('print');

M_Ia = cellfun(@(i) i.initDensityIa(:,1), lftRes, 'unif', 0);
M_above = cellfun(@(i) i.initDensityCCP(:,1), lftRes, 'unif', 0);

ha = setupFigure(1,2, 'AxesWidth', 2.25, 'AxesHeight', 3.5,...
    'XSpace', [1.5 0.85 0.5], 'YSpace', [2 0.5 0.5]);
set(ha, 'FontSize', 8);


opts = {'LineWidth', 1, 'GroupDistance', 0, 'BarWidth', 0.6, 'XTickLabel', xlabels,...
    'AdjustFigure', false, 'BorderWidth', 0.5, 'Angle', 0, 'DetectOutliers', false,...
    'FaceColor', cf, 'EdgeColor', cv, 'ErrorbarColor', cv};

boxplot2(M_Ia, ip.Results.SigLink{1}, opts{:}, 'Parent', ha(1));
boxplot2(M_above, ip.Results.SigLink{2}, opts{:}, 'Parent', ha(2));

ylabel(ha(1), ['Initiations (' char(181) 'm^{-2} min^{-1})'], fset.lfont{:});

YLim = ip.Results.YLim;
if isempty(YLim)
    YLim = arrayfun(@(i) get(i, 'YLim'), ha, 'unif', 0);
    YLim = cellfun(@(i) [0 i(2)], YLim, 'unif', 0);
    set(ha(1), 'YLim', YLim{1});
    set(ha(2), 'YLim', YLim{2});
end

YTick = ip.Results.YTick;
if ~isempty(YTick)
    set(ha(1), 'YTick', YTick{1});
    set(ha(2), 'YTick', YTick{2});
end

rotateXTickLabels(ha(1), 'Angle', 45, 'AdjustFigure', false);
rotateXTickLabels(ha(2), 'Angle', 45, 'AdjustFigure', false);

formatTickLabels(ha);
