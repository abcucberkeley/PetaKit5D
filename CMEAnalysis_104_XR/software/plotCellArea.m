function plotCellArea(lftRes, labels)

fset = loadFigureSettings('print');

mu = cellfun(@(i) mean(i.cellArea), lftRes);
%sigma = cellfun(@(i) std(i.cellArea), lftRes);
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
sigma = cellfun(@(i) std(i.cellArea)/sqrt(numel(i.cellArea)), lftRes);

figure(fset.fOpts{:}, 'Position', [5 5 5.5 7]);
axes(fset.axOpts{:}, 'Position', [2 3 3 3.5], 'TickLength', fset.TickLength/3.5*6);
barplot2(mu', sigma', 'Angle', 0, 'BarWidth', 1, 'GroupDistance', 1,...
    'FaceColor', 0.8*[1 1 1], 'EdgeColor', 0.4*[1 1 1], ... %'AxisFontSize', 8,...
    'LineWidth', 1);
ylabel(['Cell area (' char(181) 'm^2)'], fset.lfont{:})

set(gca, 'XTickLabel', labels);
rotateXTickLabels(gca, 'AdjustFigure', false);
