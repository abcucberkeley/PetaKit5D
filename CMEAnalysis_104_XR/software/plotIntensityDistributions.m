%plotIntensityDistributions(data, varargin) displays the average intensity of 
% the first 10 frames of all trajectories as a function of lifetime cohort
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

% Francois Aguet, 2012

function plotIntensityDistributions(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('CohortBounds', [0 10 15 20 40 60 120]);
ip.addParamValue('ShowPct', false, @islogical);
ip.addParamValue('FigureName', '');
ip.addParamValue('LifetimeData', 'lifetimeData.mat');
ip.parse(varargin{:});
lb = ip.Results.CohortBounds(1:end-1)+1;
ub = ip.Results.CohortBounds(2:end);

lftData = getLifetimeData(data, 'LifetimeData', ip.Results.LifetimeData,...
    'Cutoff_f', 5, 'ExcludeVisitors', false, 'Scale', true);

A = vertcat(lftData.A);
lft = vertcat(lftData.lifetime_s);

ivec = linspace(0, prctile(max(A,[],2), 95), 40); % intensity range
tvec = (0:10)*data(1).framerate; % time vector
di = ivec(2)-ivec(1);

fset = loadFigureSettings('print');

if ~ip.Results.ShowPct
    ha = setupFigure(2,3, 'AxesWidth', 1.8, 'AxesHeight', 1.6, 'Box', 'on',...
        'XSpace', [1.5 0.3 0.5], 'YSpace', [1.5 0.3 0.5], 'SameAxes', true,...
        'Name', ip.Results.FigureName);
    
else
    ha = setupFigure(2,3, 'AxesWidth', 1.8, 'AxesHeight', 1.6, 'Box', 'on',...
        'XSpace', [1.5 0.3 0.5], 'YSpace', [2.5 0.3 0.5], 'SameAxes', true,...
        'Name', ip.Results.FigureName);
end
colormap(jet(256));

nc = numel(lb);
cv = jet(nc);
M = cell(1,nc);
for c = 1:6
    M{c} = A(lb(c)<=lft&lft<=ub(c),1:numel(tvec));
    T = repmat(tvec, [size(M{c},1),1]);
    mv = M{c}(:);
    tv = T(:);
    rmIdx = mv>ivec(end);
    mv(rmIdx) = [];
    tv(rmIdx) = [];
    hm = hist3([mv tv], {ivec, tvec});
    hm = hm./repmat(sum(hm,1), [numel(ivec) 1]);
    imagesc(tvec, ivec, hm, 'Parent', ha(c));
    if ip.Results.ShowPct
        stairsXT(tvec, prctile(M{c},50,1), 'bounds', 'open', 'EdgeColor', cv(c,:), 'Parent', ha(c));
    end
    text(tvec(end), 0.975*ivec(end), ['[' num2str(lb(c)) '...' num2str(ub(c)) '] s'],...
        'Color', 'w', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',...
        fset.tfont{:}, 'FontWeight', 'bold', 'Parent', ha(c));
    
    if c==1
        text(tvec(end), 1.2*ivec(end), 'Lifetime cohort', 'Color', 'k',...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', fset.tfont{:}, 'Parent', ha(c))
    end
end
set(ha, 'XLim', [tvec(1)-data(1).framerate/2 tvec(end)+data(1).framerate/2],...
        'YLim', [ivec(1)-di/2 ivec(end)+di/2], 'FontSize', 8);

yp = ylabel(ha(4), 'Fluo. intensity (A.U.)', fset.lfont{:});
ypos = get(yp, 'Position');
ypos(2) = ivec(end);
set(yp, 'Position', ypos);
xlabel(ha(5), 'Time (s)', fset.lfont{:});

if ip.Results.ShowPct
    pos = get(ha(1), 'Position');
    pos(2) = 0.1;
    pos(4) = pos(4)/2; 
    axes('Position', pos, 'XLim', [tvec(1)-data(1).framerate/2 tvec(end)+data(1).framerate/2],...
        'FontSize', 8, 'TickDir', 'out', 'TickLength', get(ha(1), 'TickLength'));
    for c = 1:nc
            stairsXT(tvec, prctile(M{c},50,1), 'bounds', 'open', 'EdgeColor', cv(c,:));
    end
end
