%[mu, sigma, hf] = plotTrackClasses(v) generates a bar graph of track categories
%
% Inputs:
%          v: vector of track categories (values in [1...8])
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

% Francois Aguet, 02/01/2012

function [mu, sigma, hf] = plotTrackClasses(v, varargin)

fset = loadFigureSettings('screen');

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('v');
ip.addOptional('c', []);
ip.addParamValue('Handle', []);
ip.addParamValue('YLim', []);
ip.addParamValue('YTick', 0:10:100);
ip.addParamValue('FaceColor', fset.cfTrackClasses);
ip.addParamValue('EdgeColor', fset.ceTrackClasses);
ip.parse(v, varargin{:});

xlabels = {'single tracks', 'single tracks, rejected', 'single tracks, cut', 'single tracks, persistent',...
    'comp. tracks', 'comp. tracks, rejected', 'comp. tracks, cut', 'comp. tracks, persistent'};


% Setup figure window
if ~isempty(ip.Results.Handle)
    ha = ip.Results.Handle;
    hf = get(ha, 'Parent');
else
    [ha,~,hf] = setupFigure('DisplayMode', 'screen');
end

if iscell(v)
    %vpos = cellfun(@(i,j) i(j==1), v, ip.Results.c, 'unif', 0);
    %vneg = cellfun(@(i,j) i(j==0), v, ip.Results.c, 'unif', 0);
    %v = arrayfun(@(i) hist(vpos{i}, 1:8)/numel(v{i}), 1:numel(v), 'unif', 0);
    v = arrayfun(@(i) hist(v{i}, 1:8)/numel(v{i}), 1:numel(v), 'unif', 0);
    v = vertcat(v{:});
    mu = mean(v,1);
    %rstd = @(x) 1/norminv(0.75, 0, 1) * mad(x, 1, 1);
    %sigma = rstd(v);
    sigma = std(v,[],1);
else
    mu = hist(v, 1:8)/numel(v);
    sigma = [];
end
mu = 100*mu;
sigma = 100*sigma;

YLim = ip.Results.YLim;
if isempty(YLim)
    tmp = mu;
    if ~isempty(sigma)
        tmp = tmp+sigma;
    end
    YLim = [0 ceil(max(tmp)/20)*20];
end

barplot2(mu, sigma, 'Handle', ha, 'BarWidth', 0.6, 'GroupDistance', 0.8,...
    'FaceColor', ip.Results.FaceColor, 'EdgeColor', ip.Results.EdgeColor,...
    'XTickLabel', xlabels, 'YTick', 0:20:100, 'YLim', YLim);
ylabel('% tracks', fset.lfont{:});
