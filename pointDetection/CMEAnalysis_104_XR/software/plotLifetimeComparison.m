%plotLifetimeComparison(lftRes, varargin) plots CCP lifetime distributions from different conditions
%
% Input:
%         lftRes : cell array of structures returned by runLifetimeAnalysis()
%         legend : string array of the same size as 'lftRes'
%
% Options ('specifier', value):
%      'Frequency' : 'relative' normalizes the distributions by proportion of CCPs
%        'PlotAll' : true|{false} also displays the lifetime distributions for all CCSs
%          'Color' : color matrix, Nx3
%
% Note: the control condition should be first in the 'lftRes' array
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

% Francois Aguet, 06/03/2013

function [ha] = plotLifetimeComparison(lftRes, varargin)

N = numel(lftRes);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('lftRes', @iscell);
ip.addOptional('legend', arrayfun(@(i) [' Condition ' num2str(i)], 1:N, 'unif', 0));
ip.addParamValue('PlotOrder', 1:N);
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('Frequency', '', @(x) strcmpi(x, 'relative'));
ip.addParamValue('Color', []);
ip.addParamValue('SlaveName', [], @iscell);
ip.addParamValue('Control', []);
ip.addParamValue('XSpace', []);
ip.addParamValue('ShowUncertainty', false, @islogical);
ip.addParamValue('forceSingle', false, @islogical);
ip.addParamValue('DisplayMode', 'screen', @(x) any(strcmpi(x, {'print', 'screen'})));
ip.parse(lftRes, varargin{:});

fset = loadFigureSettings(ip.Results.DisplayMode);

% normalization
if strcmpi(ip.Results.Frequency, 'relative')
    w = cellfun(@(i) mean(i.pctCCP), lftRes);
    w = w/w(1);
else
    w = ones(N,1);
end
hp = zeros(N,1);

cv = ip.Results.Color;
if isempty(cv)
    cv = hsv2rgb([0 0 0;
        0.3 1 1;
        0.55 1 1;
        0.99 1 1;
        0.11 1 1]);
    if N>5
        cv = [0 0 0; jet(N-1)];
    end
end

cf = rgb2hsv(cv);
idx = cf(:,3)~=0;
cf(idx,2) = 0.4;
cf(~idx,3) = 0.6;
cf = hsv2rgb(cf);

madFactor = 1/norminv(0.75, 0, 1); % MAD -> S.D.

% 1) Plot CCP distributions
if (~all(cellfun(@(i) isfield(i, 'lftHistSlaveCCP'), lftRes))||(ip.Results.forceSingle))
    ha = setupFigure('DisplayMode', ip.Results.DisplayMode);
    for i = ip.Results.PlotOrder
        if ip.Results.ShowUncertainty
            t = lftRes{i}.t;
            mu = w(i)*mean(lftRes{i}.lftHistCCP,1);
            sd = w(i)*madFactor*mad(lftRes{i}.lftHistCCP,1,1);
            fill([t t(end:-1:1)], [mu+sd mu(end:-1:1)-sd(end:-1:1)], cf(i,:),...
                'EdgeColor', 'none');
        end
        hp(i) = plot(lftRes{i}.t, w(i)*lftRes{i}.meanLftHistCCP, '-', 'LineWidth', 1, 'Color', cv(i,:));
    end
    
    hl = legend(ha(end), hp, cellfun(@(i) [' ' i], ip.Results.legend, 'unif', 0), 'EdgeColor', 'w');
    if strcmpi(ip.Results.DisplayMode, 'print')
        set(hl, fset.sfont{:}, 'Units', 'centimeters',...
            'Position', [3.5 5-(N)*0.4 2 N*0.4]);
    end
else
    na = 2;
    % adapted from plotIntensityCohorts
    pctS = cellfun(@(i) mean(i.pctSlaveCCP,1), lftRes, 'unif', 0);
    pctS = vertcat(pctS{:});
    meanPct = mean(pctS,1);
    stdPct = std(pctS,[],1);
    
    SlaveName = ip.Results.SlaveName;
    if isempty(SlaveName)
        SlaveName = cell(1,na);
    end
    sigCombIdx = [1 0]';
    tmp = sigCombIdx;
    tmp(tmp==1) = '+';
    tmp(tmp==0) = '-';
    atext = cell(1,na);
    switch na
        case 2
            for a = 1:na
                atext{a} = [tmp(a,1) SlaveName{1} ': ' num2str(meanPct(a)*100, '%.1f') '�' num2str(stdPct(a)*100, '%.1f') '%'];
            end
        case 3
            for a = 1:na
                atext{a} = [SlaveName{1} tmp(a,1) ' / ' SlaveName{2} tmp(a,2) ': '...
                    num2str(meanPct(a)*100, '%.1f') '�' num2str(stdPct(a)*100, '%.1f') '%'];
            end
    end
   
    ha = setupFigure(1,2, 'DisplayMode', ip.Results.DisplayMode, 'SameAxes', true, 'YSpace', [1.5 1 0.75]);
    legendText = cell(1,N);
    if ~isempty(ip.Results.Control)
        plot(ha(1), ip.Results.Control.t, ip.Results.Control.meanLftHistCCP, 'k', 'LineWidth', 1);
    end
    for i = ip.Results.PlotOrder
        if ip.Results.ShowUncertainty
            t = lftRes{i}.t;
            mu = w(i)*mean(lftRes{i}.lftHistSlaveCCP{1},1);
            sd = w(i)*madFactor*mad(lftRes{i}.lftHistSlaveCCP{1},1,1);
            fill([t t(end:-1:1)], [mu+sd mu(end:-1:1)-sd(end:-1:1)], cf(i,:),...
                'EdgeColor', 'none', 'Parent', ha(1));
        end
        hp(i) = plot(ha(1), lftRes{i}.t, w(i)*mean(lftRes{i}.lftHistSlaveCCP{1},1), '-', 'LineWidth', 1, 'Color', cv(i,:));
        legendText{i} = [' +' SlaveName{1} ', ' ip.Results.legend{i} ', ' num2str(100*pctS(i,1), '%.1f') '%'];
    end
    hl = legend(ha(1), hp, legendText, 'EdgeColor', 'w');
    if strcmpi(ip.Results.DisplayMode, 'print')
        set(hl, fset.sfont{:}, 'Units', 'centimeters',...
            'Position', [2.5 5-(N-1)*0.4 2 N*0.4]);
    end
    
    if ~isempty(ip.Results.Control)
        plot(ha(2), ip.Results.Control.t, ip.Results.Control.meanLftHistCCP, 'k', 'LineWidth', 1);
    end
    for i = ip.Results.PlotOrder
        if ip.Results.ShowUncertainty
            t = lftRes{i}.t;
            mu = w(i)*mean(lftRes{i}.lftHistSlaveCCP{2},1);
            sd = w(i)*madFactor*mad(lftRes{i}.lftHistSlaveCCP{2},1,1);
            fill([t t(end:-1:1)], [mu+sd mu(end:-1:1)-sd(end:-1:1)], cf(i,:),...
                'EdgeColor', 'none', 'Parent', ha(2));
        end
        hp(i) = plot(ha(2), lftRes{i}.t, w(i)*mean(lftRes{i}.lftHistSlaveCCP{2},1), '-', 'LineWidth', 1, 'Color', cv(i,:));
        legendText{i} = [' -' SlaveName{1} ', ' ip.Results.legend{i} ', ' num2str(100*pctS(i,2), '%.1f') '%'];
    end
    hl = legend(ha(2), hp, legendText, 'EdgeColor', 'w');
    if strcmpi(ip.Results.DisplayMode, 'print')
        set(hl, fset.sfont{:}, 'Units', 'centimeters',...
            'Position', [9.5 5-(N-1)*0.4 2 N*0.4], 'EdgeColor', 'w');
    end
end

XLim = [0 160];
YLim = [0 0.03];
set(ha, 'XLim', XLim, 'YLim', YLim, 'XTick', 0:20:200, 'YTick', 0:0.01:0.1);
arrayfun(@(i) xlabel(i, 'Lifetime (s)', fset.lfont{:}), ha);
if strcmpi(ip.Results.Frequency, 'relative')
    ylabel(ha(1), 'Relative frequency', fset.lfont{:});
else
    ylabel(ha(1), 'Frequency', fset.lfont{:});
end


if ip.Results.ShowUncertainty
    dt = lftRes{1}.t(2) - lftRes{1}.t(1);
    comb = pcombs(1:N);
    
    fprintf('Comparison of CCP lifetime distributions:\n');
    for i = 1:size(comb,1)
        % since histograms are already corrected for observation bias
        % calculate EDFs from histograms
        d1 = cumsum(lftRes{comb(i,1)}.lftHistCCP,2) * dt;
        d2 = cumsum(lftRes{comb(i,2)}.lftHistCCP,2) * dt;
        [hval, pval] = permKSTestMeanEDF(d1, d2, 'CmpFunction', @nanmean);
        fprintf('Cond. %d vs Cond. %d: H=%d, p = %.4f\n', comb(i,1), comb(i,2), hval, pval);
    end
    
    fprintf('Comparison of CS lifetime distributions:\n');
    for i = 1:size(comb,1)
        % since histograms are already corrected for observation bias
        % calculate EDFs from histograms
        d1 = cumsum(lftRes{comb(i,1)}.lftHistCS,2) * dt;
        d2 = cumsum(lftRes{comb(i,2)}.lftHistCS,2) * dt;
        [hval, pval] = permKSTestMeanEDF(d1, d2, 'CmpFunction', @nanmean);
        fprintf('Cond. %d vs Cond. %d: H=%d, p = %.4f\n', comb(i,1), comb(i,2), hval, pval);
    end
    
end


% 2) Optional: plot CCS distributions (all objects)
if ip.Results.PlotAll
    setupFigure('DisplayMode', ip.Results.DisplayMode);
    hold on;
    for i = ip.Results.PlotOrder
        hp(i) = plot(lftRes{i}.t, nanmean(lftRes{i}.lftHist_Ia,1), '-', 'LineWidth', 1, 'Color', cv(i,:));
    end
    hl = legend(hp, ip.Results.legend);
    set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [4 5-N*0.3 2 N*0.3]);
    axis([0 120 0 0.05]);
    set(gca, 'YTick', 0:0.01:0.05);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
end
