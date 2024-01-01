%BOXPLOT2 Box plot grouping multiple sets/categories of data, with error bars and SEM.
%
% Inputs:   prm : matrix or cell array of matrices that contain the box properties:
%                 row 1: mean or median
%                 row 2: optional, SEM
%                 row 2/3: 25th percentile, bottom of box
%                 row 3/4: 75th percentile, top of box
%                 row 4/5: optional, bottom whisker
%                 row 5/6: optional, top whisker
% Options:
%  
%     FaceColor : Nx3 matrix of colors, where N is the number of bars or groups
%     EdgeColor : "
%       xLabels : cell array of strings, labels for each bar
%        yLabel : string, y-axis label
%
% Examples:
%
% 1) Simple box plot
% prm = [3 4; 0.2 0.2; 2 3; 4 5; 0.5 0.5; 0.5 0.5];
% figure; boxplot2(prm, 'BarWidth', 0.8, 'XLabel', 'x label', 'YLabel', 'y label', ...
%     'XTickLabel', arrayfun(@(k) ['S' num2str(k)], 1:2, 'UniformOutput', false),...
%     'Angle', 0);
% 
% 
% 2) Multiple groups
% prm = {[3 4; 0.2 0.2; 2 3; 4 5; 0.5 0.5; 0.5 0.5],...
%     [3 4; 0.2 0.2; 2 3; 4 5; 0.5 0.5; 0.5 0.5]};
% figure; boxplot2(prm, 'BarWidth', 0.8, 'XLabel', 'x label', 'YLabel', 'y label', ...
%     'XTickLabel', arrayfun(@(k) ['Group label ' num2str(k)], 1:2, 'UniformOutput', false),...
%     'Angle', 45, 'FaceColor', [1 0.5 0.5; 0.5 1 0.5], 'EdgeColor', [0.8 0 0; 0 0.8 0]);
% 
% 
% Francois Aguet, 22 Feb 2011 (Last modified: 06/14/2013)
% Andrew R. Jamieson, Dec 2016 -- setErrorbars deprecated now.
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

function h = boxplot2(prm, varargin)

if isnumeric(prm)
    prm = {prm};
end

% input can either be a cell array of vectors, or a cell array of 5/6xnb matrices
[s1,s2] = cellfun(@(i) size(i), prm);
s1 = unique(s1);
if all(cellfun(@isvector, prm))
    % input is cell array of vectors
    [nb, ng] = size(prm); % #bars in each group, #groups
    rawdata = true;
elseif numel(s1)==1 && any(s1==[5 6]) && numel(unique(s2))==1
    ng = numel(prm);
    nb = size(prm{1},2);
    rawdata = false;
else
    error('Input must be a cell array of vectors, or a cell array of 5xN matrices');
end

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('prm');
ip.addOptional('Annotations', [], @(x) isempty(x) || size(x,2)==2);
ip.addParamValue('FaceColor', jet(nb), @(x) any(size(x,1)==[1 nb ng]));
ip.addParamValue('EdgeColor', []);
ip.addParamValue('BorderWidth', [], @isscalar); 
ip.addParamValue('XLabel', [], @ischar);
ip.addParamValue('YLabel', ' ', @ischar);
ip.addParamValue('YLim', [], @(x) numel(x)==2);
ip.addParamValue('YTick', []);
ip.addParamValue('XTickLabel', []);%, @(x) isempty(x) || (iscell(x) && (numel(x)==sum(nb) || numel(x)==ng)));
ip.addParamValue('BarWidth', 0.8, @isscalar);
ip.addParamValue('GroupDistance', 0.8, @isscalar);
ip.addParamValue('LineWidth', 1, @isscalar);
ip.addParamValue('Angle', 45, @(x) isscalar(x) && (0<=x && x<=90));
ip.addParamValue('ErrorBarWidth', 0.2, @(x) 0<x && x<=1);
ip.addParamValue('Parent', gca, @ishandle);
ip.addParamValue('Interpreter', 'tex', @(x) any(strcmpi(x, {'tex', 'latex', 'none'})));
ip.addParamValue('X', [], @(x) numel(x)==nb); % cell array of x-coordinates (groups only)
ip.addParamValue('AdjustFigure', true, @islogical);
ip.addParamValue('ErrorbarColor', []);
ip.addParamValue('PlotDensity', false, @islogical);
ip.addParamValue('DetectOutliers', true, @islogical);
ip.parse(prm, varargin{:});

% identify outliers
if ip.Results.DetectOutliers && rawdata
    [prm, outliers] = cellfun(@(i) detectOutliers(i), prm, 'unif', 0);
else
    outliers = cell(nb,ng);
end

if rawdata
    if ip.Results.PlotDensity
        [f,xi] = cellfun(@ksdensity, prm, 'unif', 0);
    end
    % calculate percentiles
    prm = cellfun(@(i) [prctile(i(:), [50 25 75]) min(i) max(i)]', prm, 'unif', 0);
    % concatenate
    prm = mat2cell(cat(2,prm{:}), 5, nb*ones(1,ng));
end

faceColor = ip.Results.FaceColor;
if size(faceColor,1)==1
    faceColor = repmat(faceColor, [nb 1]);
end

edgeColor = ip.Results.EdgeColor;
if size(edgeColor,1)==1
    edgeColor = repmat(edgeColor, [nb 1]);
elseif isempty(edgeColor)
    edgeColor = zeros(size(faceColor));
end

errorbarColor = ip.Results.ErrorbarColor;
if isempty(errorbarColor)
    errorbarColor = zeros(size(faceColor));
end

ha = ip.Results.Parent;
bw = ip.Results.BarWidth;
dg = ip.Results.GroupDistance; % distance between groups, in bar widths

xa = cell(1,ng);
if isempty(ip.Results.X)
    xa{1} = 1:nb;
    for g = 2:ng
        xa{g} = xa{1} + xa{g-1}(end) + dg;
    end
else
    dx = min(diff(ip.Results.X));
    for g = 1:ng
        w = (nb-1)/2;
        xa{g} = ip.Results.X(g) + (g-1)*dg + (-w:w)*dx/nb;
    end
end

if isempty(ip.Results.BorderWidth)
    if ng>1
        border = bw/2+dg/2;
    else
        border = 1-bw/2;
    end
else
    border = ip.Results.BorderWidth;
end


hold on;
h = zeros(1,nb); % handles for legend
topval = cellfun(@(i) max(i,[],1), prm, 'unif', 0); % keep track of highest value
topval = [topval{:}];
for g = 1:ng
    
    if size(prm{g},1)==6 % contains SEM
        p25 = prm{g}(3,:);
        p75 = prm{g}(4,:);
    else
        p25 = prm{g}(2,:);
        p75 = prm{g}(3,:);
    end
    
    % whiskers (plot first to mask bar at '0')
    if size(prm{g},1)==6
        w1 = prm{g}(5,:);
        w2 = prm{g}(6,:);
    else
        w1 = prm{g}(4,:);
        w2 = prm{g}(5,:);
    end
    
    % the box
    lb = xa{g} - bw/2;
    rb = xa{g} + bw/2;
    xv = [lb; rb; rb; lb; lb; rb];
    yv = [p75; p75; p25; p25; p75; p75];
    
    for b = 1:nb
        if size(faceColor,1)==nb
            ci = b;
        else
            ci = g;
        end

        if ip.Results.PlotDensity
            f{b,g} = f{b,g}/max(f{b,g})/2.5;
            fill([xa{g}(b)+f{b,g} xa{g}(b)-f{b,g}(end:-1:1)], [xi{b,g} xi{b,g}(end:-1:1)], 0.8*[1 1 1],...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end

        % plot whiskers
        wopts = {'Color', errorbarColor(ci,:), 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off'};
        he = errorbar(ha, xa{g}(b), p25(b), w1(b)-p25(b), 0, 'LineStyle', 'none', wopts{:});
%         setErrorbarStyle(he, ip.Results.ErrorBarWidth, 'Position', 'bottom');
        he = errorbar(ha, xa{g}(b), p75(b), 0, w2(b)-p75(b), 'LineStyle', 'none', wopts{:});
%         setErrorbarStyle(he, ip.Results.ErrorBarWidth, 'Position', 'top');
             
        hp = patch(xv(:,b), yv(:,b), faceColor(ci,:), 'EdgeColor', edgeColor(ci,:),...
            'LineWidth', ip.Results.LineWidth, 'Parent', ha);
        if g==1 % for legend
            h(b) = hp;
        end
        
        % mean/median line
        line([lb(b); rb(b)], prm{g}(1,b)*[1; 1], wopts{:}, 'Parent', ha);
        
        % plot outliers
        if ~isempty(outliers{b,g})
            plot(xa{g}(b), outliers{b,g}, 'kp', 'Color', 0.4*[1 1 1], 'MarkerFaceColor', 0.4*[1 1 1]);
        end
        
        % replot border
        hp = patch(xv(:,b), yv(:,b), faceColor(ci,:), 'EdgeColor', edgeColor(ci,:),...
            'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off', 'Parent', ha);
        set(hp, 'FaceColor', 'none');
    end
    
    % SEM
    if size(prm{g},1)==6
        he = errorbar(xa{g}, prm{g}(1,:), prm{g}(2,:), 'k', 'LineStyle', 'none', 'LineWidth',...
            ip.Results.LineWidth, 'HandleVisibility', 'off');
        setErrorbarStyle(he, 0.15);
    end
end

if ng>1
    XTick = arrayfun(@(k) (xa{k}(1) + xa{k}(end))/2, 1:ng);
else
    XTick = [xa{:}];
end

% position of the bars
xa = [xa{:}];

XTickLabel = ip.Results.XTickLabel;
if isempty(XTickLabel)
    if ng>1
        XTickLabel = 1:ng;
    else
        XTickLabel = 1:nb;
    end
end

XLim = [xa(1)-border xa(end)+border];
set(ha, 'XTick', XTick, 'XTickLabel', XTickLabel, 'XLim', XLim);
if ~isempty(ip.Results.YLim);
    YLim = ip.Results.YLim;
    set(ha, 'YLim', YLim);
else
    YLim = get(gca, 'YLim');
end
if ~isempty(ip.Results.YTick)
    set(ha, 'YTick', ip.Results.YTick);
end

% add annotation links (for significance etc.)
av = ip.Results.Annotations;
if ~isempty(av)
    pos = get(ha, 'Position');
    fpos = get(gcf, 'Position');
    dy = 0.25 * diff(YLim) / (diff(XLim) * (pos(4)*fpos(4))/(pos(3)*fpos(3)));

    maxposCount = zeros(nb*ng,1);
    for k = 1:size(av,1)
        y0 = max(topval(av(k,1):av(k,2)));
        maxpos = find(topval==y0, 1, 'first');
        maxposCount(maxpos) = maxposCount(maxpos)+1;
        plot(ha, xa(av(k,[1 1 2 2])), y0+dy+2.5*dy*(maxposCount(maxpos)-1)+[0 dy dy 0], 'k',...
            'LineWidth', 0.75*ip.Results.LineWidth);
    end
end
hold off;

% x labels
if ip.Results.Angle~=0 && ~isempty(ip.Results.XTickLabel)
    rotateXTickLabels(ha, 'Angle', ip.Results.Angle, 'Interpreter', ip.Results.Interpreter,...
        'AdjustFigure', ip.Results.AdjustFigure);
end


function [inliers, outliers] = detectOutliers(obs)

% Same methods as boxplot(). Outliers: values that are
% larger than q3 + w(q3 - q1) or smaller than q1 - w(q3 - q1),
% where q1 and q3 are the 25th and 75th percentiles, respectively.
% The default w = 1.5 corresponds to approximately +/- 2.7 SD and 99.3
% coverage if the data are normally distributed.
    
q = prctile(obs, [25 75]);
d = 1.5*(q(2)-q(1));
% outliers:
idx = obs > q(2)+d | obs < q(1)-d;

inliers = obs(~idx);
outliers = obs(idx);
