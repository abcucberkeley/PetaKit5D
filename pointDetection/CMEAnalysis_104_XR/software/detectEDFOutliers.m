%[outlierIdx] = detectEDFOutliers(samples, varargin) identifies outlier data sets based on their empirical distribution function
%
% Inputs:
%         samples : cell array of sample vectors
%
% Options:
%          offset : percent of missing data for each data set (vector of size numel(samples))
%          refIdx : index of the reference data set
%
% Outputs:
%      outlierIdx : index of outlier data sets
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

% Francois Aguet, 06/08/2012 (modified on 09/05/2012)

function outlierIdx = detectEDFOutliers(samples, varargin)

ns = numel(samples);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('samples', @iscell);
ip.addOptional('offset', zeros(1,ns));
ip.addOptional('refIdx', []);
ip.addParamValue('Display', true, @islogical);
ip.addParamValue('FigureName', 'Outlier analysis');
ip.parse(samples, varargin{:});
offset = ip.Results.offset;
medIdx = ip.Results.refIdx;
xEDF = cell(1,ns);
fEDF = cell(1,ns);
for i = 1:ns
    [fEDF{i}, xEDF{i}] = ecdf(samples{i});
end

% perform calculation over range common to all input EDFs
minV = max(cellfun(@(i) min(i), xEDF));
maxV = min(cellfun(@(i) max(i), xEDF));
x = vertcat(xEDF{:});
x = unique(x(minV<=x & x<=maxV))';

% interpolate all EDFs to common coordinates
for i = 1:ns
    fEDF{i} = interp1(xEDF{i}(2:end), offset(i) + (1-offset(i))*fEDF{i}(2:end), x);
end
M = vertcat(fEDF{:});

% Median EDF, approach (1): majority of median indexes
% medianEDF = median(M,1);
% if isempty(medIdx)
%     J = nansum((M-repmat(medianEDF, [ns 1])).^2, 2);
%     medIdx = find(J==min(J),1,'first');
% end
% medEDF = M(medIdx,:);

% Median EDF, approach (2): median at each coordinate
medEDF = median(M,1);


% MAD: median(abs(X-median(X)))
% sigma = 1/norminv(0.75) * mad(M, 1, 1); % equiv. for approach (2)
sigma = 1/norminv(0.75) * median(abs(M-repmat(medEDF,[ns 1])),1);

% smoothen the envelope (somewhat arbitrary)
xi = linspace(x(1), x(end), 1000);
si = filterGauss1D(interp1(x, sigma, xi, 'pchip'), 15);
si = interp1(xi, si, x);

% 3*sigma bounds
ub = medEDF + 3*si;
lb = medEDF - 3*si;

% for each EDF, test whether 95% of its points fall within median ± 3*sigma
outlierIdx = zeros(1,ns);
for i = 1:ns
    tmp = M(i,:)<lb | ub<M(i,:);
    outlierIdx(i) = sum(tmp)/numel(tmp) > 0.05;    
end
outlierIdx = find(outlierIdx);

if ~isempty(outlierIdx) && ip.Results.Display
    fset = loadFigureSettings();
    hp = NaN(1,4);
    figure('Name', ip.Results.FigureName);
    hold on;
    hp(4) = fill([x x(end:-1:1)], [lb ub(end:-1:1)], 0.6*[1 1 1], 'EdgeColor', 'none');
    
    if ~isempty(outlierIdx)
        hx = plot(x, M(outlierIdx,:), 'Color', [1 0 0], 'LineWidth', 2);
        hp(2) = hx(1);
    end
    hx = plot(x, M(setdiff(1:ns, [outlierIdx medIdx]),:), 'Color', [0 1 0], 'LineWidth', 2);
    hp(1) = hx(1);
    hp(3) = plot(x, medEDF, 'k', 'LineWidth', 2);

    set(gca, fset.sfont{:}, fset.axOpts{:},...
        'YTick', 0:0.1:1, 'YLim', [0 1.01]);
    xlabel('Time (s)');
    ylabel('F(t)');
    if isnan(hp(2))
        hp(2) = [];
        hl = legend(hp, ' Inlier distribution', ' Median of distributions', ' Critical area', 'Location', 'SouthEast');
    else
        hl = legend(hp, ' Inlier distribution', ' Outlier distribution', ' Median of distributions', ' Critical area', 'Location', 'SouthEast');
    end
    set(hl, fset.tfont{:}, 'Box', 'off');
    pos = get(hl, 'Position');
    pos(2) = 0.2;
    set(hl, 'Position', pos);
    drawnow;
end
