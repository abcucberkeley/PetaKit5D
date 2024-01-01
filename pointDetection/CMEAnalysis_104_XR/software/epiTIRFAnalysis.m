% epiTIRFAnalysis(data, labels, MaxIntensityThreshold, varargin) runs epi-TIRF analysis on the input data sets
%
% Inputs:
%       data : cell array of data sets
%     labels : cell array of labels for the data sets (e.g., {'Control', 'Perturbation'})
%     MaxIntensityThreshold : threshold from runLifetimeAnalysis.m
%
% Optional inputs:
%    Display : true|{false} show intensity scaling
%
% This function reproduces the epi:TIRF analysis from Aguet et al., Dev Cell, 2013
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

% Author: Francois Aguet

function epiTIRFAnalysis(data, labels, MaxIntensityThreshold, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @iscell);
ip.addRequired('labels', @iscell)
ip.addRequired('MaxIntensityThreshold', @isscalar);
ip.addParamValue('Display', false, @islogical);
ip.addParamValue('CurvaRatioThresh', 1.5, @isnumeric);
ip.parse(data, labels, MaxIntensityThreshold, varargin{:});

fset = loadFigureSettings('print');
rv = 0:0.05:3;
lv = 0:2:120;
T = ip.Results.MaxIntensityThreshold;

cv = hsv2rgb([0 0 0;
              0.99 1 1;              
              0.55 1 1]);

nexp = numel(data);
res = struct([]);
for ei = 1:nexp

	lftData = getLifetimeData(data{ei}, 'ReturnValidOnly', true, 'Cutoff_f', 5, 'ExcludeVisitors', true);

	maxA = arrayfun(@(i) max(i.A(:,:,1),[],2), lftData, 'unif', 0);
	a = scaleEDFs(maxA, 'Display', ip.Results.Display);

	% apply scaling to both channels
	for i = 1:numel(lftData)
	    lftData(i).A =   a(i) * lftData(i).A;
	    lftData(i).sbA = a(i) * lftData(i).sbA;
	    lftData(i).ebA = a(i) * lftData(i).ebA;
	end

	A = vertcat(lftData.A);
	maxTIR = nanmax(A(:,:,1),[],2);
	maxEpi = nanmax(A(:,:,2),[],2); 
	
	% per-cell arrays for lifetime analysis
	res(ei).iMaxTIR = arrayfun(@(i) nanmax(i.A(:,:,1),[],2), lftData, 'unif', 0);
	res(ei).iMaxEpi = arrayfun(@(i) nanmax(i.A(:,:,2),[],2), lftData, 'unif', 0);

	%%
	%==================================================
	% Panel D: Epi/TIRF ratio vs. lifetime
	%==================================================
	% ratio and lifetime for each CCS
	rLS = sum(A(:,1:5,1).*A(:,1:5,2),2)./sum(A(:,1:5,1).^2,2);
	ratioVec = maxEpi ./ maxTIR / mean(rLS);
	lftVec = vertcat(lftData.lifetime_s);

	% CCPs > R=1+ci and in [1-ci,1+ci]
	nd = numel(lftData);
	cellIdx = arrayfun(@(i) i*ones(numel(lftData(i).lifetime_s),1), 1:nd, 'unif', 0);
	cellIdx = vertcat(cellIdx{:});
	res(ei).iRatio = arrayfun(@(i) ratioVec(cellIdx==i), 1:max(cellIdx), 'unif', 0);

	% remove CSs and CCSs with negative ratio, if any
	idx = maxTIR>=T & ratioVec>0;
	lftVec = lftVec(idx);
	ratioVec = ratioVec(idx);
	res(ei).ratioVec = ratioVec;
	res(ei).label = labels(ei);

	% 99% confidence interval for ratio=1, calculated from lower half (ratio <1) of the distribution
	tmp = ratioVec(ratioVec<1 & lftVec>20)-1;
	ciC = norminv(1-0.01/2,0,1)*std(tmp)/sqrt(1-2/pi); % correction w/ half-normal s.d.

	na = arrayfun(@(i) sum(ratioVec(cellIdx(idx)==i)>1+ciC), 1:nd);
	nb = arrayfun(@(i) sum(ratioVec(cellIdx(idx)==i)>1-ciC & ratioVec(cellIdx(idx)==i)<1+ciC), 1:nd);
	fprintf('Curved CCPs (ratio > %.1f) WT:  %.1f±%.1f, flat CCPs: %.1f±%.1f (%.1f±%.1f%%; %d cells)\n',...
	    1+ciC, mean(na), std(na), mean(nb), std(nb), mean(na./nb)*100, std(na./nb)*100, nd);

	figure(fset.fOpts{:}, 'Position', [15 5 6.5 6.5]);
	axes(fset.axOpts{:}, 'Position', [1.5 1.5 4.5 4.5], 'TickLength', fset.TickLength/4.5*6);
	hold on;
	[~, rangeCtrl] = densityplot(lftVec(lftVec<=lv(end)), ratioVec(lftVec<=lv(end)), lv, rv,...
		'DisplayFunction', @sqrt);
	plot([0 120], 1-ciC*[1 1], 'w--', 'LineWidth', 1);
	plot([0 120], 1+ciC*[1 1], 'w--', 'LineWidth', 1);
	set(gca, 'XTick', 0:20:120);
	formatTickLabels(gca);
	xlabel('Lifetime (s)', fset.lfont{:});
	ylabel('Epi:TIR ratio', fset.lfont{:});
    
    numFLatSL=sum((lftVec<20)&(ratioVec>0.5)&(ratioVec<1.5))/length(lftVec);
    disp(['Flat (0.5 < Ratio < 1.5) and Short Lived CCPs (<20s) ' num2str(numFLatSL)]);
end


%==================================================
% Compare lifetimes for the first two conditions
%==================================================
rRange = ip.Results.CurvaRatioThresh;
% rRange = 1+max(ciC,ciD);
for ri = 1:numel(rRange)
    rf = rRange(ri);
    
    % create cell arrays for track indexes
    ratio1GT = cellfun(@(x,y) x>=rf & y>=T, res(1).iRatio, res(1).iMaxTIR, 'unif', 0);
    ratio1LT = cellfun(@(x,y) x<rf & y>=T, res(1).iRatio, res(1).iMaxTIR, 'unif', 0);
    ratio2GT = cellfun(@(x,y) x>=rf & y>=T, res(2).iRatio, res(2).iMaxTIR, 'unif', 0);
    ratio2LT = cellfun(@(x,y) x<rf & y>=T, res(2).iRatio, res(2).iMaxTIR, 'unif', 0);
    opts = {'MaxIntensityThreshold', T, 'Display', 'off', 'Rescale', false,...
        'ExcludeVisitors', true, 'RemoveOutliers', false, 'Overwrite', false,'CorrectObservationBias',true};
    lftRes1GT = runLifetimeAnalysis(data{1}, opts{:}, 'SelectIndex', ratio1GT);
    lftRes1LT = runLifetimeAnalysis(data{1}, opts{:}, 'SelectIndex', ratio1LT);
    lftRes2GT = runLifetimeAnalysis(data{2}, opts{:}, 'SelectIndex', ratio2GT);
    lftRes2LT = runLifetimeAnalysis(data{2}, opts{:}, 'SelectIndex', ratio2LT);
	
    figure(fset.fOpts{:}, 'Position', [5 5 8 7]);
    axes(fset.axOpts{:});
    hold on;
    h(1) = plot(lftRes1GT.t, lftRes1GT.meanLftHistCCP, 'k', 'LineWidth', 1);
    h(2) = plot(lftRes1LT.t, lftRes1LT.meanLftHistCCP, '-', 'Color', 0.6*[1 1 1], 'LineWidth', 1);
    
    h(3) = plot(lftRes2GT.t, lftRes2GT.meanLftHistCCP, 'Color', cv(3,:), 'LineWidth', 1);
    h(4) = plot(lftRes2LT.t, lftRes2LT.meanLftHistCCP, 'Color', cv(2,:), 'LineWidth', 1);

    axis([0 120 0 0.06]);
    set(gca, 'YTick', 0:0.01:0.06);
    
    na1 = cellfun(@(x) sum(x), ratio1GT);
    na2 = cellfun(@(x) sum(x), ratio1LT);
    na3 = cellfun(@(x) sum(x), ratio2GT);
    na4 = cellfun(@(x) sum(x), ratio2LT);
    rstd = @(x) 1/norminv(0.75, 0, 1) * mad(x, 1);
    rmean = @(x) round(mean(x)/100)*100;
    
    str = num2str(rf, '%.1f');
    fmt = '%.1f';
    fmt2 = '%.2f';
    hl = legend(h,...
        [' ' labels{1} ', Epi:TIR > ' str ' (~' num2str(mean(na1./(na1+na2))*100,fmt) '% CCPs)'],...
        [' ' labels{1} ', Epi:TIR < ' str ],...
        [' ' labels{2} ', Epi:TIR > ' str ' (~' num2str(mean(na3./(na3+na4))*100,fmt) '% CCPs)'],...
        [' ' labels{2} ', Epi:TIR < ' str]);
    set(hl, 'Box', 'off', 'Position', [1.5 5.25 2 1.5]);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
end
