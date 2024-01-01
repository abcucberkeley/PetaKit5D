%analyzeBleaching(data) generates plots of mean intensity values over time 
% to characterize bleaching
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

% Francois Aguet, 08/31/2011

function analyzeBleaching(data)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.parse(data);

nd = length(data);
% res(1:nd) = struct('t', [], 'intensity', [], 'model', [], 'k', [], 'a', [], 'd', []);
for i = 1:nd

    % load data (much faster than loading frames individually)
    frames = readtiff(data.framePaths{1});
    
    % generate/load cell mask
    cellMask = getCellMask(data)~=0;
    
    % load detection mask, if it exists
    ccpMasks = [];
    if exist(data.maskPaths, 'file')==2
        ccpMasks = readtiff(data.maskPaths)~=0;
    end


    t = (0:data.movieLength-1)*data.framerate;
    rawInt = squeeze(mean(mean(frames)));

    
    cellMaskInt = zeros(1,data.movieLength);
    cytosolInt = zeros(1,data.movieLength);
    for f = 1:data.movieLength
        iframe = double(frames(:,:,f));
        cellMaskInt(f) = mean(iframe(cellMask));
        if ~isempty(ccpMasks)
            imask = cellMask & ~ccpMasks(:,:,f);
            cytosolInt(f) = mean(iframe(imask));
        end
    end

    if ~isempty(ccpMasks)
        [k, a, d] = fitDoubleExp(t, cytosolInt);
    else
        [k, a, d] = fitDoubleExp(t, cellMaskInt);
    end
        
    figure;
    hold on;
    plot(t, rawInt, 'Color', 0.6*[1 1 1]);
    plot(t, cellMaskInt, 'g');
    if ~isempty(ccpMasks)
        plot(t, cytosolInt, 'r');
    end
    plot(t, expModel(k,a,d,t), 'k');
    % YLim = get(gca, 'YLim');
    % set(gca, 'YLim', [0, YLim(2)]);
    xlabel('Time (s)');
    ylabel('Mean intensity');
    
    if ~isempty(ccpMasks)
        hl = legend(' Mean frame intensity', ' Mean cell intensity',...
            ' Mean cytosol intensity', ' Exp. fit');
    else
        hl = legend('Mean frame intensity', ' Mean cell intensity', ' Exp. fit');
    end
    set(hl, 'EdgeColor', 'w');
    
    %res.t = t;
    %res.model = expModel(k,a,d,t);
    %res.k = k;
    %res.a = a;
    %res.d = d;
end
    

function [k, a, d] = fitDoubleExp(t, f)

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

dy = f(end);
mu = sum(t.*(f-dy))/sum(f-dy);
init = [1/mu f(1)-dy 1/mu f(1)-dy dy];

lb = zeros(1,5);
ub = Inf(1,5);

% parameters: k1 a1 k2 a2 dy
p = lsqnonlin(@(p) expModel(p([1 3]), p([2 4]), p(5), t) - f,...
    init, lb, ub, opts);
k = p([1 3]);
a = p([2 4]);
d = p(5);


% double exponential model (with constant offset for background);
function m = expModel(k, a, d, t)
m = a(1)*exp(-k(1)*t) + a(2)*exp(-k(2)*t) + d;
