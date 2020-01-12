%[psnr] = getPSNRDistribution(data, varargin) plots the PSNR distribution for the detections in 'data'.
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

% Francois Aguet, 12/18/12

function varargout = plotPSNRDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('ha', [], @isnumeric);
ip.addParamValue('Color', hsv2rgb([0.33 1 0.9]));
ip.addParamValue('XLim', [1 300]);
ip.addParamValue('Pool', true, @islogical);
ip.addParamValue('DisplayMode', 'screen', @(x) any(strcmpi(x, {'screen', 'print'})));
ip.addParamValue('Mode', 'all', @(x) any(strcmpi(x, {'all', 'max'})));
ip.addParamValue('Channel', 1);
ip.parse(data, varargin{:});
ha = ip.Results.ha;

nd = numel(data);
psnr = cell(1,nd);
parfor i = 1:nd
    psnr{i} = getPSNRDistribution(data(i), 'Mode', ip.Results.Mode, 'Channel', ip.Results.Channel); %#ok<PFBNS>
end

fset = loadFigureSettings(ip.Results.DisplayMode);
if isempty(ha)
    figure(fset.fOpts{:}, 'Name', 'PSNR distribution');
    ha = axes(fset.axOpts{:}, 'xscale', 'log');
    hold on;
end
% PSNR values, in [dB]
dx = 0.2;
xi = 10.^((0:dx:30)/10);

cmap = jet(nd);
[~,idx] = sort(cellfun(@median, psnr), 'ascend');
cmap = cmap(idx,:);

if ~ip.Results.Pool
    
    hp = zeros(1,nd);
    for i = 1:numel(data)
        fi = hist(psnr{i}, xi)/numel(psnr{i});
        hp(i) = plot(xi, fi, '-', 'Color', cmap(i,:), 'LineWidth', 1);
    end
    ltext = arrayfun(@getMovieName, data, 'unif', 0);
    hl = legend(hp(idx), ltext(idx));
    set(hl, 'Interpreter', 'none', fset.tfont{:}, 'Box', 'on', 'EdgeColor', 'w');
else
    psnr = [psnr{:}];
    fi = hist(psnr, xi)/numel(psnr);
    plot(xi, fi, '-', 'Color', ip.Results.Color, 'LineWidth', 1);
end

set(ha, 'XLim', ip.Results.XLim, 'XTick', 10.^(0:2), 'XTickLabel', {'1', '10', '100'},...
    fset.sfont{:});
xlabel('PSNR', fset.sfont{:});
ylabel('Frequency', fset.sfont{:});
formatTickLabels();

if nargout>0
    varargout{1} = psnr;
end
if nargout>1
    varargout{2} = cmap;
end
