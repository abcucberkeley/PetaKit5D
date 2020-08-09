%[sigma] = getGaussianPSFsigmaFromData(imageList) returns the s.d. of the Gaussian PSF estimated from the input data.
% The estimation is performed by running pointSourceDetection() with 'sigma' as a free parameter.
%
% Inputs:
%    imageList : single image, cell array of images, or cell array of image path strings
%
% Options:
%    'Display' : {true}|false displays the distribution of 'sigma' values
%
% Output:
%        sigma : standard deviation of the Gaussian PSF estimated from point sources in input data
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

% Francois Aguet, September 2010

function sigma = getGaussianPSFsigmaFromData(imageList, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imageList');
ip.addOptional('frameRange', []);
ip.addParamValue('Display', true, @islogical);
ip.parse(imageList, varargin{:});

if ~iscell(imageList)
    imageList = {imageList};
end
frameRange = ip.Results.frameRange;

if isempty(frameRange)
    nd = numel(imageList);
else
    nd = numel(frameRange);
end
svect = cell(1,nd);
parfor i = 1:nd
    if isempty(frameRange)
        if ischar(imageList{i})
            img = double(imread(imageList{i}));
        else
            img = double(imageList{i});
        end
    else
        img = double(readtiff(imageList{1}, frameRange(i)));
    end
    % First pass with fixed sigma
    pstruct = pointSourceDetection(img, 1.5, 'Mode', 'xyac');
    if ~isempty(pstruct)
        pstruct = fitGaussians2D(img, pstruct.x, pstruct.y, pstruct.A, 1.5*ones(1,length(pstruct.x)), pstruct.c, 'xyasc');
        isPSF = ~[pstruct.hval_AD] & [pstruct.pval_Ar] < 0.05;
        svect{i} = pstruct.s(~isnan(pstruct.s) & isPSF);
    end
end
svect = [svect{:}];

opts = statset('maxIter', 200);
try
    w = warning('off', 'stats:gmdistribution:FailedToConverge');
    obj = cell(1,3);
    for n = 1:3
        obj{n} = gmdistribution.fit(svect', n, 'Options', opts);
    end
    [~,idx] = min(cellfun(@(i) i.BIC, obj));
    obj = obj{idx};
    [mu,idx] = sort(obj.mu);
    svec = sqrt(squeeze(obj.Sigma(:,:,idx)));
    amp = obj.PComponents(idx)';
    %[~,idx] = max(amp);
    [~,idx] = max(amp./(sqrt(2*pi)*svec));
    sigma = mu(idx);

    warning(w);

    if ip.Results.Display
        si = linspace(0, prctile(svect,99.5), 101);
        ni = histc(svect, si);
        ds = si(2)-si(1);
        ni = ni/sum(ni)/ds;
        
        setupFigure('DisplayMode', 'screen');
        h = bar(si,ni,'histc');
        
        set(h, 'FaceColor', 0.6*[1 1 1], 'EdgeColor', 'none');
        
        hold on;
        plot(si, pdf(obj,si'), 'r');
        for i = setdiff(1:3, idx)
            plot(si, amp(i)*normpdf(si, mu(i), svec(i)), 'k', 'HandleVisibility', 'off');
        end
        plot(si, amp(idx)*normpdf(si, mu(idx), svec(idx)), 'g');
        plot(mu(idx)*[1 1], [0 amp(idx)/(sqrt(2*pi)*svec(idx))], 'g', 'LineWidth', 3);
        xlabel('Standard deviation of 2D Gaussian model (pixels)');
        ylabel('Frequency');
        hl = legend(' Measured distribution', ' Mixture model', 'Primary mode', ['\sigma = ' num2str(mu(idx), '%.2f')]); 
        set(hl, 'Box', 'off');
        formatTickLabels();
        drawnow;
    end
catch
    fprintf('Could not determine distribution, potentially due to insufficient samples.');
    sigma = mean(svect);
end
