%runDetection(data) detects CCPs using a combination of model-based (PSF) fitting and statistical tests
%
% Inputs:      data : data/movie structure
%         {'Sigma'} : standard deviation of the Gaussian used for fitting
%     {'Overwrite'} : true | {false}
%
% Notes: All input data sets must have the same channels.
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

% Francois Aguet, April 2011 (last modified 08/18/2013)

function runDetection(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Sigma', [], @(x) isempty(x) || numel(x)==length(data(1).channels));
ip.addParamValue('SigmaSource', 'data', @(x) any(strcmpi(x, {'data', 'model'})));
ip.addParamValue('RemoveRedundant', true, @islogical);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('Master', [], @isnumeric);
ip.addParamValue('Alpha', 0.05, @isnumeric);
ip.addParamValue('CellMask', [],  @isnumeric);
ip.parse(data, varargin{:});
overwrite = ip.Results.Overwrite;
mCh = ip.Results.Master;
if isempty(mCh)
    mCh = unique(arrayfun(@(i) find(strcmpi(i.channels, i.source)), data));
    if numel(mCh)>1
        error('Master channel index is inconsistent accross data sets.');
    end
end

hasDet = arrayfun(@(i) exist([i.channels{mCh} 'Detection' filesep 'detection_v2.mat'], 'file')==2, data);

nd = numel(data);
sigma = ip.Results.Sigma;
if isempty(sigma) && (~all(hasDet) || overwrite)
    % verify that all data sets have same channels
    nCh = unique(arrayfun(@(i) numel(i.channels), data));
    markers = arrayfun(@(c) unique(arrayfun(@(i) i.markers{c}, data, 'unif', 0)), 1:nCh, 'unif', 0);
    if numel(nCh)==1 && all(cellfun(@(x) numel(x)==1, markers))
        if strcmpi(ip.Results.SigmaSource, 'model')
            sigma = getGaussianPSFsigma(data(1).NA, data(1).M, data(1).pixelSize, data(1).markers);
        else
            fprintf('Determining Gaussian PSF parameters from data ... ');
            % evenly sample all data sets
            sigma = zeros(1,nCh);
            % use ~20 frames distributed accross data sets
            nf = round(40/nd);
            for c = 1:nCh
                frames = cell(nd, nf);
                for i = 1:nd
                    fidx = round(linspace(1,data(i).movieLength,nf));
                    if iscell(data(i).framePaths{c})
                        frames(i,:) = arrayfun(@(f) double(imread(data(i).framePaths{c}{f})), fidx, 'unif', 0);
                    else
                        frames(i,:) = arrayfun(@(f) double(readtiff(data(i).framePaths{c}, f)), fidx', 'unif', 0); 
                    end
                end
                sigma(c) = getGaussianPSFsigmaFromData(vertcat(frames(:)), 'Display', true);
            end
            fprintf('done.\n');
        end
        fprintf('Gaussian PSF s.d. values: ');
        fprintf(' %.2f', sigma);
        fprintf('\n');            
    else
        error('runDetection error: mismatch between the channel fluorophores in ''data''. Could not estimate ''sigma''.');
    end
end
if any(sigma<1.1)
    fprintf(2, 'Sigma values < 1.1 were rounded to 1.1 to avoid poor localization performance.\n');
    sigma(sigma<1.1) = 1.1;
end

for i = 1:nd
    if ~hasDet(i) || overwrite
        fprintf('Running detection for %s ...', getShortPath(data(i)));
        main(data(i), sigma, mCh, ip.Results);
        fprintf(' done.\n');
    else
        fprintf('Detection has already been run for %s\n', getShortPath(data(i)));
    end
end



function main(data, sigma, mCh, opts)

% master channel
nCh = length(data.channels);

frameInfo(1:data.movieLength) = struct('x', [], 'y', [], 'A', [], 's', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'A_pstd', [], 'c_pstd', [],...
    'x_init', [], 'y_init', [], 'maskA', [], 'maskN', [],...
    'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_Ar', [],  'mask_Ar', [], 'hval_Ar', [],  'hval_AD', [], 'isPSF', [],...
    'xCoord', [], 'yCoord', [], 'amp', [], 'dRange', []);

[~,~] = mkdir([data.channels{mCh} 'Detection']);
if iscell(data.framePaths{mCh})
    [~,~] = mkdir([data.channels{mCh} 'Detection' filesep 'Masks']);
end

% double fields, multi-channel
dfields = {'x', 'y', 'A', 'c', 'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'RSS', 'pval_Ar'};
% logical fields, multi-channel
lfields = {'hval_Ar', 'hval_AD', 'isPSF'};
% slave channel fields
sfields = [dfields {'hval_Ar', 'hval_AD'}]; % 'isPSF' is added later

rmfields = [dfields lfields {'x_init', 'y_init', 'maskA', 'maskN', 'mask_Ar'}];

fmt = ['%.' num2str(ceil(log10(data.movieLength+1))) 'd'];
mask = false([data.imagesize data.movieLength]);
parfor k = 1:data.movieLength
    if ~iscell(data.framePaths{mCh}) %#ok<PFBNS>
        img = double(readtiff(data.framePaths{mCh}, k));
    else
        img = double(imread(data.framePaths{mCh}{k}));
    end
    
    [pstruct, mask(:,:,k)] = pointSourceDetection(img, sigma(mCh), 'Alpha', opts.Alpha,...
        'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant); %#ok<PFBNS>
    
    if ~isempty(pstruct)
        pstruct.s = sigma;
        pstruct = rmfield(pstruct, 's_pstd');
        
        pstruct.dRange{mCh} = [min(img(:)) max(img(:))];
        np = numel(pstruct.x);
        
        % expand structure for slave channels
        for f = 1:length(dfields)
            tmp = NaN(nCh, np);
            tmp(mCh,:) = pstruct.(dfields{f});
            pstruct.(dfields{f}) = tmp;
        end
        for f = 1:length(lfields)
            tmp = false(nCh, np);
            tmp(mCh,:) = pstruct.(lfields{f});
            pstruct.(lfields{f}) = tmp;
        end
        
        % get component size and intensity for each detection
        CC = bwconncomp(mask(:,:,k));
        labels = labelmatrix(CC);
        % mask label for each detection
        loclabels = labels(sub2ind(size(img), pstruct.y_init, pstruct.x_init));
        compSize = cellfun(@(i) numel(i), CC.PixelIdxList);
        pstruct.maskN = compSize(loclabels);
        compInt = cellfun(@(i) sum(img(i))/numel(i), CC.PixelIdxList);
        pstruct.maskA = compInt(loclabels);
               
        for ci = setdiff(1:nCh, mCh)
            if ~iscell(data.framePaths{ci})
                img = double(readtiff(data.framePaths{ci}, k));
            else
                img = double(imread(data.framePaths{ci}{k}));
            end
            pstruct.dRange{ci} = [min(img(:)) max(img(:))];
            pstructSlave = fitGaussians2D(img, pstruct.x(mCh,:), pstruct.y(mCh,:), [], sigma(ci)*ones(1,np), [], 'Ac');
            
            % localize, and compare intensities & (x,y)-coordinates. Use localization result if it yields better contrast
            pstructSlaveLoc = fitGaussians2D(img, pstruct.x(mCh,:), pstruct.y(mCh,:), pstructSlave.A, sigma(ci)*ones(1,np), pstructSlave.c, 'xyAc');
            idx = sqrt((pstruct.x(mCh,:)-pstructSlaveLoc.x).^2 + (pstruct.y(mCh,:)-pstructSlaveLoc.y).^2) < 3*sigma(mCh) & pstructSlaveLoc.A > pstructSlave.A;
            
            % fill slave channel information
            for f = 1:length(sfields)
                pstruct.(sfields{f})(ci,~idx) = pstructSlave.(sfields{f})(~idx);
                pstruct.(sfields{f})(ci,idx) = pstructSlaveLoc.(sfields{f})(idx);
            end
            
            nanIdx = isnan(pstructSlave.x); % points within slave channel border, remove from detection results
            for f = 1:length(rmfields)
                pstruct.(rmfields{f})(:,nanIdx) = [];
            end
            np = size(pstruct.x,2);
            
            pstruct.isPSF(ci,:) = ~pstruct.hval_AD(ci,:);
        end
        
        % add fields for tracker
        pstruct.xCoord = [pstruct.x(mCh,:)' pstruct.x_pstd(mCh,:)'];
        pstruct.yCoord = [pstruct.y(mCh,:)' pstruct.y_pstd(mCh,:)'];
        pstruct.amp =    [pstruct.A(mCh,:)' pstruct.A_pstd(mCh,:)'];
        frameInfo(k) = orderfields(pstruct, fieldnames(frameInfo(k))); %#ok<PFOUS>
    else
        frameInfo(k).dRange{mCh} = [min(img(:)) max(img(:))];
        for ci = setdiff(1:nCh, mCh)
            if ~iscell(data.framePaths{ci})
                img = double(readtiff(data.framePaths{ci}, k));
            else
                img = double(imread(data.framePaths{ci}{k}));
            end
            frameInfo(k).dRange{ci} = [min(img(:)) max(img(:))];
        end
    end
    if iscell(data.framePaths{1})
        maskPath = [data.channels{mCh} 'Detection' filesep 'Masks' filesep 'dmask_' num2str(k, fmt) '.tif'];
        imwrite(uint8(255*mask(:,:,k)), maskPath, 'tif', 'compression' , 'lzw');
    end
end

% write masks
if ~iscell(data.framePaths{1})
    imwrite(uint8(255*mask(:,:,1)), data.maskPaths, 'tif', 'compression' , 'lzw');
    for k = 2:data.movieLength
        imwrite(uint8(255*mask(:,:,k)), data.maskPaths, 'tif', 'compression' , 'lzw', 'writemode', 'append');
    end
end

save([data.channels{mCh} 'Detection' filesep 'detection_v2.mat'], 'frameInfo');
