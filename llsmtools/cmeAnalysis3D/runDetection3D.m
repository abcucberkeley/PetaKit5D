%runDetection3D(data) detects CCPs using a combination of model-based (PSF) fitting and statistical tests
%
% Inputs:      data : data/movie structure
%         {'Sigma'} : standard deviation of the Gaussian used for fitting
%     {'Overwrite'} : true | {false}
%
% Notes: 
%  - 3D coordinates are in pixel space, except for the zCoord field passed
%    to the tracker, which is corrected for z-anisotropy
%  - All input data sets must have the same channels.

% Francois Aguet, 08/2013

function runDetection3D(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
% ip.addOptional('destPath', [], @ischar);
ip.addParamValue('Sigma', []);%, @(x) numel(x)==length(data(1).channels));
% ip.addParamValue('SigmaSource', 'model', @(x) any(strcmpi(x, {'data', 'model'})));
ip.addParamValue('RemoveRedundant', true, @islogical);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('Master', [], @isnumeric);
ip.addParamValue('Alpha', 0.05, @isnumeric);
ip.addParamValue('CellMask', [],  @isnumeric);
ip.addParamValue('FileName', 'Detection3D.mat');
ip.addParamValue('WindowSize', []);
ip.addParamValue('ResultsPath', arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0), @iscell);
ip.parse(data, varargin{:});
overwrite = ip.Results.Overwrite;
mCh = ip.Results.Master;
if isempty(mCh)
    mCh = unique(arrayfun(@(i) find(strcmpi(i.channels, i.source)), data));
    if numel(mCh)>1
        error('Master channel index is inconsistent accross data sets.');
    end
end

sigma = ip.Results.Sigma;
% if isempty(sigma)
%     % verify that all data sets have same channels
%     nCh = unique(arrayfun(@(i) numel(i.channels), data));
%     markers = arrayfun(@(c) unique(arrayfun(@(i) i.markers{c}, data, 'unif', 0)), 1:nCh, 'unif', 0);
%     if numel(nCh)==1 && all(cellfun(@(x) numel(x)==1, markers))
%         if strcmpi(ip.Results.SigmaSource, 'model')
%             sigma = getGaussianPSFsigma(data(1).NA, data(1).M, data(1).pixelSize, data(1).markers);
%         else
%             % evenly sample all data sets
%             sigma = zeros(1,nCh);
%             % use ~20 frames distributed accross data sets
%             nf = round(20/numel(data));
%             for c = 1:nCh
%                 fpaths = arrayfun(@(i) i.framePaths{c}(round(linspace(1,i.movieLength,nf))), data, 'unif', 0);
%                 sigma(c) = getGaussianPSFsigmaFromData(vertcat(fpaths{:}), 'Display', false);
%             end
%         end
%         fprintf('Gaussian PSF s.d. values: ');
%         fprintf(' %.2f', sigma);
%         fprintf('\n');            
%     else
%         error('runDetection error: mismatch between the channel fluorophores in ''data''. Could not estimate ''sigma''.');
%     end
% end
if any(sigma<1.1)
    fprintf(2, 'Sigma values < 1.1 were rounded to 1.1 to avoid poor localization performance.\n');
    sigma(sigma<1.1) = 1.1;
end

for i = 1:length(data)
    if ~(exist([ip.Results.ResultsPath{i} ip.Results.FileName], 'file') == 2) || overwrite
        fprintf('Running detection for %s ...', getShortPath(data(i)));
                
        main(data(i), sigma, mCh, ip.Results.ResultsPath{i}, ip.Results);
        fprintf(' done.\n');
    else
        fprintf('Detection has already been run for %s\n', getShortPath(data(i)));
    end
end



function main(data, sigma, mCh, resultsPath, opts)

% master channel
nCh = length(data.channels);

frameInfo(1:data.movieLength) = struct('x', [], 'y', [], 'z', [], 'A', [], 's', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'z_pstd', [], 'A_pstd', [], 'c_pstd', [],...
    'x_init', [], 'y_init', [], 'z_init', [],...
    'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_Ar', [], 'hval_Ar', [],  'hval_AD', [], 'isPSF', [],...
    'xCoord', [], 'yCoord', [], 'zCoord', [], 'amp', [], 'dRange', []);

fmt = ['%.' num2str(ceil(log10(data.movieLength+1))) 'd'];
[~,~] = mkdir(resultsPath);
[~,~] = mkdir([resultsPath 'Masks']);

% double fields, multi-channel
dfields = {'x', 'y', 'z', 'A', 'c', 'x_pstd', 'y_pstd', 'z_pstd', 'A_pstd', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'RSS', 'pval_Ar'};
% logical fields, multi-channel
lfields = {'hval_Ar', 'hval_AD', 'isPSF'};
% slave channel fields
sfields = [dfields {'hval_Ar', 'hval_AD'}]; % 'isPSF' is added later

rmfields = [dfields lfields {'x_init', 'y_init', 'z_init', 'mask_Ar'}];

parfor k = 1:data.movieLength
    frame = double(readtiff(data.framePathsDS{mCh}{k})); %#ok<PFBNS>
    frame(frame==0) = NaN; % to avoid border effects

    [ny,nx,nz] = size(frame);
    
    [pstruct, mask] = pointSourceDetection3D(frame, sigma(mCh,:), 'Alpha', opts.Alpha,...
        'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant,...
        'RefineMaskLoG', false, 'WindowSize', opts.WindowSize); %#ok<PFBNS>
    
    if ~isempty(pstruct)
        pstruct.s = sigma;
        pstruct = rmfield(pstruct, 's_pstd');
        
        pstruct.dRange{mCh} = [min(frame(:)) max(frame(:))];
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
        
        % retain only mask regions containing localizations
        CC = bwconncomp(mask);
        labels = labelmatrix(CC);
        loclabels = labels(sub2ind([ny nx nz], pstruct.y_init, pstruct.x_init, pstruct.z_init));
        idx = setdiff(1:CC.NumObjects, loclabels);
        CC.PixelIdxList(idx) = [];
        CC.NumObjects = length(CC.PixelIdxList);
        
        % clean mask
        labels = labelmatrix(CC);
        mask = labels~=0;
        
        for ci = setdiff(1:nCh, mCh)
            frame = double(readtiff(data.framePathsDS{ci}{k}));
            pstruct.dRange{ci} = [min(frame(:)) max(frame(:))];
            X = [pstruct.x(mCh,:)' pstruct.y(mCh,:)' pstruct.z(mCh,:)'];
            
            [A_est, c_est] = estGaussianAmplitude3D(frame, sigma(ci,:));
            % linear index of positions
            linIdx = sub2ind(size(frame), roundConstr(X(:,2),ny), roundConstr(X(:,1),nx), roundConstr(X(:,3),nz));
            pstructSlave = fitGaussians3D(frame, X, A_est(linIdx), sigma(ci,:), c_est(linIdx), 'Ac');
            
            % localize, and compare intensities & (x,y)-coordinates. Use localization result if it yields better contrast
            pstructSlaveLoc = fitGaussians3D(frame, X, pstructSlave.A, sigma(ci,:), pstructSlave.c, 'xyAc');
            idx = sqrt((pstruct.x(mCh,:)-pstructSlaveLoc.x).^2 ...
                + (pstruct.y(mCh,:)-pstructSlaveLoc.y).^2) < 3*sigma(mCh,1) & pstructSlaveLoc.A > pstructSlave.A;
            
            % fill slave channel information
            for f = 1:length(sfields)
                pstruct.(sfields{f})(ci,~idx) = pstructSlave.(sfields{f})(~idx);
                pstruct.(sfields{f})(ci,idx) = pstructSlaveLoc.(sfields{f})(idx);
            end
            
            nanIdx = isnan(pstructSlave.x); % points within slave channel border, remove from detection results
            for f = 1:length(rmfields)
                pstruct.(rmfields{f})(:,nanIdx) = [];
            end
            
            pstruct.isPSF(ci,:) = ~pstruct.hval_AD(ci,:);
        end
        
        % add fields for tracker
        pstruct.xCoord = [pstruct.x(mCh,:)' pstruct.x_pstd(mCh,:)'];
        pstruct.yCoord = [pstruct.y(mCh,:)' pstruct.y_pstd(mCh,:)'];
        pstruct.zCoord = [pstruct.z(mCh,:)' pstruct.z_pstd(mCh,:)'] * data.zAniso;
        pstruct.amp =    [pstruct.A(mCh,:)' pstruct.A_pstd(mCh,:)'];
        frameInfo(k) = orderfields(pstruct, fieldnames(frameInfo(k))); %#ok<PFOUS>
    else
        frameInfo(k).dRange{mCh} = [min(frame(:)) max(frame(:))];
        for ci = setdiff(1:nCh, mCh)
            frame = double(imread(data.framePathsDS{ci}{k}));
            frameInfo(k).dRange{ci} = [min(frame(:)) max(frame(:))];
        end
        frameInfo(k).s = sigma;
    end
    
    maskPath = [resultsPath filesep 'Masks' filesep 'dmask_' num2str(k, fmt) '.tif'];
    writetiff(uint8(255*mask), maskPath);
end

save([resultsPath opts.FileName], 'frameInfo');


function y = roundConstr(x, nx)
y = round(x);
y(y<1) = 1;
y(y>nx) = nx;
