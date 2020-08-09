%[lftRes] = runLifetimeAnalysis(data, varargin) computes the CCP lifetime distribution
%
% Inputs:
%               data : structure returned by loadConditionData()
%
% Options ('specifier', value):
%       'SlaveNames' : cell array of string(s) for plot legends
%   'RemoveOutliers' : {true}|false toggles identification of data sets with outlier distributions
%
% Outputs:
%             lftRes : structure containing lifetime distributions

% Francois Aguet (last mod. 05/29/2013)

function [lftRes, res] = runLifetimeAnalysis(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) && numel(unique([data.framerate]))==1);
ip.addOptional('CohortBounds', [5 10 15 20 40 60 120]); % pairs: (5 10], (10 15],...
ip.addParamValue('Display', 'on', @(x) any(strcmpi(x, {'on', 'off', 'all'})));
ip.addParamValue('DisplayMode', 'screen', @(x) any(strcmpi(x, {'print', 'screen'})));
ip.addParamValue('ProcessedTracks', 'ProcessedTracks.mat', @ischar);
ip.addParamValue('LifetimeData', 'lifetimeData.mat', @ischar);
ip.addParamValue('Type', 'all', @ischar);
ip.addParamValue('Cutoff_f', 5, @isscalar);
ip.addParamValue('Print', false, @islogical);
ip.addParamValue('Buffer', 5);
ip.addParamValue('MaxIntensityThreshold', []);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('ClassificationSelector', 'significantMaster');
ip.addParamValue('ShowThresholdRange', false, @islogical);
ip.addParamValue('MaxP', 3);
ip.addParamValue('YLim', []);
ip.addParamValue('Rescale', true, @islogical);
ip.addParamValue('RemoveOutliers', true, @islogical);
ip.addParamValue('ExcludeVisitors', true, @islogical);
ip.addParamValue('FirstNFrames', []);
ip.addParamValue('PoolDatasets', false, @islogical);
ip.addParamValue('ShowStatistics', false, @islogical);
ip.addParamValue('SelectIndex', [], @iscell);
ip.addParamValue('SlaveNames', []);
ip.addParamValue('Colormap', jet(numel(data)));
ip.addParamValue('CorrectObservationBias', true, @islogical);
ip.addParamValue('InitDensity', 'mean', @(x) any(strcmpi(x, {'mean', 'median'})));
ip.addParamValue('AmplitudeCorrection', []);
ip.addParamValue('MaskPath', ['Detection' filesep 'cellmask.tif'], @ischar);
ip.addParamValue('Binning', 1, @isposint);
ip.addParamValue('CellArea', [], @isnumeric);
ip.parse(data, varargin{:});
cb = ip.Results.CohortBounds;
nc = numel(cb)-1; % # cohorts
mCh = find(strcmp(data(1).source, data(1).channels));
FirstNFrames = ip.Results.FirstNFrames;
selIdx = ip.Results.SelectIndex;
slaveNames = ip.Results.SlaveNames;
if isempty(slaveNames)
    slaveNames = data(1).markers(setdiff(1:numel(data(1).channels),mCh));
end
cmap = ip.Results.Colormap;

% median absolute deviation -> standard deviation
madFactor = 1/norminv(0.75, 0, 1);

% Extend all to max. movie length, in case of mismatch
buffer = ip.Results.Buffer;
cutoff_f = ip.Results.Cutoff_f;

Nmax = max([data.movieLength])-2*buffer;

% movieLength = min([data.movieLength]);
framerate = data(1).framerate;

firstN = 3:20;

% loop through data sets, load tracks, store max. intensities and lifetimes
res = struct([]);

[lftData, outlierIdx] = getLifetimeData(data, 'Overwrite', ip.Results.Overwrite,...
    'ReturnValidOnly', false, 'ExcludeVisitors', ip.Results.ExcludeVisitors, 'Cutoff_f', cutoff_f,...
    'Scale', ip.Results.Rescale, 'DisplayScaling', any(strcmpi(ip.Results.Display, {'on','all'})),...
    'RemoveOutliers', ip.Results.RemoveOutliers, 'Mask', true, 'Colormap', ip.Results.Colormap, ...
    'ProcessedTracks', ip.Results.ProcessedTracks, 'LifetimeData', ip.Results.LifetimeData,...
    'AmplitudeCorrectionFactor', ip.Results.AmplitudeCorrection);
if ~isempty(selIdx)
    selIdx(outlierIdx) = [];
end
cmap(outlierIdx,:) = [];
data(outlierIdx) = [];

if ip.Results.PoolDatasets
    pnames = {'lifetime_s', 'start', 'catIdx', 'A', 'lifetime_s_all', 'start_all', 'catIdx_all'};
    if isfield(lftData, 'significantMaster')
        pnames = [pnames 'significantMaster'];
    end
    for k = 1:numel(pnames)
        tmp.(pnames{k}) = vertcat(lftData.(pnames{k}));
    end
    tmp2 = arrayfun(@(i) i.visitors.lifetime_s, lftData, 'unif', 0);
    tmp.visitors.lifetime_s = vertcat(tmp2{:});
    tmp.a = 1;
    lftData = tmp;
end
fprintf('=============================================================================\n');
fprintf('Lifetime analysis - processing:   0%%');
nd = numel(lftData);

bf = ip.Results.Binning; % binning factor

lftRes.t = ((cutoff_f+(bf-1)/2):bf:Nmax)*framerate;
lftRes.cellArea = zeros(nd,1);
for i = 1:nd
    
    % Category statistics
    idx_Ia = [lftData(i).catIdx_all]==1;
    idx_Ib = [lftData(i).catIdx_all]==2;
    idx_IIa = [lftData(i).catIdx_all]==5;
    
    % raw histograms, calculated over N usable frames
    N = data(i).movieLength-2*buffer;
    t = ((cutoff_f+(bf-1)/2):bf:N)*framerate;
    
    % apply correction
    % relative probabilities:
    % P(obs. lifetime==1) = N
    % P(obs. lifetime==N) = 1
    % => weighting:
    if ip.Results.CorrectObservationBias
        w = N./(N-cutoff_f+1:-1:1);
    else
        w = ones(1,N-cutoff_f+1);
    end
    if bf>1
        w = mean(reshape(w(1:floor(numel(w)/bf)*bf), [bf floor(numel(w)/bf)]),1);
    end
    pad0 = zeros(1,(Nmax-N)/bf);
    lftHist_Ia =  [hist(lftData(i).lifetime_s_all(idx_Ia), t).*w  pad0];
    lftHist_Ib =  [hist(lftData(i).lifetime_s_all(idx_Ib), t).*w  pad0];
    lftHist_IIa = [hist(lftData(i).lifetime_s_all(idx_IIa), t).*w pad0];
    
    % Normalization
    lftRes.lftHist_Ia(i,:) = lftHist_Ia / sum(lftHist_Ia) / framerate / bf;
    lftRes.lftHist_Ib(i,:) = lftHist_Ib / sum(lftHist_Ib) / framerate / bf;
    lftRes.lftHist_IIa(i,:) = lftHist_IIa / sum(lftHist_IIa) / framerate / bf;
    lftRes.nSamples_Ia(i) = sum(idx_Ia);
    
    %-------------------------------------------------------------
    % Max. intensity distribution for cat. Ia CCP tracks
    %-------------------------------------------------------------
    for k = 1:nc
        % indexes within cohorts
        cidx = cb(k)<lftData(i).lifetime_s & lftData(i).lifetime_s<=cb(k+1);
        res(i).maxA{k} = nanmax(lftData(i).A(cidx,:,mCh),[],2);
        
        % lifetimes for given cohort
        res(i).lft{k} = lftData(i).lifetime_s(cidx);
    end
    
    res(i).maxA_all = nanmax(lftData(i).A(:,:,mCh),[],2);
    if isfield(lftData, 'significantMaster')
        res(i).significantMaster = lftData(i).significantMaster;
    end
    res(i).firstN = firstN;
    
    fprintf('\b\b\b\b%3d%%', round(100*i/nd));
end
fprintf('\n');


%====================
% Threshold
%====================
% Check whether saved and whether equal
% hasThreshold = false(1,nd);
% tvec = NaN(1,nd);
% for i = 1:nd
%     prmFile = [data(i).source 'Analysis' filesep 'Parameters.mat'];
%     if exist(prmFile, 'file')==2
%         hasThreshold(i) = true;
%         prm = load(prmFile);
%         tvec(i) = prm.MaxIntensityThreshold;
%     end
% end
% tvec = unique(tvec);

% if all data sets have the same stored threshold value, and no values are given in
% input, use stored threshold
% if numel(tvec)==1 && ~isnan(tvec) && isempty(ip.Results.MaxIntensityThreshold) && isempty(FirstNFrames)
%     T = tvec(1);
%     fprintf('Max. intensity threshold on first %d frames: %.2f\n', prm.nFramesThreshold, T);
% else
if ~isempty(ip.Results.MaxIntensityThreshold)
    T = ip.Results.MaxIntensityThreshold;
else
    
    A = arrayfun(@(i) i.A(:,:,mCh), lftData, 'unif', 0);
    A = vertcat(A{:});
    lft = vertcat(lftData.lifetime_s);
    
    % check whether all cohorts contain tracks
    ntc = arrayfun(@(c) sum(cb(c)<lft & lft<=cb(c+1)), 1:nc);
    if any(ntc<5)
        warning('At least one lifetime cohorts contains less than 5 tracks.');
        fprintf('Calculating cohort bounds based on lifetime distribution percentiles\n.');
        cb = prctile(lft, linspace(0,100,7));
        %cb = linspace(min(lft), max(lft), nc+1);
        cb(1) = cb(1)-framerate;
    end
    
    if isempty(FirstNFrames)
        frameRange = 3:12;
        hval = zeros(1,frameRange(end));
        for ni = 1:numel(frameRange)
            M = max(A(:,1:frameRange(ni)),[],2);
            
            muC = zeros(1,nc);
            sC = zeros(1,nc);
            for c = 1:nc
                cidx = cb(c)<lft & lft<=cb(c+1);
                [muC(c), sC(c)] = fitGaussianModeToPDF(M(cidx,:), 'FixMode', true);
            end
            hval(frameRange(ni)) = adtest1(muC(2:end), 'mu', muC(1), 'sigma', sC(1)/sqrt(nc));
        end
        FirstNFrames = find(hval==1, 1, 'first')-1;
    end
    
    M = nanmax(A(:,1:FirstNFrames,mCh),[],2);
    
    [mu_g, sigma_g] = fitGaussianModeToPDF(M, 'FixMode', true, 'Display', false);
    T = norminv(0.99, mu_g, sigma_g);
    
    % save threshold value
    prm.MaxIntensityThreshold = T;
    prm.nFramesThreshold = FirstNFrames;
    for i = 1:nd
        if isfield(lftData(i), 'a')
            prm.a = lftData(i).a;
        else
            prm.a = ones(numel(data(i).channels),1);
        end
        save([data(i).source 'Analysis' filesep 'Parameters.mat'], '-struct', 'prm');
    end
    
    % 95th percentile of first frame intensity distribution
    T95 = prctile(A(:,1,mCh), 95);
    
    fprintf('Max. intensity threshold on first %d frames: %.2f\n', FirstNFrames, T);
    fprintf('95th percentile of 1st frame distribution: %.2f\n', T95);
end
lftRes.MaxIntensityThreshold = T;

% loop through data sets, apply max. intensity threshold
lftRes.pctCCP = zeros(nd,1);
lftRes.pctCS = zeros(nd,1);
lftRes.pctVisit = zeros(nd,1);
minMovieLength = min([data.movieLength]);
for i = 1:nd
    
    % Selection indexes for each data set
    % 1) Intensity threshold based on maximum intensity distribution
    idxMI = res(i).maxA_all >= T;
    
    if ~isempty(selIdx)
        idxMI = idxMI & selIdx{i};
    end
    
    % 2) Remove non-endocytic structures (i.e., endosomal CCSs)
    if ip.Results.ExcludeVisitors
        res(i).lftVisitors = lftData(i).visitors.lifetime_s;
        nCS = numel(lftData(i).lifetime_s) + numel(lftData(i).visitors.lifetime_s);
        lftRes.pctVisit(i) = numel(res(i).lftVisitors) / nCS;
    else
        nCS = numel(lftData(i).lifetime_s);
    end
    res(i).lftAboveT = lftData(i).lifetime_s(idxMI);
    res(i).lftBelowT = lftData(i).lifetime_s(~idxMI);
    lftRes.nCCP(i) = sum(idxMI);
    
    %res(i).maxAAboveT = res(i).maxA_all(idxMI);
    %res(i).AaboveT = lftData(i).A(idxMI,:,:);
    lftRes.pctCCP(i) = sum(idxMI)/nCS;
    lftRes.pctCS(i) = sum(~idxMI)/nCS;
    
    N = data(i).movieLength-2*buffer;
    t = ((cutoff_f+(bf-1)/2):bf:N)*framerate;
    
    % apply correction
    % relative probabilities:
    % P(obs. lifetime==1) = N
    % P(obs. lifetime==N) = 1
    % => weighting:
    if ip.Results.CorrectObservationBias
        w = N./(N-cutoff_f+1:-1:1);
    else
        w = ones(1,N-cutoff_f+1);
    end
    if bf>1
        w = mean(reshape(w(1:floor(numel(w)/bf)*bf), [bf floor(numel(w)/bf)]),1);
    end
    pad0 = zeros(1,(Nmax-N)/bf);
    lftHistCCP = [hist(res(i).lftAboveT, t).*w pad0];
    lftHistCS = [hist(res(i).lftBelowT, t).*w pad0];
    % Normalization
    lftRes.lftHistCCP(i,:) = lftHistCCP / sum(lftHistCCP) / framerate / bf;
    lftRes.lftHistCS(i,:) = lftHistCS / sum(lftHistCS) / framerate / bf;
    
    % Raw, unweighted histograms with counts/bin
    lftRes.lftHistCCP_counts(i,:) = [hist(res(i).lftAboveT, t) pad0];
    lftRes.lftHistCS_counts(i,:) = [hist(res(i).lftBelowT, t) pad0];
    lftRes.histWeight = w;
    
    if ip.Results.ExcludeVisitors
        lftHistVisit = [hist(res(i).lftVisitors, t).*w pad0];
        lftRes.lftHistVisit(i,:) = lftHistVisit / sum(lftHistVisit) / framerate / bf;
    end
    
    % Multi-channel data
    if isfield(res, 'significantMaster')
        
        % slave channels
        sCh = setdiff(1:numel(data(i).channels), mCh);
        ns = numel(sCh);
        
        % slave combinations
        scomb = dec2bin(2^ns-1:-1:0)=='1';
        lftRes.slaveCombs = scomb;
        
        for s = 1:size(scomb,1)
            sIdx = all(bsxfun(@eq, lftData(i).significantMaster(:,sCh), scomb(s,:)),2);
            lftHistSlaveAll = [hist(lftData(i).lifetime_s(sIdx), t).*w pad0];
            lftHistSlaveCCP = [hist(lftData(i).lifetime_s(sIdx & idxMI), t).*w pad0];
            lftHistSlaveCS = [hist(lftData(i).lifetime_s(sIdx & ~idxMI), t).*w pad0];
            if ~(sum(lftHistSlaveAll)==0)
                lftRes.lftHistSlaveAll{s}(i,:) = lftHistSlaveAll / sum(lftHistSlaveAll) / framerate / bf;
            else
                lftRes.lftHistSlaveAll{s}(i,:) = zeros(size(lftHistSlaveAll));
            end
            if ~(sum(lftHistSlaveCCP)==0)
                lftRes.lftHistSlaveCCP{s}(i,:) = lftHistSlaveCCP / sum(lftHistSlaveCCP) / framerate / bf;
            else
                lftRes.lftHistSlaveCCP{s}(i,:) = zeros(size(lftHistSlaveCCP));
            end
            lftRes.lftHistSlaveCS{s}(i,:) = lftHistSlaveCS / sum(sIdx & ~idxMI) / framerate / bf;
            lftRes.lftHistSlaveAll_counts{s}(i,:) = [hist(lftData(i).lifetime_s(sIdx), t) pad0];
            lftRes.lftHistSlaveCCP_counts{s}(i,:) = [hist(lftData(i).lifetime_s(sIdx & idxMI), t) pad0];
            lftRes.lftHistSlaveCS_counts{s}(i,:) = [hist(lftData(i).lifetime_s(sIdx & ~idxMI), t) pad0];
            lftRes.pctSlaveCCP(i,s) = sum(idxMI & sIdx)/numel(idxMI);
            lftRes.pctSlaveCS(i,s) = sum(~idxMI & sIdx)/numel(idxMI);
        end
    end
    % Note: sum(pctSlaveCCP,2) differs from pctCCP because it excludes visitors
    
    %-----------------------------------
    % Initiation density
    %-----------------------------------
    % Cell area
    if ~isempty(ip.Results.CellArea)
        lftRes.cellArea(i) = ip.Results.CellArea(i);
    elseif exist([data(i).source ip.Results.MaskPath], 'file')==2
        px = data(i).pixelSize / data(i).M; % pixels size in object space
        mask = logical(readtiff([data(i).source ip.Results.MaskPath]));
        lftRes.cellArea(i) = sum(mask(:)) * px^2 * 1e12; % in µm^2
    else
        lftRes.cellArea(i) = NaN;
    end
    
    % birth/death statistics
    startsPerFrameAll = hist(lftData(i).start_all, 1:minMovieLength);
    startsPerFrameAll = startsPerFrameAll(6:end-2);
    startsPerFrameIa = hist(lftData(i).start(lftData(i).catIdx==1), 1:minMovieLength);
    startsPerFrameIa = startsPerFrameIa(6:end-2);
    startsPerFrameCCP = hist(lftData(i).start(idxMI), 1:minMovieLength);
    startsPerFrameCCP = startsPerFrameCCP(6:end-2);
    lftRes.startsPerFrameAll(i,:) = startsPerFrameAll;
    
    % in µm^-2 min^-1
    dnorm = 60/data(i).framerate/lftRes.cellArea(i);
    if strcmpi(ip.Results.InitDensity, 'mean')
        lftRes.initDensityAll(i,:) = [mean(startsPerFrameAll); std(startsPerFrameAll)]*dnorm;
        lftRes.initDensityIa(i,:) =  [mean(startsPerFrameIa);  std(startsPerFrameIa)]*dnorm;
        lftRes.initDensityCCP(i,:) = [mean(startsPerFrameCCP); std(startsPerFrameCCP)]*dnorm;
    else
        lftRes.initDensityAll(i,:) = [median(startsPerFrameAll); madFactor*mad(startsPerFrameAll, 1)]*dnorm;
        lftRes.initDensityIa(i,:) =  [median(startsPerFrameIa);  madFactor*mad(startsPerFrameIa, 1)]*dnorm;
        lftRes.initDensityCCP(i,:) = [median(startsPerFrameCCP); madFactor*mad(startsPerFrameCCP,1)]*dnorm;
    end
    lftRes.persistentDensity(i,:) = sum(lftData(i).catIdx_all==4) / lftRes.cellArea(i);
end

% print lifetime distribution percentiles
lftCDF = cumsum(mean(lftRes.lftHistCS,1));
[~,uidx] = unique(lftCDF);
lftPctCS = interp1(lftCDF(uidx), lftRes.t(uidx), [0.05 0.25 0.5 0.75 0.95]);
lftCDF = cumsum(mean(lftRes.lftHistCCP,1));
[~,uidx] = unique(lftCDF);
lftPctCCP = interp1(lftCDF(uidx), lftRes.t(uidx), [0.05 0.25 0.5 0.75 0.95]);
fprintf('Lifetime distribution percentiles (5th, 25th, 50th, 75th, 95th):\n');
fprintf(['  CSs:  [' strjoin(arrayfun(@(i) num2str(i, '%.1f'), lftPctCS, 'unif', 0), ', ') '] s\n']);
fprintf(['  CCPs: [' strjoin(arrayfun(@(i) num2str(i, '%.1f'), lftPctCCP, 'unif', 0), ', ') '] s\n']);

%====================
% Initiation density
%====================
if any(strcmpi(ip.Results.Display, {'on','all'})) && ~ip.Results.PoolDatasets
    
    ha = setupFigure(2,2, 'SameAxes', false, 'AxesWidth', 10, 'AxesHeight', 7.5,...
        'XSpace', [3 3 0.5], 'YSpace', [4 1 1], 'Name', 'Track statistics');
    fset = loadFigureSettings('');
    set(ha, 'FontSize', 12);
    XTickLabel = arrayfun(@getMovieName, data, 'unif', 0);
    % 1) Density
    if ~isnan(lftRes.cellArea)
        hb = barplot2([lftRes.initDensityAll(:,1) lftRes.initDensityIa(:,1) lftRes.initDensityCCP(:,1)],...
            [lftRes.initDensityAll(:,2) lftRes.initDensityIa(:,2) lftRes.initDensityCCP(:,2)], [],[],...
            'XTickLabel', XTickLabel, 'Interpreter', 'none',...
            'FaceColor', hsv2rgb([0 0 0.5; 0.33 0.5 0.5; 0.33 0.8 1]), 'Handle', ha(1), 'AdjustFigure', false);
        ylabel(ha(1), ['Initiations (' char(181) 'm^{-2} min^{-1})'], fset.lfont{:});
        hl = legend(ha(1), hb, ' All tracks', ' Valid tracks', 'CCPs');
        set(hl, fset.tfont{:});
        
        % 2) Cell area
        colormap(cmap);
        %scatter(ha(2), 1:nd, lftRes.cellArea, 50, cmap, 'o', 'fill', 'MarkerEdgeColor', 'k');
        scatter(ha(2), 1:nd, lftRes.cellArea, 50, 1:size(cmap,1), 'o', 'fill', 'MarkerEdgeColor', 'k');
        ylabel(ha(2), ['Cell area (' char(181) 'm^2)'], fset.lfont{:});
        YLim = get(ha(2), 'YLim');
        YLim(1) = 0;
        set(ha(2), 'XTick', 1:nd, 'XTickLabel', XTickLabel, 'XLim', [0.5 nd+0.5], 'YLim', YLim);
        rotateXTickLabels(ha(2), 'Angle', 45, 'AdjustFigure', false, 'Interpreter', 'none');
    end
    
    
    % # valid tracks
    scatter(ha(3), 1:nd, lftRes.nSamples_Ia, 50, 1:size(cmap,1), 'o', 'fill', 'MarkerEdgeColor', 'k');
    ylabel(ha(3), '# valid tracks', fset.lfont{:});
    YLim = get(ha(3), 'YLim');
    YLim(1) = 0;
    set(ha(3), 'XTick', 1:nd, 'XTickLabel', XTickLabel, 'XLim', [0.5 nd+0.5], 'YLim', YLim);
    rotateXTickLabels(ha(3), 'Angle', 45, 'AdjustFigure', false, 'Interpreter', 'none');
    
    % # persistent tracks
    scatter(ha(4), 1:nd, arrayfun(@(i) sum([i.catIdx_all]==4), lftData), 50, fset.ceTrackClasses(4,:), 'o', 'fill', 'MarkerEdgeColor', 'k');
    scatter(ha(4), 1:nd, arrayfun(@(i) sum([i.catIdx_all]==8), lftData), 50, fset.ceTrackClasses(8,:), 'o', 'fill', 'MarkerEdgeColor', 'k');
    ylabel(ha(4), '# persistent tracks', fset.lfont{:});
    YLim = get(ha(4), 'YLim');
    YLim(1) = 0;
    set(ha(4), 'XTick', 1:nd, 'XTickLabel', XTickLabel, 'XLim', [0.5 nd+0.5], 'YLim', YLim);
    hl = legend(ha(4), 'Single tracks', 'Compound tracks', 'Location', 'NorthWest');
    set(hl, 'Box', 'on');
    rotateXTickLabels(ha(4), 'Angle', 45, 'AdjustFigure', false, 'Interpreter', 'none');
    
    formatTickLabels(ha(1:2));
    
    fprintf('Initiation density, all detected tracks:                  %.3f ± %.3f [µm^-2 min^-1]\n', mean(lftRes.initDensityAll(:,1)), std(lftRes.initDensityAll(:,1)));
    fprintf('Initiation density, valid tracks (CCPs + CSs + visitors): %.3f ± %.3f [µm^-2 min^-1]\n', mean(lftRes.initDensityIa(:,1)), std(lftRes.initDensityIa(:,1)));
    fprintf('Initiation density, CCPs:                                 %.3f ± %.3f [µm^-2 min^-1]\n', mean(lftRes.initDensityCCP(:,1)), std(lftRes.initDensityCCP(:,1)));
    fprintf('Density of persistent structures:                         %.3f ± %.3f [µm^-2]\n', mean(lftRes.persistentDensity), std(lftRes.persistentDensity));
    fprintf('  Valid tracks/cell: %.1f ± %.1f (total valid tracks: %d)\n', mean(lftRes.nSamples_Ia), std(lftRes.nSamples_Ia), sum(lftRes.nSamples_Ia));
    
    % gap statistics
    ha = setupFigure(1,2, 'SameAxes', false, 'AxesWidth', 10, 'AxesHeight', 7.5,...
        'XSpace', [3 3 0.5], 'YSpace', [2 1 0.5], 'Name', 'Gap statistics');
    set(ha, 'FontSize', 12);
    
    % 1) #gaps/track
    ngaps = cell(1,nd);
    for k = 1:nd
        ngaps{k} = sum(lftData(k).gapMat_Ia, 2);
    end
    maxGaps = max(vertcat(ngaps{:}));
    xi = 1:maxGaps;
    for k = 1:nd
        plot(ha(1), xi, hist(ngaps{k}, xi)/numel(ngaps{k})*100, 'Color', cmap(k,:), 'LineWidth', 1);
    end
    hl = legend(ha(1), XTickLabel, 'Interpreter', 'none', 'Location', 'NorthEast');
    set(hl, 'Box', 'off');
    xlabel(ha(1), '# gaps', fset.lfont{:});
    ylabel(ha(1), '% of tracks', fset.lfont{:});
    
    % 2) gap length
    gapLengths = cell(1,nd);
    for k = 1:nd
        cg = cumsum(lftData(k).gapMat_Ia, 2);
        dx = [diff(lftData(k).gapMat_Ia, [], 2)==-1 zeros(size(lftData(k).gapMat_Ia,1),1)];
        M = cg .* dx;
        tmp = arrayfun(@(i) M(i,M(i,:)~=0), 1:size(M,1), 'unif', 0);
        tmp = cellfun(@(i) diff([0 i]), tmp, 'unif', 0);
        gapLengths{k} = [tmp{:}];
    end
    maxg = max([gapLengths{:}]);
    xi = 1:maxg+1;
    for k = 1:nd
        plot(ha(2), xi, hist(gapLengths{k}, xi)/numel(gapLengths{k})*100, '.-',...
            'Color', cmap(k,:), 'LineWidth', 1, 'MarkerSize', 12);
    end
    set(ha(2), 'XLim', [0.5 maxg+1.5], 'XTick', 1:maxg+1);
    xlabel(ha(2), 'Gap length (frames)', fset.lfont{:});
    ylabel(ha(2), '% of gaps', fset.lfont{:});
    
    % 3) #gaps as % of track vs. lifetime
    %for k = 1:nd
    %    ngaps{k} = ngaps{k} ./ lftData(k).trackLengths;
    %    scatter(ha(2), lftData(k).lifetime_s, ngaps{k}, 50, cmap(k,:), 'o', 'fill', 'MarkerEdgeColor', 'k');
    %end
    %xi = linspace(min(vertcat(ngaps{:})), max(vertcat(ngaps{:})), 10);
    %for k = 1:nd
    %    plot(ha(2), xi, hist(ngaps{k}, xi)/numel(ngaps{k}), 'Color', cmap(k,:));
    %end
    
    
    
    % plot cumulative lifetime distributions
    t = lftRes.t;
    med = zeros(nd,1);
    
    ha = setupFigure(1,2, 'SameAxes', true, 'AxesWidth', 10, 'AxesHeight', 7.5,...
        'XSpace', [2 1 1], 'YSpace', [2 1 1], 'Name', 'Cumulative lifetime distributions');
    set(ha, 'FontSize', 12);
    hp = zeros(nd,1);
    
    % all structures
    edf = cumsum(lftRes.lftHist_Ia,2)*framerate*bf;
    % plot(ha(1), [0 200], 0.75*[1 1], 'k--');
    for k = 1:nd
        hp(k) = plot(ha(1), lftRes.t, edf(k,:), 'Color', cmap(k,:));
        [~,idx] = unique(edf(k,:));
        med(k) = interp1(edf(k,idx), t(idx), 0.75);
    end
    [~,idx] = sort(med);
    hl = legend(hp(idx), XTickLabel(idx), 'Interpreter', 'none', 'Location', 'SouthEast');
    set(hl, 'Box', 'off');
    
    % CCPs
    edf = cumsum(lftRes.lftHistCCP,2)*framerate*bf;
    for k = 1:nd
        hp(k) = plot(ha(2), lftRes.t, edf(k,:), 'Color', cmap(k,:));
        [~,idx] = unique(edf(k,:));
        med(k) = interp1(edf(k,idx), t(idx), 0.75);
    end
    [~,idx] = sort(med);
    hl = legend(hp(idx), XTickLabel(idx), 'Interpreter', 'none', 'Location', 'SouthEast');
    set(hl, 'Box', 'off');
    
    axis(ha, [0 200 0 1]);
    formatTickLabels(ha);
    ylabel(ha(1), 'Cumulative frequency', fset.lfont{:});
    xlabel(ha(1), 'Lifetime (s)', fset.lfont{:});
    xlabel(ha(2), 'Lifetime (s)', fset.lfont{:});
    title(ha(1), 'All valid trajectories', fset.lfont{:});
    title(ha(2), 'CCPs', fset.lfont{:});
end

% Mean distributions
lftRes.meanLftHistCCP = nanmean(lftRes.lftHistCCP,1);
lftRes.meanLftHistCS = nanmean(lftRes.lftHistCS,1);
if ip.Results.ExcludeVisitors
    lftRes.meanLftHistVisit = mean(lftRes.lftHistVisit,1);
end


%---------------------------------------
% Raw lifetime distributions + average
%---------------------------------------
if any(strcmpi(ip.Results.Display, {'all'}))
    
    a = [lftData.a];
    colorV = hsv(nd);
    [~,idxa] = sort(a(1,:));
    [~,idxa] = sort(idxa);
    colorV = colorV(idxa,:);
    
    fset = loadFigureSettings('');
    setupFigure('DisplayMode', 'screen', 'Name', 'Raw lifetime distribution');
    hold on;
    for i = nd:-1:1
        plot(lftRes.t, lftRes.lftHist_Ia(i,:), '-', 'Color', colorV(i,:), 'LineWidth', 1);
    end
    plot(lftRes.t, mean(vertcat(lftRes.lftHist_Ia), 1), 'k', 'LineWidth', 2);
    ya = 0:0.02:0.1;
    axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
    set(gca, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'unif', 0)]);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    
    % Inset with zoom
    zf = 0.6;
    aw = fset.axPos(3);
    ah = fset.axPos(4);
    axes(fset.axOpts{:}, 'Units', fset.units, 'Position', [fset.axPos(1)+(1-zf)*aw fset.axPos(2)+(1-zf)*ah zf*aw zf*ah]);
    hold on;
    for i = nd:-1:1
        plot(lftRes.t, lftRes.lftHist_Ia(i,:), '-', 'Color', colorV(i,:), 'LineWidth', 1);
    end
    plot(lftRes.t, mean(vertcat(lftRes.lftHist_Ia), 1), 'k', 'LineWidth', 2);
    axis([0 60 0 0.035]);
    ya = 0:0.01:0.04;
    set(gca, 'FontSize', 7, 'TickLength', fset.TickLength/zf, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'unif', 0)]);
    
    
    lftCDF = cumsum(mean(vertcat(lftRes.lftHist_Ia),1))*framerate;
    [uCDF, idx] = unique(lftCDF);
    lft50 = interp1(uCDF, lftRes.t(idx), 0.5);
    
    figure(fset.fOpts{:}, 'Name', 'Raw lifetime distribution');
    axes(fset.axOpts{:});
    hold on;
    meanHist = mean(vertcat(lftRes.lftHist_Ia), 1);
    plot(lftRes.t, meanHist, 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    ya = 0:0.02:0.1;
    axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
    set(gca, fset.axOpts{:}, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'unif', 0)]);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    [mu,~,Aexp] = fitExpToHist(lftRes.t, meanHist);
    plot(lftRes.t, Aexp/mu*exp(-1/mu*lftRes.t), 'r-', 'LineWidth', 1);
    %hl = legend(' Best exponential fit', 'Location', 'SouthEast');
    %set(hl, 'Box', 'off', 'Position', [4.5 1.5 2 1]);
    
    % Inset with zoom
    zf = 0.6;
    aw = fset.axPos(3);
    ah = fset.axPos(4);
    axes(fset.axOpts{:}, 'Units', fset.units, 'Position', [fset.axPos(1)+(1-zf)*aw fset.axPos(2)+(1-zf)*ah zf*aw zf*ah]);
    hold on;
    idx = find(lftRes.t==round(lft50/framerate)*framerate);
    fill([lftRes.t(1:idx) lftRes.t(idx:-1:1)], [lftCDF(1:idx) zeros(1,idx)], fset.ceB, 'EdgeColor', 'none');
    plot(lftRes.t, lftCDF, 'k', 'LineWidth', 1.5);
    plot([0 lft50], [0.5 0.5], 'k--', 'LineWidth', 1);
    ya = 0:0.25:1;
    set(gca, 'FontSize', 7, 'TickLength', fset.TickLength/zf, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'unif', 0)]);
    axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
    ylabel('Cumulative freq.', fset.sfont{:});
    
    plotMaxIntensityDistribution(data);
end


if any(strcmpi(ip.Results.Display, {'on','all'}))
    %printPath = [getExpDir(data) 'Figures' filesep];
    %[~,~] = mkdir(printPath);
    plotLifetimes(lftRes, 'ShowStatistics', ip.Results.ShowStatistics,...
        'DisplayMode', ip.Results.DisplayMode, 'PlotAll', true,...
        'SlaveNames', slaveNames, 'SingleChannel', numel(data(1).channels)==1);
    %print(h(1), '-depsc2', '-loose', [printPath 'lifetimeDistributions.eps']);
    %if ip.Results.ShowStatistics
    %    print(h(2), '-depsc2', '-loose', [printPath 'lifetimeDistributionsStats.eps']);
    %end
end

if ip.Results.Print
    fpath = cell(1,nd);
    for k = 1:nd
        [~,fpath{k}] = getCellDir(data(k));
    end
    fpath = unique(fpath);
    if numel(fpath)>1
        fprintf('Figures could not be printed.');
    else
        fpath = [fpath{1} 'Figures' filesep];
        [~,~] = mkdir(fpath);
        %print(hf(1), '-depsc2', [fpath 'trackClassDistribution.eps']);
        %print(hf(2), '-depsc2', [fpath 'meanLftHist_classes.eps']);
        %print(hf(3), '-depsc2', [fpath 'gapStatistics.eps']);
    end
end
