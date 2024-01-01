%[lftData, rmIdx] = getLifetimeData(data, varargin) returns track information in compact form for lifetime analysis

% Francois Aguet, 05/2012

function [lftData, rmIdx] = getLifetimeData(data, varargin)

nd = numel(data);
nCh = numel(data(1).channels);

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('ProcessedTracks', [], @ischar);
ip.addParamValue('LifetimeData', 'lifetimeData.mat', @ischar);
ip.addParamValue('ReturnValidOnly', true, @islogical);
ip.addParamValue('Cutoff_f', [], @isscalar);
ip.addParamValue('ExcludeVisitors', true, @islogical);
ip.addParamValue('Scale', false, @islogical);
ip.addParamValue('DisplayScaling', false, @islogical);
ip.addParamValue('RemoveOutliers', false, @islogical);
ip.addParamValue('AmplitudeCorrectionFactor', [], @(x) isempty(x) || (size(x,1)==nd && size(x,2)==nCh));
ip.addParamValue('Mask', false, @islogical);
ip.addParamValue('Colormap', []);
ip.addParamValue('AnalysisPath', 'Tracking', @ischar);
ip.parse(varargin{:});

rescale = ip.Results.Scale;
if numel(rescale)==1
    rescale = repmat(rescale, [nCh 1]);
end

fnames = {'lifetime_s', 'trackLengths', 'start', 'catIdx', 'index', 'A', 'A_pstd',...
    'sigma_r', 'SE_sigma_r', 'sbA', 'ebA', 'sbSigma_r', 'ebSigma_r', 'gapMat_Ia'};
lftData(1:nd) = cell2struct(cell(size(fnames)), fnames, 2);
vnames = fnames(1:5);
mnames = fnames(6:end);

parfor i = 1:nd
    fpath = [data(i).source 'Analysis' filesep ip.Results.LifetimeData]; %#ok<PFBNS>
    if ~(exist(fpath, 'file')==2) || ip.Results.Overwrite
        
        tracks = loadTracks(data(i), 'Mask', ip.Results.Mask, 'Category', 'all', 'Cutoff_f', 2,...
            'AnalysisPath', ip.Results.AnalysisPath, 'FileName', ip.Results.ProcessedTracks);
        
        % concatenate amplitudes of master channel into matrix
        trackLengths = [tracks.end]-[tracks.start]+1;
        
        lftData(i).lifetime_s = [tracks.lifetime_s]';
        lftData(i).trackLengths = trackLengths';
        lftData(i).start = [tracks.start]';
        lftData(i).catIdx = [tracks.catIdx]';
        lftData(i).index = (1:numel(tracks))';
        if isfield(tracks, 'significantMaster')
            lftData(i).significantMaster = [tracks.significantMaster]';
            lftData(i).significantSlave = [tracks.significantSlave]';
        end
        
        % store intensities of cat. Ia tracks
        idx_Ia = find([tracks.catIdx]==1);
        tracks = tracks(idx_Ia);
        
        nt = numel(tracks);
        if nt>0
            nsb = numel(tracks(1).startBuffer.t);
            neb = numel(tracks(1).endBuffer.t);
            
            % store intensity matrices
            nf = data(i).movieLength;
            lftData(i).A = NaN(nt,nf,nCh);
            lftData(i).A_pstd = NaN(nt,nf,nCh);
            lftData(i).sigma_r = NaN(nt,nf,nCh);
            lftData(i).SE_sigma_r = NaN(nt,nf,nCh);
            
            lftData(i).sbA = NaN(nt,nsb,nCh);
            lftData(i).ebA = NaN(nt,neb,nCh);
            lftData(i).sbSigma_r = NaN(nt,nsb,nCh);
            lftData(i).ebSigma_r = NaN(nt,neb,nCh);
            %lftData(i).RSS = NaN(nt,nf,nCh);
            lftData(i).gapMat_Ia = false(nt,nf);
            for k = 1:nt
                range = 1:trackLengths(idx_Ia(k));
                lftData(i).A(k,range,:) = tracks(k).A';
                lftData(i).A_pstd(k,range,:) = tracks(k).A_pstd';
                lftData(i).sigma_r(k,range,:) = tracks(k).sigma_r';
                lftData(i).SE_sigma_r(k,range,:) = tracks(k).SE_sigma_r';
                lftData(i).sbA(k,:,:) = tracks(k).startBuffer.A';
                lftData(i).ebA(k,:,:) = tracks(k).endBuffer.A';
                lftData(i).sbSigma_r(k,:,:) = tracks(k).startBuffer.sigma_r';
                lftData(i).ebSigma_r(k,:,:) = tracks(k).endBuffer.sigma_r';
                %lftData(i).RSS(k,range,:) = tracks(k).RSS';
                lftData(i).gapMat_Ia(k,range) = tracks(k).gapVect';
            end
        end
    else
        tmp = load(fpath);
        if isfield(tmp, 'significantMaster')
            lftData(i).significantMaster = [];
            lftData(i).significantSlave = [];
        end
        if isfield(tmp, 'RSS')
            tmp = rmfield(tmp, 'RSS');
        end
        if isfield(tmp, 'significantSignal') % deprecated
            tmp2 = tmp.significantSignal;
            tmp = rmfield(tmp, 'significantSignal');
            tmp.significantMaster = tmp2;
            tmp.significantSlave = NaN(size(tmp2));
            lftData(i).significantMaster = [];
            lftData(i).significantSlave = [];
        end
        if ~isfield(tmp, 'index');
            tmp.index = NaN(size(tmp.catIdx));
            if isfield(tmp, 'significantMaster')
                tmp = orderfields(tmp, [fnames, 'significantMaster', 'significantSlave']);
            else
                tmp = orderfields(tmp, fnames);
            end
        end
        lftData(i) = tmp;
    end
end

% amplitude fields
afields = {'A', 'A_pstd', 'sigma_r', 'SE_sigma_r', 'sbA', 'ebA', 'sbSigma_r', 'ebSigma_r'};

% apply amplitude correction
acorr = ip.Results.AmplitudeCorrectionFactor;
if ~isempty(acorr)
    for c = 1:nCh
        for i = 1:nd
            for f = 1:numel(afields)
                if ~isempty(lftData(i).A)
                    lftData(i).(afields{f})(:,:,c) = acorr(i,c)*lftData(i).(afields{f})(:,:,c);
                end
            end
        end
    end
end

% save outside of parfor
for i = 1:nd
    [~,~] = mkdir([data(i).source 'Analysis']);
    fpath = [data(i).source 'Analysis' filesep ip.Results.LifetimeData];
    if ~(exist(fpath, 'file')==2) || ip.Results.Overwrite
        iData = lftData(i); %#ok<NASGU>
        save(fpath, '-struct', 'iData');
    end
end

if isfield(lftData(1), 'significantMaster')
    vnames = [vnames 'significantMaster' 'significantSlave'];
    fnames = [vnames mnames];
end

maxA = cell(nCh,nd);
for i = 1:nd
    for c = 1:nCh
        if rescale(c) && size(lftData(i).A,1) > 0 % in case of empty input
            maxA{c,i} = nanmax(lftData(i).A(:,:,c),[],2);
        end
    end
    
    % apply frame cutoff to all fields
    if ~isempty(ip.Results.Cutoff_f)
        idx = lftData(i).trackLengths(lftData(i).catIdx==1)>=ip.Results.Cutoff_f;
        for f = 1:numel(mnames)
            lftData(i).(mnames{f}) = lftData(i).(mnames{f})(idx,:,:);
        end
        
        idx = lftData(i).trackLengths>=ip.Results.Cutoff_f;
        for f = 1:numel(vnames)
            lftData(i).(vnames{f}) = lftData(i).(vnames{f})(idx,:);
        end
    end
    
    if ~ip.Results.ReturnValidOnly
        for f = 1:numel(vnames)
            lftData(i).([vnames{f} '_all']) = lftData(i).(vnames{f});
        end
    end
    
    % remaining fields: retain category==1
    idx = lftData(i).catIdx==1;
    for f = 1:numel(vnames)
        lftData(i).(vnames{f}) = lftData(i).(vnames{f})(idx,:);
    end
    
    % remove visitors
    if ip.Results.ExcludeVisitors && size(lftData(i).A,1) > 0
        vidx = getVisitorIndex(lftData(i));
        for f = 1:numel(fnames)
            lftData(i).visitors.(fnames{f}) = lftData(i).(fnames{f})(vidx{1},:,:);
            lftData(i).(fnames{f}) = lftData(i).(fnames{f})(~vidx{1},:,:);
        end
    end
end


av = zeros(nCh,nd);
rmIdx = [];
for c = 1:nCh
    if rescale(c)
        %maxA(c,:) = arrayfun(@(i) nanmax(i.A(:,:,c),[],2), lftData, 'UniformOutput', false);
        [a, offset, refIdx] = scaleEDFs(maxA(c,:), 'Display', ip.Results.DisplayScaling,...
            'FigureName', ['Ch. ' num2str(c) ' max. intensity scaling'],...
            'Colormap', ip.Results.Colormap, 'Legend', getMovieName(data));
        
        av(c,:) = a;
        movieLength = min([data.movieLength]);
        for i = 1:nd
            if size(lftData(i).A,1) > 0
                lftData(i).A = lftData(i).A(:,1:movieLength,:);
                maxA{c,i} = a(i) * maxA{c,i};
                for f = 1:numel(afields)
                    lftData(i).(afields{f})(:,:,c) = a(i)*lftData(i).(afields{f})(:,:,c);
                end
            end
        end
    end
    
    if ip.Results.RemoveOutliers && nd>5
        outlierIdx = detectEDFOutliers(maxA(c,:), offset, refIdx);
        if ~isempty(outlierIdx)
            fprintf('Outlier data sets:\n');
            for i = 1:numel(outlierIdx)
                fprintf('Index %d: %s\n', outlierIdx(i), getShortPath(data(outlierIdx(i))));
            end
            rmv = input('Remove outliers? (y/n) ', 's');
            if strcmpi(rmv, 'y') || isempty(rmv)
                rmIdx = [rmIdx outlierIdx]; %#ok<AGROW>
            end
        end
    end
end
if ~isempty(rmIdx)
    lftData(rmIdx) = [];
    av(:,rmIdx) = [];
end

if rescale(1)
    a = mat2cell(av,nCh,ones(1,numel(lftData)));
    [lftData.a] = deal(a{:});
    for i = 1:numel(lftData)
        lftData(i).maxA = squeeze(nanmax(lftData(i).A(:,:,:),[],2));
    end
end
