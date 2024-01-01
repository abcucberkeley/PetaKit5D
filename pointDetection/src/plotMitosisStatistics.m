%plotMitosisStatistics(data, varargin) plots volume, surface area and initiation
% density statistics for movies from cells undergoing mitosis. Individual sections of
% the time-lapse are autmatically recognized based on the time stamp (create time) of
% the files.

% Francois Aguet, 2014

function plotMitosisStatistics(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Cutoff_f', 4, @isscalar);
ip.addParamValue('MaxIntensityThreshold', 0, @isscalar);
ip.addParameter('ExcludeIndex', []);
ip.addParameter('Frame', []);
ip.addParameter('PrintPath', [], @ischar);
ip.parse(data, varargin{:});

t = data.timestamps;
% start index of recording periods
startIdx = find([true diff(t)>5*data.framerate]);
startTimes = t(startIdx);
% time points in each partial acquisition
npoints = diff([startIdx data.movieLength+1]);

idx = arrayfun(@(k) startIdx(k) + (0:npoints(k)-1), 1:numel(startIdx), 'unif', 0);
t = arrayfun(@(k) startTimes(k) + (0:npoints(k)-1)*data.framerate, 1:numel(startIdx), 'unif', 0);
t = [t{:}];
tmin = t/60;

frameIndex = ip.Results.Frame;
if isempty(frameIndex)
    frameIndex = numel(tmin);
end


% load tracks
tfile = [data.source 'Analysis' filesep 'ProcessedTracks.mat'];
tracks = load(tfile);
tracks = tracks.tracks;
tracks = tracks([tracks.catIdx]==1);

rfile = [data.source 'Analysis' filesep 'voldata.mat'];
M = load(rfile);

lftData = getLifetimeData(data, 'ExcludeVisitors', false, 'ReturnValidOnly', false);

%selIdx = [tracks.lifetime_s]>ip.Results.Cutoff_f*data.framerate & [tracks.catIdx]<4;
selIdx = [tracks.lifetime_s]>ip.Results.Cutoff_f*data.framerate & ...
    max(lftData.A,[],2)'>=ip.Results.MaxIntensityThreshold;
initDens = hist([tracks(selIdx).start], 1:data.movieLength);
% set first frame to NaN
endIdx = [startIdx(2:end)-1 data.movieLength];
initDens([startIdx endIdx]) = NaN;
initDens(initDens==0) = NaN;
dens = initDens./(M.area)*60/data.framerate;

% smoothen density
for k = 1:numel(idx)
    tmp = dens(idx{k});
    dens(idx{k}) = conv(padarrayXT(tmp, [0 1], 'symmetric'), ones(1,3)/3, 'valid');
end


ha = setupFigure(3,1, 'AxesWidth', 8, 'AxesHeight', 2.5,...
    'XSpace', [2 0.5 0.5], 'YSpace', [1.5 1 1], 'SameAxes', true);

labels = arrayfun(@(i) i*ones(1,npoints(i)), 1:numel(npoints), 'unif', 0);
labels = [labels{:}];

c = hsv2rgb([0 1 0.9]);

opts = {'Color', c, 'LineWidth', 1.5};

set(ha, 'XTick', 0:10:100, 'XLim', [0 tmin(end)]);
set(ha(1), 'YLim', [0 5000], 'YTick', 0:1000:10000);
set(ha(2), 'YLim', [0 2000], 'YTick', 0:500:5000);
set(ha(3), 'YLim', [0 0.4], 'YTick', 0:0.1:1);
xlabel('Time (min)');
yh(1) = ylabel(ha(1), ['Volume (' char(181) 'm^3)']);
yh(2) = ylabel(ha(2), ['Surface area (' char(181) 'm^2)']);
yh(3) = ylabel(ha(3), ['# Events (' char(181) 'm^{-2} min^{-1})']);
ypos = get(yh, 'Position');

% ymin = round(min(cellfun(@(k) k(1), ypos)));
ymin = -tmin(end)/6.2817;
arrayfun(@(k) set(yh(k), 'Position', [ymin ypos{k}(2:3)]), 1:3);
linkaxes(ha, 'x');


if ~isempty(ip.Results.PrintPath)
    [~,~] = mkdir(ip.Results.PrintPath);
end

ltext = {'(i)', '(ii)', '(iii)', '(iv)', '(v)', '(vi)'};
for f = frameIndex % loop through frames
    
    % plot complete
    for k = setdiff(1:labels(f)-1, ip.Results.ExcludeIndex)
        plot(ha(1), tmin(idx{k}), M.vol(idx{k}), opts{:});
        plot(ha(2), tmin(idx{k}), M.area(idx{k}), opts{:});
        plot(ha(3), tmin(idx{k}), dens(idx{k}), opts{:});
    end
    
    for k = setdiff(1:numel(idx), ip.Results.ExcludeIndex)
        text(median(tmin(idx{k})), 6000, ltext{k}, 'Color', 'k',...
            'HorizontalAlignment', 'center', 'Parent', ha(1));
    end
    
    % plot current acquisition window
    k = labels(f);
    if ~ismember(k,ip.Results.ExcludeIndex)
        b = f - sum(npoints(1:labels(f)-1));
        plot(ha(1), tmin(idx{k}(1:b)), M.vol(idx{k}(1:b)), opts{:});
        plot(ha(2), tmin(idx{k}(1:b)), M.area(idx{k}(1:b)), opts{:});
        plot(ha(3), tmin(idx{k}(1:b)), dens(idx{k}(1:b)), opts{:});
    end
    
    if ~isempty(ip.Results.PrintPath)
        print('-depsc2', '-loose', [ip.Results.PrintPath filesep 'tmp.eps']);
        cmd = ['export DYLD_LIBRARY_PATH=""; convert -density 144 ' ip.Results.PrintPath filesep ...
            'tmp.eps ' ip.Results.PrintPath filesep 'frame_' num2str(f, fmt) '.png'];
        system(cmd);
    end
    
    if f<frameIndex(end)
        cla(ha(1));
        cla(ha(2));
        cla(ha(3));
    end
end

% % mean +- s.d.
% cellfun(@(k) plot(ha(1), mean(tmin(k))*[1 1], mean(vol(k))+std(vol(k))*[-1 1], 'k',...
%     'LineWidth', 1.5), idx);
% cellfun(@(k) plot(ha(1), mean(tmin(k)), mean(vol(k)), 'ko',...
%     'MarkerSize', 4, 'MarkerFaceColor', 'r', 'LineWidth', 1.5), idx);
% 
% cellfun(@(k) plot(ha(2), mean(tmin(k))*[1 1], mean(area(k))+std(area(k))*[-1 1], 'k',...
%     'LineWidth', 1.5), idx);
% cellfun(@(k) plot(ha(2), mean(tmin(k)), mean(area(k)), 'ko',...
%     'MarkerSize', 4, 'MarkerFaceColor', 'r', 'LineWidth', 1.5), idx);
