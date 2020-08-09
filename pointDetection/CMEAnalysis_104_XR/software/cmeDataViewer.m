%cmeDataViewer(data, varargin) displays movies with associated detection and tracking results.
%
% Inputs:    
%             data : single movie structure returned by loadConditionData.m
%
% Options (specifier/value format):
%     'LoadTracks' : {true}|false specifies whether tracking data is loaded
%     'LoadFrames' : {true}|false specifies whether the raw movie data is loaded
%                    This is required for visualization of tracks overlaid on data
%       'LoadMask' : {true}| false specifies whether the cell outline mask is loaded
%       'Cutoff_f' : Minimum length of tracks to load, in frames. Default: 3
%
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

% Francois Aguet, 2011 (last modified 03/30/2014)

function cmeDataViewer(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParameter('LoadTracks', true, @islogical);
ip.addParameter('LoadFrames', true, @islogical);
ip.addParameter('LoadMask', true, @islogical);
ip.addParameter('RelativePath', 'Tracking', @ischar);
ip.addParameter('Cutoff_f', 3);
ip.parse(data, varargin{:});

% Handles/settings are stored in 'appdata' of the figure handle
handles.data = data;

% detect number of channels (up to 4)
nCh = length(data.channels);
if nCh>4
    error('Max. 4 channels supported.');
end

handles.nCh = nCh; % # channels
handles.mCh = find(strcmp(data.source, data.channels)); % master channel index

fidx = 1; % selected frame
trackIndex = [];%.current = 1; % selected track; absolute index (of all detected tracks)

nx = data.imagesize(2);
ny = data.imagesize(1);
nf = data.movieLength;

xs = round(nx/2); % position of (x,y) cross-section
ys = round(ny/2);
yshift = 0; % vertical shift between R,G channels in RGB display
lcolor = hsv2rgb([0.55 0.5 0.8]); % color of cross-section selector lines

tlFile = [data.source 'Analysis' filesep 'TrackLabels.mat'];

%===============================================================================
% Setup main GUI window/figure
%===============================================================================
bgColor = get(0,'defaultUicontrolBackgroundColor');
hfig = figure('Units', 'normalized', 'Position', [0.025 0.2 0.95 0.8],...
    'PaperPositionMode', 'auto', 'Toolbar', 'figure',...
    'Color', bgColor,...
    'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels', 'Name', [getDirFromPath(getExpDir(data)) filesep getCellDir(data)]);

pos = get(hfig, 'Position'); % [pixels]

%-------------------------------------------------------------------------------
% Control panels at bottom of GUI
%-------------------------------------------------------------------------------
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', '', 'Position', [5 5 650 70], 'BackgroundColor', bgColor);

uicontrol(ph, 'Style', 'text', 'String', 'Data display: ',...
    'Position', [5 40 90 20], 'HorizontalAlignment', 'left');

maskCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Cell mask',...
    'Position', [5 25 100 15], 'HorizontalAlignment', 'left',...
    'Callback', @frameCheck_Callback);
% plot on top
frameChoice = uicontrol(ph, 'Style', 'popup',...
    'String', {'Raw', 'Detections', 'RGB'},...
    'Position', [95 42 100 20], 'Callback', @frameChoice_Callback);

labelCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Channel labels',...
    'Position', [5 5 140 15], 'HorizontalAlignment', 'left',...
    'Callback', @chlabel_Callback);

detectionCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Detections',...
    'Position', [200 50 100 15], 'HorizontalAlignment', 'left',...
    'Callback', @frameCheck_Callback);
trackCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Tracks:', 'Value', true,...
    'Position', [200 30 80 15], 'HorizontalAlignment', 'left',...
    'Callback', @frameCheck_Callback);
trackChoice = uicontrol('Style', 'popup',...
    'String', {'Category', 'Lifetime', 'EAP Status', 'Object Type', 'Random'},...
    'Position', [280 33 110 20], 'Callback', @trackChoice_Callback);
trackRangeButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Settings',...
    'Position', [280 5 80 20], 'HorizontalAlignment', 'left',...
    'Callback', @trackSettings_Callback);


gapCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Gaps',...
    'Position', [390 45 140 15], 'HorizontalAlignment', 'left',...
    'Callback', @frameCheck_Callback);
trackEventCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Births/Deaths',...
    'Position', [390 25 140 15], 'HorizontalAlignment', 'left',...
    'Callback', @frameCheck_Callback);
eapCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'EAP status',...
    'Position', [390 5 140 15], 'HorizontalAlignment', 'left',...
    'Callback', @frameCheck_Callback);


trackSelectButton = uicontrol(ph, 'Style', 'togglebutton', 'String', 'Select track',...
    'Position', [540 40 100 20], 'HorizontalAlignment', 'left',...
    'Callback', @trackSelectButton_Callback);
statsButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Show statistics',...
    'Position', [540 10 100 20], 'HorizontalAlignment', 'left',...
    'Callback', @statsButton_Callback);


%---------------------
% Tracks
%---------------------
% Track plot panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Plot options', 'Position', [pos(3)-220-150-5-210 5 210 70]);
tplotText = uicontrol(ph, 'Style', 'text', 'String', 'Units: ',...
    'Position', [5 35 60 20], 'HorizontalAlignment', 'left');

tplotUnitChoice = uicontrol(ph, 'Style', 'popup',...
    'String', {'Seconds', 'Frames'},...
    'Position', [40 40 100 15], 'Callback', {@unitChoice_Callback, hfig});
tplotBackgroundCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Subtract background',...
    'Position', [5 20 150 15], 'HorizontalAlignment', 'left', 'Value', true, 'Callback', @trackCheck_Callback);
tplotOverlayCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Overlay',...
    'Position', [140 20 60 15], 'HorizontalAlignment', 'left', 'Value', false, 'Callback', @tplotOverlay_Callback);
tplotScaleCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Autoscale',...
    'Position', [5 5 90 15], 'HorizontalAlignment', 'left', 'Value', false, 'Callback', @trackCheck_Callback);
tplotRangeCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Total time',...
    'Position', [95 5 90 15], 'HorizontalAlignment', 'left', 'Value', false, 'Callback', @trackCheck_Callback);
handles.tplotPanel = ph;

% Montage panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Montage plot', 'Position', [pos(3)-220-150 5 220 70]);
montageAlignCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Align to track',...
    'Position', [90 38 115 15], 'HorizontalAlignment', 'left', 'Value', true, 'Callback', @restoreFocus);
montageMarkerCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Show markers',...
    'Position', [90 23 115 15], 'HorizontalAlignment', 'left', 'Callback', @restoreFocus);
montageDetectionCheckbox = uicontrol(ph, 'Style', 'checkbox', 'String', 'Show detection',...
    'Position', [90 8 120 15], 'HorizontalAlignment', 'left', 'Callback', @restoreFocus);
montageButton = uicontrol(ph, 'Style', 'pushbutton','String','Generate',...
    'Units', 'pixels', 'Position', [5 20 80 20],...
    'Callback', @montageButton_Callback);
handles.montagePanel = ph;


% Output panel
ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Output', 'Position', [pos(3)-145 5 140 70]);

handles.printButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Print figures',...
    'Units', 'normalized', 'Position', [0.1 0.5 0.8 0.45],...
    'Callback', {@printButton_Callback, hfig});

handles.movieButton = uicontrol(ph, 'Style', 'pushbutton', 'String', 'Make movie',...
    'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.45],...
    'Callback', @movieButton_Callback);
handles.outputPanel = ph;


setappdata(hfig, 'handles', handles); % write 'handles' to hfig

%===============================================================================
% Set up frame display
%===============================================================================
handles.frameLabel = uicontrol('Style', 'text', 'String', ['Frame ' num2str(fidx)], ...
    'Position', [10 pos(4)-20 100 15], 'HorizontalAlignment', 'left');

setupFrameAxes(nCh);

% Frame slider
w = 320; % fixed width of the track plots, in pixels
lspace = 10;
rspace = w+30+50;
if data.movieLength>1
    handles.frameSlider = uicontrol('Style', 'slider', 'Units', 'pixels',...
        'Value', fidx, 'SliderStep', [1/(nf-1) 0.05], 'Min', 1, 'Max', nf,...
        'Position', [lspace 77 pos(3)-rspace-lspace 18], 'Callback', @frameSliderRelease_Callback);
end
% this definition (instead of regular callback) enable continuous sliding
addlistener(handle(handles.frameSlider), 'Value', 'PostSet', @frameSlider_Callback);


%===============================================================================
% Set up track display
%===============================================================================
handles.trackIndex.label = uicontrol('Style', 'text', 'String', 'Track 1',...
    'Units', 'pixels', 'Position', [pos(3)-70 pos(4)-18 100 15], 'HorizontalAlignment', 'right');

setupTrackAxes(nCh);

% Track slider
h_tot = pos(4) - 140;
handles.trackSlider = uicontrol('Style', 'slider',...
    'Value', 1, 'SliderStep', [0.01 0.05], 'Min', 1, 'Max', 1000,...
    'Position', [pos(3)-24 120 10 h_tot], 'Callback', @trackSliderRelease_Callback);
% this definition (instead of regular callback) enable continuous sliding
addlistener(handle(handles.trackSlider), 'Value', 'PostSet', @trackSlider_Callback);


%===============================================================================
% Set up menu
%===============================================================================
hmenu = uimenu('Label','Options');
uimenu(hmenu, 'Label', 'Annotate', 'Callback', @annotation_Callback);
doAnnotate = false;


% colormap/contrast settings
p = 1;
invertContrast = false;
cmap = gray(256).^p;
colormap(cmap);


setappdata(hfig, 'handles', handles);
set(hfig, 'ResizeFcn', @figResize);

isRGB = false;

hxy = [];
hxz = [];
hyz = [];
hl = [];

% handles for track plot objects in frames window
hpt = []; % tracks
hpd = []; % detections
hpg = []; % gaps
hps = []; % starts/ends

hst = []; % selected track marker

hms = []; % cell mask

% handles for track plots
% ht = [];
hChLabel = [];

trackColormap = [];

displayType = 'raw';
pUnitType = 's';


%-------------------------------------------------------------------------------
% Allocate arrays for RGB display
%-------------------------------------------------------------------------------
rframe = zeros(ny,nx,3,'uint8');
tframe = zeros(nf,nx,3,'uint8');
lframe = zeros(ny,nf,3,'uint8');
idxRGB = getRGBindex(data.markers);

% Change background color (looks bad)
%bgcolor = [1 1 1];
%set(hfig, 'Color', bgcolor);
%set(findobj(hfig,'-property', 'BackgroundColor'), 'BackgroundColor', bgcolor);


%===============================================================================
% Load movie and associated analysis results
%===============================================================================
% readfct = @(path, i) imread(path, i);

stack = cell(1,nCh);
if ip.Results.LoadFrames
    fprintf('Loading frames ... ');
    if ~iscell(data.framePaths{1})
        for c = 1:nCh
            stack{c} = readtiff(data.framePaths{c});
        end
    else
        for c = 1:nCh
            stack{c} = zeros([data.imagesize data.movieLength], 'uint16');
            for i = 1:data.movieLength
                stack{c}(:,:,i) = imread(data.framePaths{c}{i});
            end
        end
    end
    fprintf('done.\n');
end

%-------------------------------------------------------------------------------
% Load detection masks
%-------------------------------------------------------------------------------
if ip.Results.LoadFrames && exist(data.maskPaths, 'file')==2
    fprintf('Loading detection masks ... ');
    dmask = readtiff(data.maskPaths);
    fprintf('done.\n');
elseif exist([data.source 'Detection' filesep 'Masks'], 'dir')==7 && ip.Results.LoadFrames
    fprintf('Loading detection masks ... ');
    dmask = zeros(ny,nx,nf, 'uint8');
    for i = 1:data.movieLength
        dmask(:,:,i) = imread(data.maskPaths{i});
    end
    fprintf('done.\n');
else
    dmask = [];
end

if exist([data.source 'Detection' filesep 'cellmask.tif'], 'file')==2 && ip.Results.LoadMask
    cellMask = imread([data.source 'Detection' filesep 'cellmask.tif']);
else
    cellMask = [];
end

%-------------------------------------------------------------------------------
% Load detection files
%-------------------------------------------------------------------------------
if ip.Results.LoadFrames
    % for c = 1:nCh
    detectionFile = [data.channels{1} 'Detection' filesep 'detection_v2.mat'];
    if (exist(detectionFile, 'file')==2)
        frameInfo = load(detectionFile);
        frameInfo = frameInfo.frameInfo;
    else
        frameInfo = [];
    end
    % end
end

%-------------------------------------------------------------------------------
% Load tracks
%-------------------------------------------------------------------------------
tracks = [];
%bgA = [];
maxA = [];
if ip.Results.LoadTracks
    fprintf('Loading tracks ... ');
    
    % identify track file
    fileList = dir([data.source ip.Results.RelativePath filesep 'ProcessedTracks*.mat']);
    fileList = {fileList.name};
    if numel(fileList)>1
        idx = 0;
        while ~(idx>=1 && idx<=numel(fileList) && round(idx)==idx)
            fprintf('Tracking results found for this data set:\n');
            for i = 1:numel(fileList)
                fprintf('[%d] %s\n', i, fileList{i});
            end
            idx = str2double(input('Please enter the number of the set to load: ', 's'));
        end
        fileName = fileList{idx};
    elseif numel(fileList)==1
        fileName = fileList{1};
    else
        fileName = [];
    end
end

if ip.Results.LoadTracks && exist([data.source ip.Results.RelativePath filesep fileName], 'file')==2
    tmp = load([data.source ip.Results.RelativePath filesep fileName]);
    tracks = tmp.tracks;
    %if isfield(tmp, 'bgA')
    %    bgA = cellfun(@(i) prctile(i, 95, 2), tmp.bgA, 'unif', 0);
    %    bgA = [bgA{:}];
    %end
    clear tmp;
    
    [~, sortIdx] = sort([tracks.lifetime_s], 'descend');
    tracks = tracks(sortIdx);
    nt = numel(tracks);
    trackIndex.all = 1:nt;
    trackIndex.lft = [tracks.lifetime_s];
    trackIndex.cat = [tracks.catIdx];
    
    % apply frame cutoff, cell mask
    idx = [tracks.lifetime_s] >= data.framerate*ip.Results.Cutoff_f;
    if ~isempty(cellMask)
        x = NaN(1,nt);
        y = NaN(1,nt);
        for t = 1:nt
            x(t) = round(nanmean(tracks(t).x(1,:)));
            y(t) = round(nanmean(tracks(t).y(1,:)));
        end
        idx = idx & (cellMask(sub2ind([ny nx], y, x))==1);
    end
    
    % track labels (GUI annotation)
    if exist(tlFile,'file')==2
        tmp = load(tlFile);
        trackIndex.label = tmp.trackIndex.label;
    else
        % label vector (for all tracks)
        trackIndex.label = zeros(1,nt);
        trackIndex.label(idx & [tracks.catIdx]==1) = 1; % all loaded & valid
    end
    
    % index of loaded tracks
    trackIndex.loaded = find(idx);
    tracks = tracks(idx);
    nt = numel(tracks);
    
    
    minLft = min([tracks.lifetime_s]);
    maxLft = max([tracks.lifetime_s]);
    % slider values; minLft, maxLft are const
    minVal = max(data.framerate*5, minLft); % Display tracks >= 5 frames by default
    maxVal = maxLft;
    
    % initially tracks with lifetime >= minVal are selected
    trackIndex.selected = [tracks.lifetime_s]>=minVal;
    % current track: pointer to 'loaded' array
    trackIndex.current = find(trackIndex.selected, 1, 'first');
    
    nseg = [tracks.nSeg];
    
    np = sum(nseg);
    X = NaN(nf, np);
    Y = NaN(nf, np);
    G = false(nf, np);
    % for significance values, store vectors
    mvec = [tracks.hval_Ar];
    if isfield(tracks, 'significantVsBackground')
        svec = [tracks.significantVsBackground];
    else
        svec = [];
    end
    fvec = [tracks.f];
    xvec = [tracks.x];
    yvec = [tracks.y];
   
    % vector of start indexes since multiple segments/track
    tidx = cumsum([1 nseg(1:end-1)]);
    
    trackStarts = [tracks.start];
    trackEnds = [tracks.end];
    mu_x = NaN(1,nt);
    mu_y = NaN(1,nt);
    
    for t = 1:nt
        if nseg(t)==1
            X(tracks(t).f, tidx(t)) = tracks(t).x(1,:);
            Y(tracks(t).f, tidx(t)) = tracks(t).y(1,:);
            G(tracks(t).f, tidx(t)) = tracks(t).gapVect;
            mu_x(t) = nanmean(X(:,tidx(t)));
            mu_y(t) = nanmean(Y(:,tidx(t)));           
        else
            sep = find(isnan(tracks(t).t));
            sep = [0 sep numel(tracks(t).f)+1]; %#ok<AGROW>
            for s = 1:tracks(t).nSeg
                sidx = sep(s)+1:sep(s+1)-1;
                X(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).x(1,sidx);
                Y(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).y(1,sidx);
                G(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).gapVect(sidx);
            end
            mu_x(t) = nanmean(nanmean(X(:,tidx(t):tidx(t)+s-1)));
            mu_y(t) = nanmean(nanmean(Y(:,tidx(t):tidx(t)+s-1)));
        end
    end
    
    % segment index
    % Example: [1 1 2 3 4 4 ... ] first two cols are from same track
    idx = diff([tidx size(X,2)+1]);
    idx = arrayfun(@(i) i+zeros(1, idx(i)), 1:numel(idx), 'unif', 0);
    trackIndex.segmentIndex = [idx{:}];
    
    % min/max track intensities
    maxA = arrayfun(@(t) max(t.A, [], 2), tracks, 'unif', 0);
    maxA = [maxA{:}];
    maxInt = prctile(maxA, 99, 2);
    da = floor(log10(maxInt));
    % y-axis unit
    yunit = round(maxInt ./ 10.^da) .* 10.^(da-1);
    maxInt = ceil(maxInt./yunit) .* yunit;
end
fprintf('done.\n');

% Settings
if ~isempty(tracks)
    catCheckVal = ones(1,8); % category selection
    eapCheckVal = ones(1,3); % EAP selection
    mitCheckVal = 0; % max. intensity threshold
    setRadioVal = 1;
    paramFile = [data.source 'Analysis' filesep 'parameters.mat'];
    if exist(paramFile, 'file')==2
        prm = load(paramFile);
        maxIntT = prm.MaxIntensityThreshold / prm.a(1);
    else
        maxIntT = 0;
    end
end


% dynamic range for each channel
dRange = cell(1,nCh);
for c = 1:nCh
    %dRange{c} = double([min(stack{c}(:)) max(stack{c}(:))]);
    tmp = cat(3, stack{c}(:,:,[1 end]));
    dRange{c} = prctile(double(tmp(:)), [0.001 99.99]);
end

hues = getFluorophoreHues(data.markers);
rgbColors = arrayfun(@(x) hsv2rgb([x 1 1]), hues, 'unif', 0);


%===============================================================================
% Set visibility for sliders and checkboxes
%===============================================================================
if ~isempty(tracks)
    ntSel = sum([trackIndex.selected]);
    if ntSel>1
        set(handles.trackSlider, 'Min', 1);
        set(handles.trackSlider, 'Max', ntSel);
        set(handles.trackSlider, 'SliderStep', [1/(ntSel-1) 0.05]);
    end
    setTrackColormap('Category');
    setColorbar('Category');
else
    set(handles.trackSlider, 'Visible', 'off');
    set(handles.trackIndex.label, 'Visible', 'off');
    set(handles.tAxes, 'Visible', 'off');
    set(trackSelectButton, 'Enable', 'off');
    set(statsButton, 'Enable', 'off');
    set(trackCheckbox, 'Value', false);
    set([trackCheckbox trackChoice trackRangeButton gapCheckbox trackEventCheckbox eapCheckbox], 'Enable', 'off');
    set([montageAlignCheckbox montageMarkerCheckbox montageDetectionCheckbox montageButton], 'Enable', 'off');
    %set(handles.montagePanel, 'Visible', 'off');
    
    set(hLegend, 'Visible', 'off');
    set([tplotText tplotUnitChoice tplotBackgroundCheckbox tplotOverlayCheckbox tplotScaleCheckbox tplotRangeCheckbox], 'Enable', 'off');
end

if ~ip.Results.LoadFrames
    set(handles.fAxes, 'Visible', 'off');
    set(hLegend, 'Visible', 'off');
    set(handles.frameSlider, 'Visible', 'off');
end

if isempty(cellMask)
    set(maskCheckbox, 'Enable', 'off');
end

if nCh==1
    set(eapCheckbox, 'Enable', 'off');
end

%===============================================================================
% populate with data, plotting functions are called only here, afterwards change data
%===============================================================================
if ip.Results.LoadFrames
    initStackviewer();
end

% current frame and track markers
hf = zeros(nCh,1);
if ~isempty(tracks)
    for c = 1:nCh
        % plot current frame marker
        hf(c) = plot(handles.tAxes(c), ([fidx fidx]-1)*data.framerate,...
            get(handles.tAxes(c), 'YLim'), '--', 'Color', 0.7*[1 1 1]);
    end
    updateTrack();
end


%===============================================================================
% Set listeners
%===============================================================================
set(hfig, 'WindowScrollWheelFcn', @scroll_Callback);
set(hfig, 'KeyPressFcn', @key_Callback);
set(hfig, 'CloseRequestFcn', @close_Callback);

hpan = pan;
set(hpan,'ActionPreCallback',@panstart);
set(hpan,'ActionPostCallback',@panstop);

hz = zoom;
set(hz, 'ActionPostCallback', @czoom);
% setAxesZoomMotion(hz, handles.tAxes, 'horizontal');
% linkaxes(handles.tAxes, 'x');


%===============================================================================
% Listener/display functions
%===============================================================================
    function axesClick_Callback(varargin)
        if get(trackSelectButton, 'Value')==0
            switch gca
                case num2cell(handles.fAxes(:,1))
                    a = get(gca, 'CurrentPoint');
                    xs = round(a(1,1));
                    ys = round(a(1,2));
                    updateProj(); % required for clicking w/o dragging
                    set(gcf, 'WindowButtonMotionFcn', @dragProj, 'WindowButtonUpFcn', @stopDragging);
                case num2cell(handles.fAxes(:,2))
                    a = get(gca,'CurrentPoint');
                    fidx = round(a(1,1));
                    updateSlice();
                    set(gcf, 'WindowButtonMotionFcn', @dragSlice, 'WindowButtonUpFcn', @stopDragging);
                case num2cell(handles.fAxes(:,3))
                    a = get(gca,'CurrentPoint');
                    fidx = round(a(1,2));
                    updateSlice();
                    set(gcf, 'WindowButtonMotionFcn', @dragSlice, 'WindowButtonUpFcn', @stopDragging);
            end
        elseif gca==handles.fAxes(1,1) && ~isempty(tracks) && get(trackCheckbox, 'Value')
            % track selection mode
            a = get(gca, 'CurrentPoint');
            x0 = a(1,1);
            y0 = a(1,2);
            
            % track segments visible in current frame
            cidx = find([tracks.start]<=fidx & fidx<=[tracks.end] & trackIndex.selected);
            if ~isempty(cidx)
                % distance to mean of tracks
                d = sqrt((x0-mu_x(cidx)).^2 + (y0-mu_y(cidx)).^2);
                [~,d] = nanmin(d);
                
                % update selected (current) track (pointer to 'loaded')
                trackIndex.current = cidx(d);
                set(handles.trackSlider, 'Value', find(trackIndex.current==find(trackIndex.selected)));% calls updateTrack
                set(hst, 'XData', X(fidx, trackIndex.segmentIndex==trackIndex.current), 'YData', Y(fidx, trackIndex.segmentIndex==trackIndex.current));
            end
        end
    end


    function dragProj(varargin)
        a = get(gca, 'CurrentPoint');
        xs = round(a(1,1));
        ys = round(a(1,2));
        updateProj();
    end


    function dragSlice(varargin)
        a = get(gca, 'CurrentPoint');
        switch gca
            case num2cell(handles.fAxes(:,2))
                fidx = min(max(1,round(a(1,1))),nf);
            case num2cell(handles.fAxes(:,3))
                fidx = min(max(1,round(a(1,2))),nf);
        end
        set(handles.frameSlider, 'Value', fidx);
        updateSlice();
    end


    function stopDragging(varargin)
        set(gcf, 'WindowButtonMotionFcn', '');
    end


    % scroll through stack slices
    function scroll_Callback(~, eventdata)
        if eventdata.VerticalScrollCount < 0
            if fidx < nf
                fidx = fidx + 1;
                set(handles.frameSlider, 'Value', fidx);
                updateSlice();
            end
        elseif eventdata.VerticalScrollCount > 0
            if fidx > 1
                fidx = fidx - 1;
                set(handles.frameSlider, 'Value', fidx);
                updateSlice();
            end
        end
    end


    function key_Callback(~, eventdata)
        switch eventdata.Key
            case 'leftarrow'
                if fidx > 1
                    fidx = fidx - 1;
                    set(handles.frameSlider, 'Value', fidx);
                    updateSlice();
                end
            case 'rightarrow'
                if fidx < nf
                    fidx = fidx + 1;
                    set(handles.frameSlider, 'Value', fidx);
                    updateSlice();
                end
            case 'downarrow'
                cidx = get(handles.trackSlider, 'Value');
                if cidx>1
                    % invokes trackSlider_callback:
                    set(handles.trackSlider, 'Value', cidx-1);
                end
            case 'uparrow'
                cidx = get(handles.trackSlider, 'Value');
                if cidx < get(handles.trackSlider, 'Max')
                    % invokes trackSlider_callback:
                    set(handles.trackSlider, 'Value', cidx+1);
                end
            case 'comma'
                if strcmpi(displayType, 'RGB')
                    yshift = yshift-1;
                    updateSlice();
                end
            case 'period'
                if strcmpi(displayType, 'RGB')
                    yshift = yshift+1;
                    updateSlice();
                end
            case 'slash'
                if strcmpi(displayType, 'RGB')
                    yshift = 0;
                    updateSlice();
                end
            case 'c'
                if p==1
                    p = 0.5;
                else
                    p = 1;
                end
                setColormap();
                
            case 'i'
                invertContrast = ~invertContrast;
                setColormap();
        end
        if all(isstrprop(eventdata.Key, 'digit')) && doAnnotate
            nkey = str2double(eventdata.Key);
            if nkey>=0 && nkey<=9
                trackIndex.label(trackIndex.loaded(trackIndex.current)) = nkey;
                set(handles.trackIndex.label, 'String', ['Track: ' num2str(trackIndex.loaded(trackIndex.current))...
                    ' (Label: ' num2str(trackIndex.label(trackIndex.loaded(trackIndex.current))) ')']);
            end
        end
    end


    function setColormap()
        if ~strcmpi(displayType, 'RGB')
            if ~invertContrast
                cmap = gray(256).^p;
            else
                cmap = gray(256).^(1/p);
                cmap = cmap(end:-1:1,:,:);
            end
            colormap(cmap);
        else % if RGB
            updateSlice();
            updateProj();
        end
    end


    function annotation_Callback(varargin)
        if strcmp(get(gcbo, 'Checked'),'on')
            set(gcbo, 'Checked', 'off');
        else
            set(gcbo, 'Checked', 'on');
        end
        
        switch get(gcbo, 'Label')
            case 'Annotate'
                doAnnotate = strcmp(get(gcbo, 'Checked'),'on');
                % update track annotation display
                if doAnnotate
                    set(handles.trackIndex.label, 'String', ['Track: ' num2str(trackIndex.loaded(trackIndex.current))...
                        ' (Label: ' num2str(trackIndex.label(trackIndex.loaded(trackIndex.current))) ')']);
                else
                    set(handles.trackIndex.label, 'String', ['Track: ' num2str(trackIndex.loaded(trackIndex.current))]);
                end
        end
    end


    function close_Callback(varargin)
        % save track labels, if annotation active
        if doAnnotate
            save(tlFile, 'trackIndex');
        end
        delete(hfig);
    end
        

    function initStackviewer()
        N = numel(handles.fPanels);
        hxy = zeros(1,N);
        hyz = zeros(1,N);
        hxz = zeros(1,N);
        hl = zeros(N,4); % cross-section selectors

        switch displayType
            case 'raw'
                for ci = 1:N
                    hxy(ci) = imagesc(stack{ci}(:,:,fidx), 'Parent', handles.fAxes(ci,1), 'HitTest', 'off');
                    hyz(ci) = imagesc(squeeze(stack{ci}(:,xs,:)), 'Parent', handles.fAxes(ci,2), 'HitTest', 'off');
                    hxz(ci) = imagesc(squeeze(stack{ci}(ys,:,:))', 'Parent', handles.fAxes(ci,3), 'HitTest', 'off');
                end
            case 'mask'
                hxy(1) = imagesc(rgbOverlay(stack{1}(:,:,fidx), dmask(:,:,fidx), [1 0 0], dRange{1}),...
                    'Parent', handles.fAxes(1,1), 'HitTest', 'off');
                for ci = 2:nCh
                    hxy(ci) = imagesc(stack{ci}(:,:,fidx), 'Parent', handles.fAxes(ci,1), 'HitTest', 'off');
                end
                % detection mask not shown in side projections
                for ci = 1:N
                    hyz(ci) = imagesc(squeeze(stack{ci}(:,xs,:)), 'Parent', handles.fAxes(ci,2), 'HitTest', 'off');
                    hxz(ci) = imagesc(squeeze(stack{ci}(ys,:,:))', 'Parent', handles.fAxes(ci,3), 'HitTest', 'off');
                end
            case 'RGB'
                % single panel
                rframe(:,:,idxRGB(2)) = uint8(scaleContrast(...
                    [zeros(max(yshift,0),nx); double(stack{2}(1+max(-yshift,0):end-max(yshift,0),:,fidx)).^p; zeros(max(-yshift,0),nx)], dRange{2}.^p));
                for ci = setdiff(1:nCh,2)
                    rframe(:,:,idxRGB(ci)) = uint8(scaleContrast(double(stack{ci}(:,:,fidx)).^p, dRange{ci}.^p));
                end
                hxy(1) = imagesc(rframe, 'Parent', handles.fAxes(1,1), 'HitTest', 'off');
                for ci = 1:nCh
                    tframe(:,:,idxRGB(ci)) = uint8(scaleContrast(double(squeeze(stack{ci}(ys,:,:))').^p, dRange{ci}.^p));
                    lframe(:,:,idxRGB(ci)) = uint8(scaleContrast(double(squeeze(stack{ci}(:,xs,:))).^p, dRange{ci}.^p));
                end
                hyz(1) = imagesc(lframe, 'Parent', handles.fAxes(1,2), 'HitTest', 'off');
                hxz(1) = imagesc(tframe, 'Parent', handles.fAxes(1,3), 'HitTest', 'off');
        end
        
        for ci = 1:N
            % x,y view
            hold(handles.fAxes(ci,1), 'on');
            set(handles.fAxes(ci,1), 'ButtonDownFcn', @axesClick_Callback);
            hl(ci,1) = plot(handles.fAxes(ci,1), [xs xs], [0.5 ny+0.5], 'Color', lcolor, 'HitTest', 'off', 'DisplayName', 'FrameMarker');
            hl(ci,2) = plot(handles.fAxes(ci,1), [0.5 nx+0.5], [ys ys], 'Color', lcolor, 'HitTest', 'off', 'DisplayName', 'FrameMarker');
            
            % y,z view
            hold(handles.fAxes(ci,2), 'on');
            % line in y,z view
            hl(ci,3) = plot(handles.fAxes(ci,2), fidx*[1 1], [0.5 ny+0.5], 'Color', lcolor, 'HitTest', 'off');
            hold(handles.fAxes(ci,2), 'off');
            
            % x,z view
            hold(handles.fAxes(ci,3), 'on');
            % line in x,z view
            hl(ci,4) = plot(handles.fAxes(ci,3), [0.5 nx+0.5], fidx*[1 1], 'Color', lcolor, 'HitTest', 'off');
            hold(handles.fAxes(ci,3), 'off');
            
            arrayfun(@(i) caxis(i, dRange{ci}), handles.fAxes(ci,:), 'unif', 0);
        end
        
        set(handles.fAxes, 'XTick', [], 'YTick', []);
        axis(handles.fAxes(:,1), 'equal');
        % this fixes a bug with axis 'equal' that allows panning beyond boundaries
        set(handles.fAxes(:,1), 'XLim', [0.5 nx+0.5], 'YLim', [0.5 ny+0.5]);
        
        set(handles.fAxes, 'ButtonDownFcn', @axesClick_Callback);
        
        % plot current track marker
        if ~isempty(tracks)
            hst = zeros(1,N);
            for ci = 1:N
                hst(ci) = plot(handles.fAxes(ci,1), X(fidx, trackIndex.segmentIndex==trackIndex.current),...
                    Y(fidx, trackIndex.segmentIndex==trackIndex.current), 'ws', 'DisplayName', 'TrackMarker', 'MarkerSize', 12, 'HitTest', 'off');
            end
        end
        
        
        if ~isRGB
            set(labelCheckbox, 'Enable', 'on');
            dx = 0.03;
            hChLabel = zeros(1,nCh);
            for ci = 1:nCh
                hChLabel(ci) = text(1-dx*ny/nx, dx, data.markers{ci},...
                    'Color', rgbColors{ci}, 'Units', 'normalized',...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
                    'Parent', handles.fAxes(ci,1), 'HitTest', 'off');
            end
            if nCh<=2 && ~get(labelCheckbox, 'Value')
                set(hChLabel, 'Visible', 'off');
            end
        else
            set(labelCheckbox, 'Enable', 'off');
        end
    end


    function updateSlice(varargin)       
        switch displayType
            case 'raw'                
                for ci = 1:nCh
                    set(hxy(ci), 'CData', stack{ci}(:,:,fidx));
                end
            case 'mask'
                set(hxy(1), 'CData', rgbOverlay(stack{1}(:,:,fidx), dmask(:,:,fidx), [1 0 0], dRange{1}));
                for ci = 2:nCh
                    set(hxy(ci), 'CData', stack{ci}(:,:,fidx));
                end
            case 'RGB'
                rframe(:,:,idxRGB(2)) = uint8(scaleContrast(...
                    [zeros(max(yshift,0),nx); double(stack{2}(1+max(-yshift,0):end-max(yshift,0),:,fidx)).^p; zeros(max(-yshift,0),nx)], dRange{2}.^p));
                for ci = setdiff(1:nCh,2)
                    rframe(:,:,idxRGB(ci)) = uint8(scaleContrast(double(stack{ci}(:,:,fidx)).^p, dRange{ci}.^p));
                end
                set(hxy(1), 'CData', rframe);
        end
        
        set(hl(:,3), 'XData', fidx*[1 1]);
        set(hl(:,4), 'YData', fidx*[1 1]);        
        set(handles.frameLabel, 'String', ['Frame ' num2str(fidx)]);
        
        if all(ishandle(hms))
            delete(hms);
        end
        hms = [];
        if ~isempty(cellMask) && get(maskCheckbox, 'Value')
            B = bwboundaries(cellMask);
            for ci = 1:numel(handles.fPanels)
                hms = [hms; cellfun(@(i) plot(handles.fAxes(ci,1), i(:,2), i(:,1), 'Color', 'r', 'LineWidth', 1), B, 'UniformOutput', false)]; %#ok<AGROW>
            end
        end
        
        if ~isempty(tracks)
            % update current frame marker in track plots
            for ci = 1:max(1, nCh*~get(tplotOverlayCheckbox, 'Value'))
                set(hf(ci), 'XData', ([fidx fidx]-1)*data.framerate,...
                    'YData', get(handles.tAxes(ci), 'YLim'));
            end
        end
        if all(ishandle(hpt))
            delete(hpt);
        end
        if all(ishandle(hpg))
            delete(hpg);
        end
        if all(ishandle(hps))
            delete(hps);
        end
        hpt = [];
        hpg = [];
        hps = [];

        if ~isempty(tracks) && fidx~=1 && get(trackCheckbox, 'Value') && any(~isnan(X(fidx,:)) & trackIndex.selected(trackIndex.segmentIndex))
            vidx = ~isnan(X(fidx,:)) & trackIndex.selected(trackIndex.segmentIndex);
            delete(hpt);
%             set(handles.fAxes(1,1), 'ColorOrder', trackColormap(trackIndex.segmentIndex(vidx),:));
            set(handles.fAxes(1,1), 'ColorOrder', trackColormap(trackIndex.segmentIndex(vidx),:));
            hpt = plot(handles.fAxes(1,1), X(1:fidx,vidx), Y(1:fidx,vidx), 'HitTest', 'off');
            if get(gapCheckbox, 'Value')
                hpg = plot(handles.fAxes(1,1), X(fidx,vidx & G(fidx,:)), Y(fidx,vidx & G(fidx,:)), 'o', 'Color', 'w', 'MarkerSize', 6, 'LineWidth', 1);
            end
            if get(trackEventCheckbox, 'Value')
                % Births
                bcoord = arrayfun(@(i) [i.x(1,1) i.y(1,1)], tracks(trackStarts==fidx & trackIndex.selected), 'unif', 0);
                bcoord = vertcat(bcoord{:});
                if~isempty(bcoord)
                    hps = plot(handles.fAxes(1,1), bcoord(:,1), bcoord(:,2), '*', 'Color', 'g', 'MarkerSize', 8, 'LineWidth', 1);
                end
                
                % Deaths
                dcoord = arrayfun(@(i) [i.x(1,1) i.y(1,1)], tracks(trackEnds==fidx & trackIndex.selected), 'unif', 0);
                dcoord = vertcat(dcoord{:});
                if ~isempty(dcoord)
                    hps = [hps; plot(handles.fAxes(1,1), dcoord(:,1), dcoord(:,2), 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 1)];
                end
            end
            set(hst, 'Visible', 'on', 'XData', X(fidx, trackIndex.segmentIndex==trackIndex.current), 'YData', Y(fidx, trackIndex.segmentIndex==trackIndex.current));
        else
            set(hst, 'Visible', 'off');
        end
        
        if ~isempty(tracks) && get(eapCheckbox, 'Value') && numel(handles.fPanels)>1
            for ci = 2:nCh
                sel = fvec==fidx & mvec(ci,:)==1;
                hp1 = plot(handles.fAxes(ci,1), xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', hsv2rgb([1/3 1 0.9]), 'MarkerSize', 8);
                sel = fvec==fidx & mvec(ci,:)==0 & svec(ci,:)==1;
                hp2 = plot(handles.fAxes(ci,1), xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', hsv2rgb([0.55 1 0.9]), 'MarkerSize', 8);
                sel = fvec==fidx & mvec(ci,:)==0 & svec(ci,:)==0;
                hp3 = plot(handles.fAxes(ci,1), xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', 0.8*[1 1 1], 'MarkerSize', 8);
                hps = [hps; hp1; hp2; hp3]; %#ok<AGROW>
            end
        end
        
        if all(ishandle(hpd))
            delete(hpd); % clear previous plots
        end
        hpd = [];
        if get(detectionCheckbox, 'Value') && ~isempty(frameInfo) && ~isempty(frameInfo(fidx).x)
            isPSF = frameInfo(fidx).isPSF(1,:)==1;
            if any(isPSF)
                hpd(1) = plot(handles.fAxes(1,1), frameInfo(fidx).x(1,isPSF), frameInfo(fidx).y(1,isPSF), 'o', 'Color', [0 0.6 0], 'MarkerSize', 8);
            end
            if any(~isPSF)
                hpd(2) = plot(handles.fAxes(1,1), frameInfo(fidx).x(1,~isPSF), frameInfo(fidx).y(1,~isPSF), 'o', 'Color', [0.6 0 0], 'MarkerSize', 8);
            end
        end
    end


    function updateProj() % (x,y) projections of stack
        % plot lines
        set(hl(:,1), 'XData', xs*[1 1]);
        set(hl(:,2), 'YData', ys*[1 1]);
        
        % update data
        xi = min(max(xs,1), nx);
        yi = min(max(ys,1), ny);
        
        switch displayType
            case 'RGB'
                for ci = 1:nCh
                    tframe(:,:,idxRGB(ci)) = uint8(scaleContrast(double(squeeze(stack{ci}(yi,:,:))').^p, dRange{ci}.^p));
                    lframe(:,:,idxRGB(ci)) = uint8(scaleContrast(double(squeeze(stack{ci}(:,xi,:))).^p, dRange{ci}.^p));
                end
                set(hxz(1), 'CData', tframe);
                set(hyz(1), 'CData', lframe);
            otherwise
                for ci = 1:nCh
                    set(hyz(ci), 'CData', squeeze(stack{ci}(:,xi,:)));
                    set(hxz(ci), 'CData', squeeze(stack{ci}(yi,:,:))');
                end
        end
    end


    function frameChoice_Callback(varargin)
        contents = cellstr(get(frameChoice,'String'));
        newAxes = false;
        isRGB = false;
        switch contents{get(frameChoice,'Value')}
            case 'Raw'
                displayType = 'raw';
                if numel(handles.fPanels)~=nCh
                    setupFrameAxes(nCh);
                    newAxes = true;
                end
            case 'RGB'
                displayType = 'RGB';
                if numel(handles.fPanels)>1
                    setupFrameAxes(1);
                    newAxes = true;
                end
                isRGB = true;
            case 'Detections'
                displayType = 'mask';
                if numel(handles.fPanels)~=nCh
                    setupFrameAxes(nCh);
                    newAxes = true;
                end
        end
        if newAxes
            initStackviewer(); % re-populate axes with data
            if isempty(tracks) || ~ip.Results.LoadFrames
                set(hLegend, 'Visible', 'off');
            end
            updateSlice();
        else
            updateSlice();
            updateProj();
        end
        str = cellstr(get(trackChoice, 'String'));
        str = str{get(trackChoice,'Value')};
        setColorbar(str);        
        restoreFocus();       
    end


    function restoreFocus(varargin)
        if nargin < 1
            hi = gcbo;
        else
            hi = varargin{1};
        end
        set(hi, 'enable', 'off');
        drawnow;
        set(hi, 'enable', 'on');
    end


    function frameCheck_Callback(varargin)
        restoreFocus();
        updateSlice();
    end


    function trackCheck_Callback(varargin)
        restoreFocus();
        updateTrack();
    end


    function tplotOverlay_Callback(varargin)
        if get(tplotOverlayCheckbox, 'Value')
            setupTrackAxes(1);
        else
            setupTrackAxes(nCh);
        end
        updateTrack();
        restoreFocus();
    end


    function czoom(~, eventdata)
        % identify panel
        ci = handles.fAxes(:,1) == eventdata.Axes;
        if any(ci) %&& nCh>1 % x,y axes zoomed
            XLim = get(handles.fAxes(ci,1), 'XLim');
            YLim = get(handles.fAxes(ci,1), 'YLim');
            set(handles.fAxes(:,1), 'XLim', XLim, 'YLim', YLim);
            set(handles.fAxes(:,2), 'YLim', YLim);
            set(handles.fAxes(:,3), 'XLim', XLim);
        end
        
    end

    % Pan functions
    function panstart(~, eventdata)
        set(hfig, 'WindowButtonMotionFcn', {@dopan, eventdata});
    end

    function panstop(varargin)
        set(hfig, 'WindowButtonMotionFcn', '');
    end

    function dopan(~,~,eventdata)
        % get limits of current axes
        XLim = get(eventdata.Axes, 'XLim');
        YLim = get(eventdata.Axes, 'YLim');
        
        switch find(any(handles.fAxes == eventdata.Axes,1))
            case 1
                set(handles.fAxes(:,1), 'XLim', XLim, 'YLim', YLim);
                set(handles.fAxes(:,2), 'YLim', YLim);
                set(handles.fAxes(:,3), 'XLim', XLim);
            case 2
                set(handles.fAxes(:,[1 2]), 'YLim', YLim);
            case 3
                set(handles.fAxes(:,[1 3]), 'XLim', XLim);
        end
    end


    function frameSlider_Callback(~, eventdata)
        fidx = round(eventdata.AffectedObject.Value);
        updateSlice();
    end


    function frameSliderRelease_Callback(varargin)
        restoreFocus(handles.frameSlider);
    end


    function trackSlider_Callback(~, eventdata)
%         obj = eventdata.AffectedObject;
        t0 = round(eventdata.AffectedObject.Value); % slider index
        tmp = find(trackIndex.selected);
        trackIndex.current = tmp(t0);
        updateTrack();
        
        % if track not visible in frame plot, jump to first frame
        % t = handles.tracks{1}(t);
        % if fidx < t.start || fidx > t.end
        %     fidx = t.start;
        %     % set frame number
        %     set(handles.frameLabel, 'String', ['Frame ' num2str(fidx)]);
        %     % set frame slider
        %     set(handles.frameSlider, 'Value', fidx);
        % end
    end


    function trackSliderRelease_Callback(varargin)
        restoreFocus(handles.trackSlider);
    end


    function updateTrack(varargin)
        if ~isempty(tracks) && ~isempty(trackIndex.current)
            % update label
            if ~doAnnotate
                set(handles.trackIndex.label, 'String', ['Track: ' num2str(trackIndex.loaded(trackIndex.current))]);
            else
                set(handles.trackIndex.label, 'String', ['Track: ' num2str(trackIndex.loaded(trackIndex.current))...
                    ' (Label: ' num2str(trackIndex.label(trackIndex.loaded(trackIndex.current))) ')']);
            end
            
            % update selected track marker position
            set(hst, 'XData', X(fidx, trackIndex.segmentIndex==trackIndex.current), 'YData', Y(fidx, trackIndex.segmentIndex==trackIndex.current));
            
            
            itrack = tracks(trackIndex.current);
            if get(tplotBackgroundCheckbox, 'Value')
                bgMode = 'zero';
            else
                bgMode = 'data';
            end
            if strcmpi(pUnitType, 'f')
                itrack.t = itrack.f;
                if ~isempty(itrack.startBuffer)
                    itrack.startBuffer.t = itrack.f(1) - (numel(itrack.startBuffer.t):-1:1);
                    itrack.endBuffer.t = itrack.f(end) + (1:numel(itrack.startBuffer.t));
                end
            end
            
            % plot options
            topts0 = {'Time', 'Movie', 'BackgroundValue', bgMode};
            if get(tplotRangeCheckbox, 'Value')
                topts0 = [topts0, 'XLim', [-2 data.movieLength+1]*data.framerate];
            end
            
            % axes handles
            if ~get(tplotOverlayCheckbox, 'Value')
                ha = handles.tAxes(1:nCh);
            else
                ha = handles.tAxes(1)+zeros(1,nCh);
                % overlay-specific options
                topts0 = [topts0, 'Background', 'off'];
                if get(tplotScaleCheckbox, 'Value')
                    [~,k] = max(maxInt);
                    topts0 = [topts0, 'YTick', -yunit(k):yunit(k):maxInt(k)];
                end
            end
            
            arrayfun(@(i) hold(i, 'off'), ha); % reset contents
            for ci = 1:nCh
                topts = [topts0, 'Handle', ha(ci)];
                if get(tplotScaleCheckbox, 'Value') && ~get(tplotOverlayCheckbox, 'Value')
                    topts = [topts, 'YTick', -yunit(ci):yunit(ci):maxInt(ci)]; %#ok<AGROW>
                end
                
                % plot background confidence level (based on secondary channel)
                %if ~isempty(bgA) && itrack.catIdx<5
                %   conf = bgA(ci, itrack.f);
                %    if ~isempty(itrack.startBuffer)
                %        conf = [bgA(ci, itrack.startBuffer.f) conf]; %#ok<AGROW>
                %    end
                %    if ~isempty(itrack.endBuffer)
                %        conf = [conf bgA(ci, itrack.endBuffer.f)]; %#ok<AGROW>
                %    end
                %    topts = [topts 'BackgroundConfidence', conf]; %#ok<AGROW>
                %end
                
                plotTrack(data, itrack, ci, topts{:});
                
                % plot current frame indicator
                if ci==1 || ~get(tplotOverlayCheckbox, 'Value')
                    hf(ci) = plot(ha(ci), ([fidx fidx]-1)*data.framerate,...
                        get(ha(ci), 'YLim'), '--', 'Color', 0.7*[1 1 1]);
                end
                set(ha(ci), 'Box', 'on');
            end
            set(handles.tAxes(1:end-1), 'XTickLabel', []);
            
            if strcmpi(pUnitType, 's')
                xlabel(handles.tAxes(end), 'Time (s)');
            else
                xlabel(handles.tAxes(end), 'Frames');
            end
        end
    end

    function trackChoice_Callback(~,~)
        str = cellstr(get(trackChoice, 'String'));
        str = str{get(trackChoice,'Value')};
        setTrackColormap(str);
        setColorbar(str);
        updateSlice();
        restoreFocus();
    end


    function unitChoice_Callback(varargin)       
        contents = cellstr(get(tplotUnitChoice,'String'));
        switch contents{get(tplotUnitChoice,'Value')}
            case 'Seconds'
                pUnitType = 's';
            case 'Frames'
                pUnitType = 'f';
        end
        updateTrack();
        restoreFocus();
    end


    function setTrackColormap(mode)
        switch mode
            case 'Category'
                trackColormap = [0 1 0; 1 1 0; 1 0.5 0; 1 0 0; 0 1 1; 0 0.5 1; 0 0 1; 0.5 0 1];
                trackColormap = trackColormap([tracks.catIdx],:);
            case 'Lifetime'
                lifetimes_f = round([tracks.lifetime_s]/data.framerate);
                df = data.movieLength-round(120/data.framerate);
                dcoord = 0.25/df;
                trackColormap = [jet(round(120/data.framerate)); (0.5:-dcoord:0.25+dcoord)' zeros(df,2)];
                trackColormap = trackColormap(lifetimes_f,:);
            case 'EAP Status'
                trackColormap = hsv2rgb([0 0 0.8; 0.55 1 0.9; 0.33 1 0.9]); % ns, slave sig., master sig.
                S = [tracks.significantSlave];
                M = [tracks.significantMaster];
                eap = ones(1,nt);
                eap(M(2,:)==1) = 3;
                eap(S(2,:)==1 & M(2,:)==0) = 2;
                trackColormap = trackColormap(eap,:);                
            case 'Object Type'
                isCCP = [tracks.isCCP];
                trackColormap = [0.8 0 0; 0 0.8 0];
                trackColormap = trackColormap(isCCP+1,:);
            case 'Random'
                trackColormap = hsv2rgb([rand(nt,1) ones(nt,2)]);
        end
    end

        
    function chlabel_Callback(~,~)
        if get(labelCheckbox, 'Value') %&& ~isRGB
            set(hChLabel, 'Visible', 'on');
        else
            set(hChLabel, 'Visible', 'off');
        end
        restoreFocus();
    end


    function statsButton_Callback(varargin)
        if ~isempty(tracks)
            plotTrackClasses([tracks.catIdx]);
        end
        restoreFocus();
    end

    %---------------------------------------------------------------------------
    % Settings window
    %---------------------------------------------------------------------------
    function trackSettings_Callback(varargin)
        
        % open window with settings panel
        tpos = get(hfig, 'Position');
        tpos = [tpos(1)+tpos(3)/2-150 tpos(2)+tpos(4)/2-75 300 255];
        pht = figure('Units', 'pixels', 'Position', tpos,...
            'PaperPositionMode', 'auto', 'Menubar', 'none', 'Toolbar', 'none',...
            'Color', get(0,'defaultUicontrolBackgroundColor'),...
            'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels',...
            'Name', 'Track display settings', 'NumberTitle', 'off', 'Resize', 'off');
        
        % Lifetime selection sliders
        b  = 195;
        uicontrol(pht, 'Style', 'text', 'String', 'Lifetimes:',...
            'Position', [5 b+35 90 20], 'HorizontalAlignment', 'left');
        
        if maxLft>minLft
            uicontrol(pht, 'Style', 'text', 'String', 'Min.:',...
                'Position', [5 b+18 35 20], 'HorizontalAlignment', 'left');
            minLftSlider = uicontrol(pht, 'Style', 'slider',...
                'Value', minVal, 'SliderStep', data.framerate/(maxLft-minLft-data.framerate)*[1 5], 'Min', minLft, 'Max', maxLft,...
                'Position', [45 b+20 200 18]);
            addlistener(handle(minLftSlider), 'Value', 'PostSet', @minSlider_Callback);
            minTxt = uicontrol(pht, 'Style', 'text', 'String', [num2str(minVal) ' s'],...
                'Position', [250 b+18 40 20], 'HorizontalAlignment', 'left');
            
            uicontrol(pht, 'Style', 'text', 'String', 'Max.:',...
                'Position', [5 b-2 35 20], 'HorizontalAlignment', 'left');
            maxLftSlider = uicontrol(pht, 'Style', 'slider',...
                'Value', maxVal, 'SliderStep', data.framerate/(maxLft-minLft-data.framerate)*[1 5], 'Min', minLft, 'Max', maxLft,...
                'Position', [45 b 200 18]);
            addlistener(handle(maxLftSlider), 'Value', 'PostSet', @maxSlider_Callback);
            maxTxt = uicontrol(pht, 'Style', 'text', 'String', [num2str(maxVal) ' s'],...
                'Position', [250 b-2 40 20], 'HorizontalAlignment', 'left');
        end
        
        % Track category selection
        b = 165;
        uicontrol(pht, 'Style', 'text', 'String', 'Track category:',...
            'Position', [5 b 165 20], 'HorizontalAlignment', 'left');
        
        % Create the button group.
        hr = uibuttongroup('Visible', 'off', 'Units', 'pixels', 'Position', [90 b+5 120 20],...
            'BorderWidth', 0);
        % Create three radio buttons in the button group.
        rSet(1) = uicontrol('Style', 'radiobutton', 'String','All',...
            'pos',[0 0 45 20],'parent',hr,'HandleVisibility','off');
        rSet(2) = uicontrol('Style','radiobutton','String','CCPs',...
            'pos',[45 0 50 20],'parent',hr,'HandleVisibility','off');
        set(hr,'SelectionChangeFcn',@radio_Callback);
        set(hr,'SelectedObject', rSet(setRadioVal));
        set(hr,'Visible','on');
        
        b = 145;
        mitCheck = uicontrol(pht, 'Style', 'checkbox', 'String', 'Max. intensity threshold:',...
            'Position', [15 b 165 20], 'HorizontalAlignment', 'left', 'Value', mitCheckVal);
        mitText = uicontrol(pht, 'Style', 'edit', 'String', num2str(maxIntT, '%.1f'),...
            'Position', [170 b 60 20], 'HorizontalAlignment', 'right');
        
        % Category selection buttons
        b = 110;
        catCheck = zeros(1,8);
        uicontrol(pht, 'Style', 'text', 'String', 'Single tracks: ',...
            'Position', [15 b+10 90 20], 'HorizontalAlignment', 'left');
        catCheck(1) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Valid',...
            'Position', [15 b 60 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(1));
        catCheck(2) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Faulty',...
            'Position', [75 b 140 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(2));
        catCheck(3) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Cut',...
            'Position', [135 b 80 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(3));
        catCheck(4) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Persistent',...
            'Position', [195 b 90 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(4));
        
        b = 75;
        uicontrol(pht, 'Style', 'text', 'String', 'Compound tracks: ',...
            'Position', [15 b+10 120 20], 'HorizontalAlignment', 'left');
        catCheck(5) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Valid',...
            'Position', [15 b 60 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(5));
        catCheck(6) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Faulty',...
            'Position', [75 b 140 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(6));
        catCheck(7) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Cut',...
            'Position', [135 b 80 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(7));
        catCheck(8) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Persistent',...
            'Position', [195 b 90 15], 'HorizontalAlignment', 'left', 'Value', catCheckVal(8));
        
        % EAP status selection buttons
        b = 35;
        eapCheck = zeros(1,3);
        uicontrol(pht, 'Style', 'text', 'String', 'EAP significance: ',...
            'Position', [5 b+10 150 20], 'HorizontalAlignment', 'left');
        eapCheck(1) = uicontrol(pht, 'Style', 'checkbox', 'String', 'Independent',...
            'Position', [5 b 110 15], 'HorizontalAlignment', 'left', 'Value', eapCheckVal(1));
        eapCheck(2) = uicontrol(pht, 'Style', 'checkbox', 'String', 'M/S',...
            'Position', [110 b 140 15], 'HorizontalAlignment', 'left', 'Value', eapCheckVal(2));
        eapCheck(3) = uicontrol(pht, 'Style', 'checkbox', 'String', 'N.S.',...
            'Position', [165 b 80 15], 'HorizontalAlignment', 'left', 'Value', eapCheckVal(3));

        
        uicontrol(pht, 'Style', 'pushbutton', 'String', 'Reset',...
            'Position', [20 5 100 20], 'HorizontalAlignment', 'left',...
            'Callback', @resetButton_Callback);
        uicontrol(pht, 'Style', 'pushbutton', 'String', 'Apply',...
            'Position', [180 5 100 20], 'HorizontalAlignment', 'left',...
            'Callback', @applyButton_Callback);
        
        if setRadioVal==2
            set(catCheck, 'enable', 'off');
        end
        
        restoreFocus();
        
        function minSlider_Callback(~, eventdata)
%             obj = get(eventdata, 'AffectedObject');
            minVal = round(eventdata.AffectedObject.Value);
            if minVal >= maxVal
                minVal = maxVal;
                set(minLftSlider, 'Value', minVal);
            end
            set(minTxt, 'String', [num2str(minVal) ' s']);
        end
        
        function maxSlider_Callback(~, eventdata)
%             obj = get(eventdata, 'AffectedObject');
            maxVal = round(eventdata.AffectedObject.Value);
            if maxVal <= minVal
                maxVal = minVal;
                set(maxLftSlider, 'Value', maxVal);
            end
            set(maxTxt, 'String', [num2str(maxVal) ' s']);            
        end
        
        function radio_Callback(~,eventdata)
            switch get(eventdata.NewValue,'String')
                case 'All' % 
                    setRadioVal = 1;
                    set(mitCheck, 'Value', mitCheckVal);
                    arrayfun(@(i) set(catCheck(i), 'Value', catCheckVal(i)), 1:8);
                    set(catCheck, 'enable', 'on');
                case 'CCPs' % CCP settings, but do not change stored values
                    setRadioVal = 2;
                    set([mitCheck catCheck(1)], 'Value', 1);
                    set(catCheck(2:end), 'Value', 0);
                    set(catCheck, 'enable', 'off');
            end
        end
        
        function resetButton_Callback(varargin)
            maxVal = maxLft;
            minVal = data.framerate*5;%minLft;
            set(maxLftSlider, 'Value', maxVal);
            set(minLftSlider, 'Value', minVal);
            setRadioVal = 1;
            set(hr,'SelectedObject', rSet(setRadioVal));
            set(mitCheck, 'Value', false, 'Enable', 'on');
            set(mitText, 'String', num2str(maxIntT, '%.1f'));
            set([catCheck eapCheck], 'Value', true, 'Enable', 'on');
        end
        
        function applyButton_Callback(varargin)
            % update track selection index
            catCheckVal = cell2mat(get(catCheck, 'Value'))==1;
            eapCheckVal = cell2mat(get(eapCheck, 'Value'))==1;
            mitCheckVal = get(mitCheck, 'Value')==1;
            maxIntT = str2double(get(mitText, 'String'));
            
            % EAP: indep: M(2,:)==1; M/S M(2,:)==0 & S(2,:)==1; n.s. S(2,:)==0
            trackIndex.selected = ismember([tracks.catIdx], find(catCheckVal)) & ...
                minVal<=[tracks.lifetime_s] & [tracks.lifetime_s]<=maxVal & maxA(1,:)>=maxIntT;
            if isfield(tracks, 'significantSlave')
                S = [tracks.significantSlave];
                M = [tracks.significantMaster];
                trackIndex.selected = trackIndex.selected & ...
                    ((eapCheckVal(1) & M(2,:)==1) | ...
                    (eapCheckVal(2) & M(2,:)==0 & S(2,:)==1) | ...
                    (eapCheckVal(3) & S(2,:)==0));
            end
            
            % update track selection
            trackIndex.current = find(trackIndex.selected, 1, 'first');
            if sum(trackIndex.selected)>1
                set(handles.trackSlider, 'Visible', 'on');
                set(handles.trackSlider, 'Min', 1);
                set(handles.trackSlider, 'Max', sum(trackIndex.selected));
                set(handles.trackSlider, 'SliderStep', [1/(sum(trackIndex.selected)-1) 0.05]);
                set(handles.trackSlider, 'Value', 1);
                set(handles.trackSlider, 'Enable', 'off');
                drawnow;
                set(handles.trackSlider, 'enable', 'on');
            else % no, or single track -> hide slider
                set(handles.trackSlider, 'Visible', 'off');
                hc = get(handles.tAxes, 'Children');
                if iscell(hc)
                    hc = [hc{:}];
                end
                set(hc, 'Visible', 'off');
                set(handles.trackIndex.label, 'String', 'Track N/A');
            end
            
            if ip.Results.LoadFrames
                updateSlice();
            end
            updateTrack();
            close(pht);
            fprintf('# tracks selected: %d\n', sum(trackIndex.selected));
        end
    end
    %---------------------------------------------------------------------------
    % End settings window
    %---------------------------------------------------------------------------
    
    
    function setColorbar(mode)        
        lfont = {'FontName', 'Helvetica', 'FontSize', 12};
        sfont = {'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'normal'};
        if ~isempty(tracks)
            switch mode
                case 'Lifetime'
                    df = 40;
                    dcoord = 0.25/df;
                    lmap = [jet(120); (0.5:-dcoord:0.25+dcoord)' zeros(df,2)];
                    imagesc(reshape(lmap, [size(lmap,1) 1 3]), 'Parent', hLegend);
                    set(hLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
                        'YTick', [1 20:20:120 160],...
                        'YTickLabel', [data.framerate 20:20:120 (nf-1)*data.framerate], sfont{:});
                    text(-0.1, 80, 'Lifetime (s)', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', hLegend, lfont{:});
                case 'Category'
                    xlabels = {' valid', ' faulty', ' cut', ' persistent',...
                        ' valid', ' faulty', ' cut', ' persistent'};
                    lmap = [0 1 0; 1 1 0; 1 0.5 0; 1 0 0; 0 1 1; 0 0.5 1; 0 0 1; 0.5 0 1];
                    imagesc(reshape(lmap, [size(lmap,1) 1 3]), 'Parent', hLegend);
                    set(hLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
                        'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0]);
                    text(-.1, 2.5, 'Single', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', hLegend, lfont{:});
                    text(-.1, 6.5, 'Compound', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Parent', hLegend, lfont{:});
                case 'EAP Status'
                    xlabels = {' N.S.', ' Signif. M/S', ' Signif. indep.'};
                    lmap = hsv2rgb([0 0 0.8; 0.55 1 0.9; 0.33 1 0.9]); % ns, slave sig., master sig.
                    imagesc(reshape(lmap, [size(lmap,1) 1 3]), 'Parent', hLegend);
                    set(hLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
                        'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0]);
                case 'Object Type'
                    xlabels = {' Diff. lim.', ' Other'};
                    lmap = [0 0.8 0; 0.8 0 0];
                    imagesc(reshape(lmap, [size(lmap,1) 1 3]), 'Parent', hLegend);
                    set(hLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
                        'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0]);
                otherwise
                    cla(hLegend);
                    set(hLegend, 'Visible', 'off');
            end
        end
    end


    function montageButton_Callback(varargin)
        
        % Creates a montage based on the master track
        if ~isempty(trackIndex.current)
            fprintf('Generating montage...');
            if get(montageAlignCheckbox, 'Value')
                ref = 'Track';
            else
                ref = 'Frame';
            end
            [istack, xa, ya] = getTrackStack(trackIndex.current, 6, ref);
            plotTrackMontage(tracks(trackIndex.current), istack, xa, ya, 'Labels', data.markers,...
                'ShowMarkers', get(montageMarkerCheckbox, 'Value')==1,...
                'ShowDetection', get(montageDetectionCheckbox, 'Value')==1);
            fprintf(' done.\n');
        else
            fprintf('Cannot create montage: no track selected.\n');
        end
        restoreFocus();
    end


    function [tstack, xa, ya] = getTrackStack(t, w, reference)
        
        sigma = frameInfo(1).s;
        w = ceil(w*sigma);
        
        % coordinate matrices
        x0 = tracks(t).x;
        y0 = tracks(t).y;
        
        % start and end buffer sizes
        if ~isempty(tracks(t).startBuffer)
            sb = numel(tracks(t).startBuffer.t);
            x0 = [tracks(t).startBuffer.x x0];
            y0 = [tracks(t).startBuffer.y y0];
        else
            sb = 0;
        end
        if ~isempty(tracks(t).endBuffer)
            eb = numel(tracks(t).endBuffer.t);
            x0 = [x0 tracks(t).endBuffer.x];
            y0 = [y0 tracks(t).endBuffer.y];
        else
            eb = 0;
        end
        
        % frame index
        tfi = tracks(t).start-sb:tracks(t).end+eb;
        tnf = length(tfi);
        
        
        if tracks(t).nSeg==1 && strcmpi(reference, 'track') % align frames to track
            xi = round(x0(handles.mCh,:));
            yi = round(y0(handles.mCh,:));
            % ensure that window falls within frame bounds
            x0 = xi - min([xi-1 w]);
            x1 = xi + min([nx-xi w]);
            y0 = yi - min([yi-1 w]);
            y1 = yi + min([ny-yi w]);
            % axes for each frame
            xa = arrayfun(@(i) x0(i):x1(i), 1:tnf, 'unif', 0);
            ya = arrayfun(@(i) y0(i):y1(i), 1:tnf, 'unif', 0);
        else
            % window around track mean
            mu_x = round(nanmean(x0,2));
            mu_y = round(nanmean(y0,2));
            x0 = max(1, min(mu_x)-w);
            x1 = min(data.imagesize(2), max(mu_x)+w);
            y0 = max(1, min(mu_y)-w);
            y1 = min(data.imagesize(1), max(mu_y)+w);
            xa = repmat({x0:x1}, [tnf 1]);
            ya = repmat({y0:y1}, [tnf 1]);
        end
        
        tstack = cell(nCh,tnf);
        for ci = 1:nCh
            for k = 1:tnf
                tstack{ci,k} = stack{ci}(ya{k}, xa{k}, tfi(k));
            end
        end
    end


    function trackSelectButton_Callback(varargin)
        if get(trackSelectButton, 'Value')==1
            set(hfig, 'Pointer', 'crosshair');
        else
            set(hfig, 'Pointer', 'arrow');
        end
        restoreFocus();
    end


    function printButton_Callback(varargin)
        fprintf('Printing figures ...');
        
        % Tracks
        if ~isempty(tracks)
            for ch = 1:nCh
                plotTrack(data, tracks(trackIndex.current), ch,...
                    'FileName', ['track_' num2str(trackIndex.loaded(trackIndex.current)) '_ch' num2str(ch) '.eps'],...
                    'Visible', 'off', 'DisplayMode', 'Print');
            end
            
            if get(montageAlignCheckbox, 'Value')
                ref = 'Track';
            else
                ref = 'Frame';
            end
            [tstack, xa, ya] = getTrackStack(trackIndex.current, 6, ref);
            fpath = [data.source 'Figures' filesep 'track_' num2str(trackIndex.loaded(trackIndex.current)) '_montage.eps'];
                plotTrackMontage(tracks(trackIndex.current), tstack, xa, ya, 'Labels', data.markers,...
                    'Visible', 'off', 'epsPath', fpath,...
                    'ShowMarkers', get(montageMarkerCheckbox, 'Value')==1,...
                    'ShowDetection', get(montageDetectionCheckbox, 'Value')==1);
        end
        
        % Frames
        if strcmp(displayType, 'RGB')
            maxCh = 1;
            chLabel = {'RGB'};
        else
            maxCh = nCh;
            chLabel = arrayfun(@(i) ['ch' num2str(i)], 1:nCh, 'unif', 0);
        end
                
        f0 = figure('PaperPositionMode', 'auto', 'Position', [20 20 nx ny], 'Visible', 'off',...
            'DefaultLineLineSmoothing', 'on', 'DefaultPatchLineSmoothing', 'on');
        colormap(gray(256));
        fpath = [data.source 'Figures' filesep];
        for ci = 1:maxCh
            h0 = copyobj(handles.fAxes(ci,1),f0);
            hx = findobj(h0, 'DisplayName', 'FrameMarker');
            delete(hx);
            hx = findobj(h0, 'DisplayName', 'TrackMarker');
            delete(hx);
            hx = findobj(h0, 'LineStyle', '-');
            set(hx, 'LineWidth', 1);
            
            hx = findobj(h0, 'Marker', 'o');
            set(hx, 'MarkerSize', 9, 'LineWidth', 1);
            
            set(h0, 'Position', [0 0 1 1]);
            print(f0, '-depsc2', '-loose', [fpath 'frame_' num2str(fidx) '_' chLabel{ci} '.eps']);
            %print(f0, '-dpng', '-loose', [fpath 'frame_' num2str(fidx) '_' chLabel{ci} '.png']);
            delete(h0);
        end
        close(f0);
        
        fprintf([' done. Figures saved in ' getShortPath(data) 'Figures.\n']);
        restoreFocus();
    end


    function movieButton_Callback(varargin)
        
        fopts = {'Visible', 'off', 'Position', [20 20 nx ny],...
            'InvertHardcopy', 'off', 'PaperUnits', 'Points', 'PaperSize', [nx ny],...
            'PaperPosition', [0 0 nx ny], 'PaperPositionMode', 'auto',...
            'DefaultLineLineSmoothing','on', 'DefaultPatchLineSmoothing','on'};
        
        if strcmp(displayType, 'RGB')
            maxCh = 1;
        else
            maxCh = nCh;
        end
        
        mpath = [data.source 'Movies' filesep];
        fpath = [mpath 'Frames' filesep];
        [~,~] = mkdir(mpath);
        [~,~] = mkdir(fpath);
        
        fmt = ['%0' num2str(ceil(log10(nf))) 'd'];
        
        f0 = figure(fopts{:});
        colormap(gray(256));
        ha = axes('Position', [0 0 1 1]);
        if ~strcmp(computer('arch'),'win64')
            movieFrames(1:nf) = getframe(f0); % Pre-allocate movie struct array
        end
        for ci = 1:maxCh
            fprintf('Generating movie frames:     ');
            for fi = 1:nf
                switch displayType
                    case 'raw'
                        imagesc(stack{ci}(:,:,fi), 'Parent', ha);
                    case 'mask'
                        if ci==1
                            imagesc(rgbOverlay(stack{1}(:,:,fi), dmask(:,:,fi), [1 0 0], dRange{1}), 'Parent', ha);
                        else
                            imagesc(stack{ci}(:,:,fi), 'Parent', ha);
                        end
                    case 'RGB'
                        for c2 = 1:nCh
                            rframe(:,:,idxRGB(c2)) = uint8(scaleContrast(double(stack{c2}(:,:,fi)), dRange{c2}));
                        end
                        imagesc(rframe, 'Parent', ha);
                end
                hold(ha, 'on');
                caxis(ha, dRange{ci});
                
                if ~isempty(tracks) && fi~=1 && get(trackCheckbox, 'Value')
                    vidx = ~isnan(X(fi,:));
                    
                    set(ha, 'ColorOrder', trackColormap(trackIndex.segmentIndex(vidx),:));
                    
                    plot(ha, X(1:fi,vidx), Y(1:fi,vidx), 'HitTest', 'off');
                    if get(gapCheckbox, 'Value')
                        hpg = plot(ha, X(fi,vidx & G(fi,:)), Y(fi,vidx & G(fi,:)), 'o', 'Color', 'w', 'MarkerSize', 6, 'LineWidth', 1);
                    end
                    if get(trackEventCheckbox, 'Value')
                        % Births
                        bcoord = arrayfun(@(i) [i.x(1,1) i.y(1,1)], tracks(trackStarts==fi), 'unif', 0);
                        bcoord = vertcat(bcoord{:});
                        if~isempty(bcoord)
                            plot(ha, bcoord(:,1), bcoord(:,2), '*', 'Color', 'g', 'MarkerSize', 8, 'LineWidth', 1);
                        end
                        
                        % Deaths
                        dcoord = arrayfun(@(i) [i.x(1,1) i.y(1,1)], tracks(trackEnds==fi), 'unif', 0);
                        dcoord = vertcat(dcoord{:});
                        if ~isempty(dcoord)
                            plot(ha, dcoord(:,1), dcoord(:,2), 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 1);
                        end
                    end
                    
                end
                if ~isempty(tracks) && get(eapCheckbox, 'Value') && ci>1
                    sel = fvec==fi & mvec(ci,:)==1;
                    plot(ha, xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', hsv2rgb([1/3 1 0.9]), 'MarkerSize', 8);
                    sel = fvec==fi & mvec(ci,:)==0 & svec(ci,:)==1;
                    plot(ha, xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', hsv2rgb([0.55 1 0.9]), 'MarkerSize', 8);
                    sel = fvec==fi & mvec(ci,:)==0 & svec(ci,:)==0;
                    plot(ha, xvec(ci,sel), yvec(ci,sel) , 'o', 'Color', 0.8*[1 1 1], 'MarkerSize', 8);
                end
                
                if get(detectionCheckbox, 'Value') && ~isempty(frameInfo)
                    isPSF = frameInfo(fi).isPSF(1,:)==1;
                    if any(isPSF)
                        plot(ha, frameInfo(fi).x(1,isPSF), frameInfo(fi).y(1,isPSF), 'o', 'Color', [0 0.6 0], 'MarkerSize', 8);
                    end
                    if any(~isPSF)
                        plot(ha, frameInfo(fi).x(1,~isPSF), frameInfo(fi).y(1,~isPSF), 'o', 'Color', [0.6 0 0], 'MarkerSize', 8);
                    end
                end
                
                axis(ha, 'off');
                
                if ~strcmp(computer('arch'),'win64')

                    print(f0, '-dpng', '-loose', ['-r' num2str(1*72)], [fpath 'frame' num2str(fi, fmt) '_ch' num2str(ci) '.png']);
                    %print(h, '-djpeg100', '-loose', ['-r' num2str(zoom*72)], [fpath 'frame' num2str(f, fmt) ext]);
                else
                    movieFrames(fi) = getframe(f0);
                end
                
                cla(ha);
                fprintf('\b\b\b\b%3d%%', round(100*fi/nf));

            end
            fprintf('\n');
        end
        fprintf(['Frames saved to ' getShortPath(data) 'Movies' filesep 'Frames.\n']);
        close(f0);
        
        if isunix
            % side-by-side frame arrangement
            for fi = 1:nf
                % channel frames
                cpath = arrayfun(@(ci) [fpath 'frame' num2str(fi, fmt) '_ch' num2str(ci) '.png '], 1:nCh, 'unif', 0);
                fname = [fpath 'montage' num2str(fi, fmt) '.png'];
                
                cmd = ['export DYLD_LIBRARY_PATH=""; montage -geometry +3+3+0+0 -background "rgb(255,255,255)" '...
                    [cpath{:}] ' -compress lzw ' fname];
                system(cmd);
                cmd = ['export DYLD_LIBRARY_PATH=""; convert ' fname ' -shave 3x3 -depth 8 ' fname];
                system(cmd);
            end
        end        
        % Generate movie, if on a unix system with ffmpeg
        if isunix && ~system('which ffmpeg >/dev/null 2>&1')
            fprintf('Generating movie ... ');
            %fr = num2str(framerate);
            fr = num2str(15);           
            
            %cmd = ['ffmpeg -y -r ' fr ' -i ' fpath 'frame' fmt '_ch' num2str(1) '.png' ' -vf "scale=' num2str(2*floor(nx/2)) ':' num2str(2*floor(ny/2))...
            %    '" -c:v libx264 -crf 22 -pix_fmt yuv420p ' mpath 'Movie_ch1.mp4'];
            cmd = ['ffmpeg -y -r ' fr ' -i ' fpath 'montage' fmt '.png' ' -vf "scale=' num2str(2*floor((nCh*nx+(nCh-1)*6)/2)) ':' num2str(2*floor(ny/2))...
                '" -c:v libx264 -crf 22 -pix_fmt yuv420p ' mpath getCellDir(data) '.mp4'];
            system(cmd);
            
            fprintf(' done.\n');
        else
           fprintf('A unix system with ffmpeg installed is required to generate ffmeg movies automatically.\n');
           disp('Creating AVI movie instead of ffmeg unix based movie!')
           v = VideoWriter([fpath 'MovieCME.avi']);
           open(v);
           writeVideo(v, movieFrames);
           close(v);         
        end
        restoreFocus();
    end


    function setupTrackAxes(nAxes)
        opts = {'Parent', hfig, 'Units', 'pixels', 'Box', 'on'};
        axPos = getTrackAxesPositions(hfig, nAxes);
        if isfield(handles, 'tAxes') && ~isempty(handles.tAxes)
            delete(handles.tAxes);
        end
        handles.tAxes = zeros(1,nAxes);
        for k = 1:nAxes
            handles.tAxes(k) = axes(opts{:}, 'Position', axPos{k});
        end
        setappdata(hfig, 'handles', handles);
    end


    function setupFrameAxes(nPanels)
        if nargin<1
            nPanels = nCh;
        end        
        opts = {'Parent', hfig, 'Units', 'pixels', 'BorderType', 'none'};
        axPos = getFrameAxesPositions(hfig, nPanels);
        if isfield(handles, 'fAxes') && ~isempty(handles.fAxes)
            delete(handles.fAxes);
        end
        if isfield(handles, 'fPanels')
            delete(handles.fPanels);
        end
        handles.fPanels = zeros(1,nPanels);
        handles.fAxes = zeros(nPanels,3);
        hLegend = zeros(1,nPanels);

        for k = 1:nPanels
            handles.fPanels(k) = uipanel(opts{:}, 'Position', axPos{k});
            % add stackviewer
            [handles.fAxes(k,:), hLegend(k)] = setupStackViewer(handles.fPanels(k), [nx ny min(nf,  max(nx,ny)/3)], k==1);
        end
        hLegend = hLegend(1);
        setappdata(hfig, 'handles', handles);
    end
end



function figResize(src,~)
handles = getappdata(src, 'handles');

pos = get(src, 'Position');

set(handles.frameLabel, 'Position', [20 pos(4)-20, 100 15]);

% tracks
set(handles.tplotPanel, 'Position', [pos(3)-585 5 210 70]);
set(handles.montagePanel, 'Position', [pos(3)-370 5 220 70]);
set(handles.outputPanel, 'Position', [pos(3)-145 5 140 70]);

% spacers:
lspace = 10;
rspace = 400;

set(handles.frameSlider, 'Position', [lspace 77 pos(3)-rspace-lspace 18]);

nPanels = numel(handles.fPanels);
axPos = getFrameAxesPositions(src, nPanels);
for i = 1:nPanels
    set(handles.fPanels(i), 'Position', axPos{i});
end

nAxes = numel(handles.tAxes);
h_tot = pos(4) - 140;
axPos = getTrackAxesPositions(src, nAxes);
for i = 1:nAxes
    set(handles.tAxes(i), 'Position', axPos{i});
end
set(handles.trackIndex.label, 'Position', [pos(3)-160 pos(4)-18 120 15]);
set(handles.trackSlider, 'Position', [pos(3)-24 120 18 h_tot]);

end



function pos = getTrackAxesPositions(hfig, nAxes)
pos = get(hfig, 'Position');

spacer = 15;
w = 320;
h_tot = pos(4) - 140;
h = min((h_tot-(nAxes-1)*spacer)/nAxes, 200);
dx = pos(3)-w-30;

pos = cell(1,nAxes);
pos{1} = [dx 120+(h_tot-h) w h];
for i = 2:nAxes
    pos{i} = [dx 120+(h_tot-i*h-(i-1)*spacer) w h];
end
end



function pos = getFrameAxesPositions(hfig, nPanels)
pos = get(hfig, 'Position');

% spacers:
tspace = 20;
bspace = 100;
lspace = 10;
rspace = 400;
spacer = 10; % space between panels

width = pos(3) - rspace - lspace;
height = pos(4) - bspace - tspace;

pos = cell(1,nPanels);
switch nPanels
    case 1
        pos{1} = [lspace bspace width height];
    case 2
        width = (width-spacer)/2;
        pos{1} = [lspace bspace width height];
        pos{2} = [lspace+width+spacer bspace width height];
    case 3
        width = (width-spacer)/2;
        height = (height-spacer)/2;
        pos{1} = [lspace bspace+spacer+height width height]; % top left
        pos{2} = [lspace+width+spacer bspace+height+spacer width height]; % top right
        pos{3} = [lspace bspace width height]; % bottom left
    case 4
        width = (width-spacer)/2;
        height = (height-spacer)/2;
        pos{1} = [lspace bspace+spacer+height width height]; % top left
        pos{2} = [lspace+width+spacer bspace+height+spacer width height]; % top right
        pos{3} = [lspace bspace width height]; % bottom left
        pos{4} = [lspace+width+spacer bspace width height]; % bottom right
end
end



function [ha, hl] = setupStackViewer(hf, dims, addLegend)
if nargin<3
    addLegend = false;
end

[apos, lpos] = getStackViewerPositions(hf, dims);
ha(1) = axes('Position', apos{1}, 'Parent', hf);
ha(2) = axes('Position', apos{2}, 'Parent', hf); % bottom left
ha(3) = axes('Position', apos{3}, 'Parent', hf);
if addLegend
    hl = axes('Position', lpos, 'Parent', hf);
else
    hl = NaN;
end

set(hf, 'ResizeFcn', @pResize);

    function pResize(~,~)
        [apos, lpos] = getStackViewerPositions(hf, dims);
        set(ha(1), 'Position', apos{1});
        set(ha(2), 'Position', apos{2});
        set(ha(3), 'Position', apos{3});
        if ~isnan(hl)
            set(hl, 'Position', lpos);
        end
    end
end



function [apos, lpos] = getStackViewerPositions(hf, dims)

spc = 6; % spacer, fixed [pixels]

nx = dims(1);
ny = dims(2);
nz = dims(3);
pos = get(hf, 'Position');
w = pos(3);
h = pos(4);

% normalized axes dimensions
fx = (w-spc)/(nx+nz);
fy = (h-spc)/(ny+nz);
f = min(fx,fy);
h = (ny+nz)*f+spc; % [pixels]
w = (nx+nz)*f+spc;

rxy = pos(3)/pos(4);
dx = spc/pos(3);
dy = spc/pos(4);

apos = cell(1,3);
if rxy > w/h % figure is too wide
    f0 = w/h / rxy;
    left = (1-f0)/2;
    apos{1} = [left+(f0*nz*f)/w+dx 0 f0*f*nx/w f*ny/h];
    apos{2} = [left 0 f0*f*nz/w f*ny/h];
    apos{3} = [left+(f0*nz*f)/w+dx (ny*f)/h+dy f0*f*nx/w f*nz/h];
else
    f0 = h/w * rxy;
    left = 0;
    apos{1} = [(nz*f)/w+dx 1-f0 f*nx/w f0*f*ny/h];
    apos{2} = [0 1-f0 f*nz/w f0*f*ny/h];
    apos{3} = [(nz*f)/w+dx 1-f0+(f0*ny*f)/h+dy f*nx/w f0*f*nz/h];
end
lpos = apos{3};
lpos([1 3]) = [left+15/pos(3) 15/pos(3)];
end

