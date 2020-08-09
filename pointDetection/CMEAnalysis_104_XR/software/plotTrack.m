% plotTrack(data, tracks, trackIdx, ch, varargin)
%
% Inputs:   data : data structure
%          track : track structure
%            {ch} : channel number. Default: 1
%
% Options
%      'Visible' : {'on'} | 'off' toggles figure visibility
%       'Handle' : axis handle (for plotting from within GUI)
%  'DisplayMode' : {'screen'}|'print'
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

% Francois Aguet, March 9 2011 (Last modified: 01/31/2012)

function ha = plotTrack(data, track, varargin)

%======================================
% Parse inputs, set defaults
%======================================
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('track', @isstruct);
ip.addOptional('ch', [], @isnumeric);
ip.addParamValue('Visible', 'on', @(x) any(strcmpi(x, {'on', 'off'})));
ip.addParamValue('FileName', [], @ischar);
ip.addParamValue('Handle', []);
ip.addParamValue('Legend', 'hide', @(x) any(strcmpi(x, {'show','hide'})));
ip.addParamValue('Background', 'on', @(x) any(strcmpi(x, {'on', 'off'})));
ip.addParamValue('BackgroundValue', 'zero', @(x) any(strcmpi(x, {'zero', 'data'})));
ip.addParamValue('Hues', []);
ip.addParamValue('TrackColor', []);
ip.addParamValue('Time', 'Track', @(x) any(strcmpi(x, {'Movie', 'Track'})) || isscalar(x));
ip.addParamValue('XTick', []);
ip.addParamValue('YTick', []);
ip.addParamValue('XLim', []);
ip.addParamValue('DisplayMode', 'Screen', @(x) any(strcmpi(x, {'Print', 'Screen'})));
ip.addParamValue('OverlayBackground', false, @islogical);
ip.addParamValue('MarkerSizes', [21 7 2]);
ip.addParamValue('PlotBuffers', true, @islogical);
ip.addParamValue('LineWidth', 1);
ip.addParamValue('BackgroundConfidence', [], @isvector);
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.parse(data, track, varargin{:});


hues = ip.Results.Hues;
if isempty(hues)
    hues = getFluorophoreHues(data.markers);
end
trackColor = ip.Results.TrackColor;

mCh = find(strcmp(data.channels, data.source));
ch = ip.Results.ch;
if isempty(ch)
    ch = mCh;
end

if strcmpi(ip.Results.Time, 'Track')
    dt = track.t(1); % first segment always contains first time point
elseif strcmpi(ip.Results.Time, 'Movie')
    dt = 0;
else
    dt = ip.Results.Time;
end

% Flags
hasStartBuffer = ~isempty(track.startBuffer);
hasEndBuffer = ~isempty(track.endBuffer);

XLim = ip.Results.XLim;

b0 = 2*data.framerate;
if hasStartBuffer
    nbs = numel(track.startBuffer.t);
else
    nbs = 0;
end
if hasEndBuffer
    nbe = numel(track.endBuffer.t);
else
    nbe = 0;
end

bt = track.end-track.start; % #frames - 1

% w = 2*b0+nbs+nbe+bt;
% u = 10;
if isempty(XLim)
    XLim = track.t(1)-dt + [-(nbs+b0) bt+nbe+b0]*data.framerate;
end


% Setup figure window
if ~isempty(ip.Results.Handle)
    ha = ip.Results.Handle;
else
    hfig = figure('Visible', ip.Results.Visible, 'PaperPositionMode', 'auto');
    ha = axes;
end

% Color definitions
if isempty(trackColor)
    trackColor = hsv2rgb([hues(ch) 1 0.8]);
end
fillLight = hsv2rgb([hues(ch) 0.4 1]);
fillDark = hsv2rgb([hues(ch) 0.2 1]);
fillLightBuffer = hsv2rgb([hues(ch) 0.4 0.85]);
fillDarkBuffer = hsv2rgb([hues(ch) 0.2 0.85]);

kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1);

% Plot track
lh = NaN(1,15);
for s = 1:track.nSeg
    
    A = track.A(ch,:);
    c = track.c(ch,:);
    if strcmpi(ip.Results.BackgroundValue, 'zero')
        bgcorr = nanmean(c);
        c = c-bgcorr;
    else
        bgcorr = 0;
    end
    sigma_r = track.sigma_r(ch,:);
    t = track.t - dt;
    
    % alpha = 0.05 level
    if strcmpi(ip.Results.Background, 'on')
        lh(1) = fill([t t(end:-1:1)], [c c(end:-1:1)+kLevel*sigma_r(end:-1:1)],...
            fillDark, 'EdgeColor', 'none', 'Parent', ha);
        hold(ha, 'on');
    end
    
    gapIdx = find(track.gapVect~=0);
    
    % plot amplitude std.
    sigma_a = track.A_pstd(ch,:);
    
    rev = c+A-sigma_a;
    lh(2) = fill([t t(end:-1:1)], [c+A+sigma_a rev(end:-1:1)],...
        fillLight, 'EdgeColor', 'none', 'Parent', ha);
    hold(ha, 'on'); % in case background was not plotted
    
    % plot track
    ampl = A+c;
    if ch==mCh
        ampl(gapIdx) = NaN;
    end
    lh(3) = plot(ha, t, ampl, '.-', 'Color', trackColor, 'LineWidth', ip.Results.LineWidth);
    
    % plot gaps separately
    if ch==mCh
        ampl = A+c;
        if ~isempty(gapIdx)
            % dashed line everywhere
            lh(4) = plot(ha, t, ampl, '-', 'Color', trackColor, 'LineWidth', ip.Results.LineWidth);
            % plot gaps as white disks
            lh(5) = plot(ha, t(gapIdx), A(gapIdx)+c(gapIdx), 'o', 'Color', trackColor, 'MarkerFaceColor', 'w', 'LineWidth', ip.Results.LineWidth);
        end
    end
    
    % plot background level
    if strcmpi(ip.Results.Background, 'on')
        lh(6) = plot(ha, t, c, '-', 'Color', trackColor);
    end    
    if ip.Results.OverlayBackground
        lh(13) = plot(ha, t, c+kLevel*sigma_r, '-', 'Color', trackColor);
    end
end

% Plot start buffer
if hasStartBuffer && ip.Results.PlotBuffers
    A = [track.startBuffer.A(ch,:) track.A(ch,1)];
    c = [track.startBuffer.c(ch,:) track.c(ch,1)]-bgcorr;
    
    sigma_a = [track.startBuffer.A_pstd(ch,:) track.A_pstd(ch,1)];
    sigma_r = [track.startBuffer.sigma_r(ch,:) track.sigma_r(ch,1)];
    t = [track.startBuffer.t track.t(1)] - dt;
    
    if strcmpi(ip.Results.Background, 'on')
        lh(12) = fill([t t(end:-1:1)], [c c(end:-1:1)+kLevel*sigma_r(end:-1:1)],...
            fillDarkBuffer, 'EdgeColor', 'none', 'Parent', ha);
    end

    rev = c+A-sigma_a;
    fill([t t(end:-1:1)], [c+A+sigma_a rev(end:-1:1)],...
        fillLightBuffer, 'EdgeColor', 'none', 'Parent', ha);
    
    lh(7) = plot(ha, t, A+c, '.-', 'Color', trackColor, 'LineWidth', ip.Results.LineWidth);
    if strcmpi(ip.Results.Background, 'on')
        lh(8) = plot(ha, t, c, '-', 'Color', trackColor);
    end
    if ip.Results.OverlayBackground
        lh(14) = plot(ha, t, c+kLevel*sigma_r, '-', 'Color', trackColor);
    end
end

% Plot end buffer
if hasEndBuffer && ip.Results.PlotBuffers
    A = [track.A(ch,end) track.endBuffer.A(ch,:)];
    c = [track.c(ch,end) track.endBuffer.c(ch,:)]-bgcorr;
    
    sigma_a = [track.A_pstd(ch,end) track.endBuffer.A_pstd(ch,:)];
    sigma_r = [track.sigma_r(ch,end) track.endBuffer.sigma_r(ch,:)];
    t = [track.t(end) track.endBuffer.t] - dt;
    
    if strcmpi(ip.Results.Background, 'on')
        fill([t t(end:-1:1)], [c c(end:-1:1)+kLevel*sigma_r(end:-1:1)],...
            fillDarkBuffer, 'EdgeColor', 'none', 'Parent', ha);
    end
    
    rev = c+A-sigma_a;
    fill([t t(end:-1:1)], [c+A+sigma_a rev(end:-1:1)],...
        fillLightBuffer, 'EdgeColor', 'none', 'Parent', ha);
    
    lh(9) = plot(ha, t, A+c, '.-', 'Color', trackColor, 'LineWidth', ip.Results.LineWidth);
    if strcmpi(ip.Results.Background, 'on')
        lh(10) = plot(ha, t, c, '-', 'Color', trackColor);
    end
    if ip.Results.OverlayBackground
        lh(15) = plot(ha, t, c+kLevel*sigma_r, '-', 'Color', trackColor);
    end
end

if ~isempty(ip.Results.BackgroundConfidence)
    c = track.c(ch,:);
    t = track.t;
    if ~isempty(track.startBuffer)
        c = [track.startBuffer.c(ch,:) c];
        t = [track.startBuffer.t t];
    end
    if ~isempty(track.endBuffer)
        c = [c track.endBuffer.c(ch,:)];
        t = [t track.endBuffer.t];
    end
    c = c-mean(track.c(ch,:));
    t = t-dt;%track.t(1)-dt;
    plot(ha, t, c+ip.Results.BackgroundConfidence, 'k-');
end

% legend
if strcmpi(ip.Results.Legend, 'show')
    lh(11) = plot(-20:-10, rand(1,11), 'o--', 'MarkerSize', 7, 'LineWidth', ip.Results.LineWidth, 'Color', trackColor, 'MarkerFaceColor', 'w');
    l = legend(lh([3 2 11 6 1 12 8]), 'Intensity', 'Intensity uncertainty', 'Gap', 'Background intensity', 'Significance level (\alpha = 0.95)', 'Buffer', 'Buffer intensity', 'Location', 'NorthEast');
end



if ~isempty(ip.Results.XTick)
    XTick = ip.Results.XTick;
    set(ha, 'XTick', XTick);
end
set(ha, 'XLim', XLim);

if ~isempty(ip.Results.YTick)
    YTick = ip.Results.YTick;
    YLim = [YTick(1) YTick(end)];
    YTick(YTick<0) = [];
    set(ha, 'YTick', YTick, 'YLim', YLim, 'Layer', 'top');
else
    YTick = get(ha, 'YTick');
    YTick(YTick<0) = [];
    set(ha, 'YTick', YTick, 'Layer', 'top');
end



box off;


% Bigger fonts, line widths etc
if isempty(ip.Results.Handle)
    if strcmpi(ip.Results.Legend, 'show')
        set(l);
    end
    
    set(gca, 'LineWidth', ip.Results.LineWidth, 'TickDir', 'out');
    xlabel('Time (s)')
    ylabel('Intensity (A.U.)');
end

if strcmpi(ip.Results.DisplayMode, 'Print')
    for k = lh([3 4 6 7:10 13:15])
        if ~isnan(k)
            set(k, 'LineWidth', ip.Results.MarkerSizes(3));
        end
    end
    
    for k = lh([3 7 9])
        if ~isnan(k)
            set(k, 'MarkerSize', ip.Results.MarkerSizes(1));
        end
    end
    
    if ~isnan(lh(5)) % gaps
        set(lh(5), 'MarkerSize', ip.Results.MarkerSizes(2), 'LineWidth', ip.Results.MarkerSizes(3));
    end
end
    


if ~isempty(ip.Results.FileName)
    fpath = [data.source 'Figures' filesep];
    if ~(exist(fpath, 'dir')==7)
        mkdir(fpath);
    end
    print(hfig, '-depsc2', '-r300', [fpath ip.Results.FileName]);
end

if strcmp(ip.Results.Visible, 'off')
    close(hfig);
end
