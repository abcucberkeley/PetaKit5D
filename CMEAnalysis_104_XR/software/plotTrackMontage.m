% plotTrackMontage(track, stack) generates a montage of the individual frames of a track
%
% Inputs:   
%    track : track data structure
%    stack : cell array of image frames corresponding to the track
%
% Example:
% [stack, xa, ya] = getTrackStack(data, track, 'WindowWidth', 6);
% plotTrackMontage(track, stack, xa, ya, 'Labels', data.markers);
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

% François Aguet (last modified 02/07/2012)

function plotTrackMontage(track, trackStack, varargin)

[nc, nf] = size(trackStack);

%======================================
% Parse inputs, set defaults
%======================================
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('track', @isstruct);
ip.addRequired('trackStack', @iscell);
ip.addOptional('xa', []);
ip.addOptional('ya', []);
ip.addParamValue('FontName', 'Helvetica', @ischar);
ip.addParamValue('FontSize', 14, @isscalar);
ip.addParamValue('Visible', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('epsPath', []);
ip.addParamValue('Labels', [], @(x) length(x)==nc && ~any(isnan(name2wavelength(x))));
ip.addParamValue('Width', 600, @isscalar);
ip.addParamValue('Mode', []);
ip.addParamValue('FramesPerRow', 20, @isscalar);
ip.addParamValue('ShowDetection', false, @islogical);
ip.addParamValue('ShowMarkers', false, @islogical);
ip.addParamValue('DynamicRange', []);
ip.parse(track, trackStack, varargin{:});
width = ip.Results.Width;
labels = ip.Results.Labels;

xa = ip.Results.xa;
ya = ip.Results.ya;
if isempty(xa)
    w = (size(trackStack{1,1},2)-1)/2;
    xa = -w:w;
end
if isempty(ya)
    w = (size(trackStack{1,1},1)-1)/2;
    ya = -w:w;
end



if ~isempty(labels)
    rgbColors = arrayfun(@(x) hsv2rgb([x 1 0.9]), getFluorophoreHues(labels), 'UniformOutput', false);
else
    rgbColors = [];
end

% dynamic range for each channel
if ~isempty(track.startBuffer)
    sb = numel(track.startBuffer.t);
else
    sb = 0;
end

% dynamic range of the entire stack
dRange = ip.Results.DynamicRange;
if isempty(dRange)
    dRange = zeros(nc,2);
    for c = 1:nc
        cstack = cat(3, trackStack{c,:});
        dRange(c,:) = [min(cstack(:)) max(cstack(:))];
    end
end

idx = find(strcmpi(varargin, 'Mode'));
if ~isempty(idx)
    rgbIdx = getRGBindex(labels);
    [ny,nx] = size(trackStack{1,1});
    rgbCells = cell(3,nf);
    rgbCells(rgbIdx,:) = trackStack;
    zcells = squeeze(num2cell(zeros(ny,nx,nf,'uint8'), [1 2]));
    for k = setdiff(1:3, rgbIdx)
        rgbCells(k,:) = zcells;
    end
    rgbCells = arrayfun(@(x) cat(3,rgbCells{:,x}), 1:nf, 'UniformOutput', false);
    
    switch varargin{idx+1}
        case 'color'
            for c = 1:nc
                trackStack(c,:) = cellfun(@(x) circshift(cat(3, x, zeros(ny,nx,2)), [0 0 rgbIdx(c)-1]), trackStack(c,:), 'UniformOutput', false);
            end
            rgbColors = repmat({zeros(1,3)}, [1 nc]);
        case 'RGB+color'
            for c = 1:nc
                trackStack(c,:) = cellfun(@(x) circshift(cat(3, x, zeros(ny,nx,2)), [0 0 rgbIdx(c)-1]), trackStack(c,:), 'UniformOutput', false);
            end
            trackStack = cat(1, rgbCells, trackStack);
            nc = nc + 1;
            labels = [' ', labels];
            rgbColors = [zeros(1,3) rgbColors];
        case 'RGB+gray'
            trackStack = cat(1, rgbCells, trackStack);
            nc = nc + 1;
            labels = [' ', labels];
            rgbColors = [zeros(1,3) rgbColors];
        case 'RGB'
            trackStack = rgbCells;
            nc = 1;
            labels = [];
    end
end

% nx is now the number of frames displayed per row
nx = ip.Results.FramesPerRow;


% display with fixed #frames/width
[wxi, dxi, dci, nr, height] = getProportions(width, nx, nf, nc);



hf = figure('Visible', 'off', 'PaperPositionMode', 'auto', 'Position', [50, 200, width, height], 'Color', 'w');
ha = axes('Position', [0 0 1 1], 'XLim', [0 width], 'YLim', [0 height]);
set(hf, 'Units', 'pixels');

% compute max. text length
textWidth = zeros(1,nc);
if ~isempty(labels)
    for c = 1:nc
        ht = text(0, 0, labels{c}, 'FontName', ip.Results.FontName, 'FontSize', ip.Results.FontSize);
        extent = get(ht, 'extent');
        textWidth(c) = extent(3);
        delete(ht);
    end
    maxTextWidth = max(textWidth);
    offset = maxTextWidth + 2*dci;
else
    offset = 0;
end
delete(ha);

% dt = track.lifetime_s/(track.end-track.start+1);
% gapIdx = round((track.t+dt)/dt);
% gapIdx = gapIdx(track.gapVect==1) - track.start + 1;

ha = zeros(nc,nf);
set(hf, 'Position', [50, 100, width+offset, height], 'Visible', ip.Results.Visible);%, 'ResizeFcn', {@resizeCallback});

% stack index si: x + (rowi-1)*nc*nx + (c-1)*nx
for si = 1:nf
    x = rem(si-1,nx)+1;
    rowi = ceil(si/nx); % each row contains all channels
    
    for c = 1:nc
        ha(c,si) = axes('Units', 'pixels',...
            'Position', [offset+(x-1)*(wxi+dxi) height-wxi-((rowi-1)*(nc*wxi+(nc-1)*dxi+dci)+(c-1)*(wxi+dxi)) wxi wxi],...
            'XLim', [0 wxi], 'YLim', [0 wxi]);
        imagesc(xa{si}, ya{si}, trackStack{c, si}); axis image off; caxis(dRange(c,:));
        hold on;
        
        % frame index, relative to movie
        fi = si-sb+track.start-1;
        
        if ip.Results.ShowDetection && track.start<=fi && fi<=track.end% && fi>sb && fi<=(track.end-track.start+1)+sb && c==1
            idx = find(track.f==fi);
            for i = idx
                %if track.gapVect(i)==1
                    %plot(track.x(c,i), track.y(c,i), 'o', 'Color', [1 1 1], 'MarkerSize', 12);
                %else
                if ~isnan(track.isPSF(c,i)) && track.isPSF(c,i)
                    plot(track.x(c,i), track.y(c,i), 'o', 'Color', [0 1 0], 'MarkerSize', 12);
                else
                    plot(track.x(c,i), track.y(c,i), 'o', 'Color', [1 0 0], 'MarkerSize', 12);
                end
            end
        end
        
        if ip.Results.ShowMarkers && track.start<=fi && fi<=track.end
            % start marker
            if fi==track.start
                if c==1
                    plot(mean(xa{si}([1 end])), ya{si}(1), 'v', 'MarkerEdgeColor', 'none', 'Markersize', 8, 'MarkerFaceColor', [0 0 0]);
                end
                if c==nc && nc>1
                    plot(mean(xa{si}([1 end])), ya{si}(end), '^', 'MarkerEdgeColor', 'none', 'Markersize', 8, 'MarkerFaceColor', [0 0 0]);
                end
            end
            
            % gaps
            if any(track.gapVect(track.f==fi)==1) && c==1
                plot(mean(xa{si}([1 end])), ya{si}(1), 'v', 'MarkerEdgeColor', 'none', 'Markersize', 10, 'MarkerFaceColor', [0.8 0 0]);
            end
            
            % end marker
            if fi==track.end
                if c==1
                    plot(mean(xa{si}([1 end])), ya{si}(1), 'v', 'MarkerEdgeColor', 'none', 'Markersize', 8, 'MarkerFaceColor', [0 0 0]);
                end
                if c==nc && nc>1
                    plot(mean(xa{si}([1 end])), ya{si}(end), '^', 'MarkerEdgeColor', 'none', 'Markersize', 8, 'MarkerFaceColor', [0 0 0]);
                end
            end
        end
        
        % Channel labels
        if x==1 && rowi==1 && ~isempty(labels)
            ht(c) = text(-dci, wxi/2, labels{c}, 'Units', 'pixels',...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline',...
                'FontUnits', 'pixels', 'FontName', ip.Results.FontName, 'FontSize', wxi/2.5, 'Color', rgbColors{c});
        end
    end
end
colormap(gray(256));

%     function resizeCallback(src, ~)
%         [nc,nf] = size(trackStack);
%         pos = get(src, 'Position');
%         width = pos(3)-offset;
%         height = pos(4);
%
%         %[wxi, dxi, dci, nx, nr] = getProportions(width-offset, height, nf, nc);
%         [wxi, dxi, dci, nr] = getProportions(width, nx, nf, nc);
%
%         for rowi = 1:nr
%             for c = 1:nc
%                 for x = 1:nx
%                     fi = x+(rowi-1)*nx;
%                     if fi<=nf
%                         set(ha(c, fi),...
%                             'Position', [offset+(x-1)*(wxi+dxi) height-wxi-((rowi-1)*(nc*wxi+(nc-1)*dxi+dci)+(c-1)*(wxi+dxi)) wxi wxi]);
%                     end
%                 end
%                 if x==1 && rowi==1 && ~isempty(labels)
%                     set(ht(c), 'Position', [-dci, wxi/2], 'FontSize', wxi/2.5);
%                 end
%             end
%         end
%     end

if ~isempty(ip.Results.epsPath)
    print(hf, '-depsc2', ip.Results.epsPath);
end

if strcmp(ip.Results.Visible, 'off')
    close(hf);
end



end



function [wxi, dxi, dci, nr, height] = getProportions(width, nx, nf, nc)

% fixed proportions
wx = 1;
dx = 1/15; % gap between frames, relative to frame width
dc = 1/5; % gap between channels

% number of rows
nr = ceil(nf/nx);

iwidth = nx*wx + (nx-1)*dx;
iheight = nr*(nc*wx+(nc-1)*dx) + (nr-1)*dc;

height = width*iheight/iwidth;

wxi = width / ((nx/dx + nx-1)*dx);
dxi = wxi*dx;
dci = wxi*dc;

end
