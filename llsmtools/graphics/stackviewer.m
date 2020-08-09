%stackviewer(stack) displays 2D projections of a 3D stack
%
% Inputs:
%      stack : 3D array
%
% Optional inputs: 
%          X : matrix of #points x 3 coordinates
%
% Parameters:
%   'ZAnisotropy' : anisotropy factor to adjust for z vs. x,y spacing differences in projections
%  'DynamicRange' : dynamic range for display. Default: [min(stack(:)) max(stack(:))]
%        'Parent' : figure handle
%
%
% Keyboard shortcuts:
% 
%              up : previous slice
%            down : next slice 
%               , : shift second channel down
%               . : shift second channel up
%               / : reset
%               c : enhance contrast
%               i : invert contrast
%               g : toggle guides
%               l : toggle legend

% Author: Francois Aguet

function varargout = stackviewer(stack, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('stack');
ip.addOptional('X', [], @(x) isstruct(x) || isnumeric(x));
ip.addParamValue('ZAnisotropy', 1);
ip.addParamValue('DynamicRange', []);
ip.addParamValue('Colormap', [], @(x) size(x,1)==nt && size(x,2)==3);
ip.addParamValue('Parent', []);
ip.addParamValue('Project', false, @islogical);
ip.addParamValue('EnhanceContrast', false, @islogical);
ip.addParamValue('ShowGuides', true, @islogical);
ip.addParamValue('ShowLegend', true, @islogical);
ip.addParamValue('SetCoord', [], @isvector);
ip.addParamValue('ScalingFactor', [], @isscalar);
ip.addParamValue('GuideColor', hsv2rgb([0.55 0.8 1]));
ip.addParamValue('GuideWidth', 1);
ip.addParamValue('GuideMode', 'full', @(x) any(strcmpi(x, {'full', 'tick'})));
ip.addParamValue('Overlay', []);
ip.parse(stack, varargin{:});


% settings
dx = 6; % spacing between axes/projections (pixels)
za = ip.Results.ZAnisotropy;

% colormap for tracks
trackColormap = ip.Results.Colormap;

dRange = ip.Results.DynamicRange;

% determine input type, convert if necessary
if ~isempty(ip.Results.Overlay) && ~iscell(stack)
    isRGB = true;
    if isempty(dRange)
        idx = isfinite(stack);
        dRange = [min(stack(idx)) max(stack(idx))];
    end
    [ny,nx,nz] = size(stack);
    stack = uint8(scaleContrast(stack, dRange));
    stackBG = stack;
    stackBG(ip.Results.Overlay~=0) = 0;
    stack(ip.Results.Overlay~=0) = 255;
    stack = {stack, stackBG, stackBG};
elseif ~iscell(stack)
    isRGB = false;
    [ny,nx,nz] = size(stack);
    if isempty(dRange)
        idx = isfinite(stack);
        dRange = [min(stack(idx)) max(stack(idx))];
        %dRange = prctile(stack(:), [0.01 99.99]);
    end
else
    isRGB = true;
    % if separate channels, convert to uint8
    if ~(size(stack{1},3)==3 && isa(stack{1}, 'uint8'))
        [ny,nx,nz] = size(stack{1});
        dRange = cellfun(@(i) [min(i(:)) max(i(:))], stack, 'unif', 0);
        for i = 1:numel(stack)
            stack{i} = uint8(scaleContrast(stack{i},dRange{i}));
        end
        if numel(stack)==2
            stack{3} = zeros(ny,nx,nz,'uint8');
        end
    else % input is cell array of RGB slices
        tmp = cat(3,stack{:});
        clear stack;
        stack{1} = tmp(:,:,1:3:end);
        stack{2} = tmp(:,:,2:3:end);
        stack{3} = tmp(:,:,3:3:end);
        clear tmp;
        [ny,nx,nz] = size(stack{1});
    end
end

% determine optional input: matrix of detection coordinates or 'tracks' structure
if ~isempty(ip.Results.X) && isstruct(ip.Results.X)
    tstruct = ip.Results.X;
    if isempty(trackColormap)
        trackColormap = hsv2rgb([rand(tstruct.n,1) ones(tstruct.n,2)]);
    end
    X = [];
else
    tstruct = [];
    X = ip.Results.X;
end

nz = nz * ip.Results.ZAnisotropy;

% setup figure properties
fpos = get(0, 'DefaultFigurePosition');
w = fpos(3); % figure dimensions in pixels
h = fpos(4);

% normalized dimensions
fx = (w-dx)/(nx+nz);
fy = (h-dx)/(ny+nz);
f = min(fx,fy);
h = (ny+nz)*f+dx; % [pixels]
w = (nx+nz)*f+dx;

% exact dimensions (for printing/scale comparisons etc.)
if ~isempty(ip.Results.ScalingFactor)
    f = ip.Results.ScalingFactor;
    w = f*(nz+nx)+dx;
    h = f*(nz+ny)+dx;    
end

hf = ip.Results.Parent;
if isempty(hf)
    hf = figure('Position', [fpos(1:2) w h], 'Color', 'w',...
        'PaperPositionMode', 'auto', 'ResizeFcn', @figResize);
else
    set(hf, 'ResizeFcn', @figResize);
end

% colormap/contrast settings
if ip.Results.EnhanceContrast
    p = 0.5;
else
    p = 1;
end
cmap = gray(256).^p;
colormap(cmap);

% frameRange = 0;
frameRange = -1:1;
% frameRange = -2:2;

% generate axes
ha(1) = axes('Position', [(nz*f+dx)/w 0 f*nx/w f*ny/h], 'Parent', hf);
ha(2) = axes('Position', [0 0 f*nz/w f*ny/h], 'Parent', hf);
ha(3) = axes('Position', [(nz*f+dx)/w (ny*f+dx)/h f*nx/w f*nz/h], 'Parent', hf);
set(ha, 'HitTest', 'off');

    function figResize(~,~)
        ipos = get(hf, 'Position');
        if strcmp(get(hf,'type'),'uipanel') && strcmp(get(hf, 'Units'), 'normalized')
            ipos = ipos.*get(gcf, 'Position');
        end
        rxy = ipos(3)/ipos(4);
        dx = 6/ipos(3);
        dy = 6/ipos(4);
        if rxy > w/h % figure is too wide
            f0 = w/h / rxy;
            set(ha(1), 'Position', [(1-f0)/2+(f0*nz*f)/w+dx 0 f0*f*nx/w f*ny/h]);
            set(ha(2), 'Position', [(1-f0)/2 0 f0*f*nz/w f*ny/h]);
            set(ha(3), 'Position', [(1-f0)/2+(f0*nz*f)/w+dx (ny*f)/h+dy f0*f*nx/w f*nz/h]);
        else
            f0 = h/w * rxy;
            set(ha(1), 'Position', [(nz*f)/w+dx 1-f0 f*nx/w f0*f*ny/h]);
            set(ha(2), 'Position', [0 1-f0 f*nz/w f0*f*ny/h]);
            set(ha(3), 'Position', [(nz*f)/w+dx 1-f0+(f0*ny*f)/h+dy f*nx/w f0*f*nz/h]);
        end
    end

hz = zoom(gcf); % needs to be figure, use gcf in case hf is a uipanel
set(hz, 'ActionPostCallback', @czoom);

% initial selection
if isempty(ip.Results.SetCoord)
    x = round(nx/2);
    y = round(ny/2);
    slice = 1;
else
    x = ip.Results.SetCoord(1);
    y = ip.Results.SetCoord(2);
    slice = ip.Results.SetCoord(3);
end
yshift = 0;

% x,y view
if ~isRGB
    if ~ip.Results.Project
        hxy = imagesc(stack(:,:,slice), 'Parent', ha(1), 'HitTest', 'off');
    else
        hxy = imagesc(max(stack,[],3), 'Parent', ha(1), 'HitTest', 'off');
    end
else
    frameRGB = cat(3, [zeros(max(yshift,0),nx); stack{1}(1+max(-yshift,0):end-max(yshift,0),:,slice); zeros(max(-yshift,0),nx)],...
        stack{2}(:,:,slice), stack{3}(:,:,slice));
    hxy = imagesc(frameRGB, 'Parent', ha(1), 'HitTest', 'off');
end
hold(ha(1), 'on');
% projection guides in x,y view
guideOpts = {'Color', ip.Results.GuideColor, 'LineWidth', ip.Results.GuideWidth, 'HitTest', 'off'};
% gw = min(nx,ny)/20;
gw = 16;
if strcmpi(ip.Results.GuideMode, 'full')
    hg(1) = plot(ha(1), [x x], [0.5 ny+0.5], guideOpts{:});
    hg(2) = plot(ha(1), [0.5 nx+0.5], [y y], guideOpts{:});
else
    hg(1) = plot(ha(1), x+zeros(1,5), [0.5 gw NaN ny-gw+0.5 ny+0.5], guideOpts{:});
    hg(2) = plot(ha(1), [0.5 gw NaN nx-gw+0.5 nx+0.5], y+zeros(1,5), guideOpts{:});
end
hl = text(0.01, 0.01, ['Slice ' num2str(slice, '%3.0d') ' (' num2str(x) ',' num2str(y) ')'], 'Color', ip.Results.GuideColor, 'Parent', ha(1),...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Units', 'normalized');

% y,z view
if ~isRGB
    if ~ip.Results.Project
        hyz = imagesc(squeeze(stack(:,x,:)), 'Parent', ha(2), 'HitTest', 'off');
    else
        hyz = imagesc(squeeze(max(stack,[],2)), 'Parent', ha(2), 'HitTest', 'off');
    end
else
    frameRGB = cat(3, squeeze(stack{1}(:,x,:)), squeeze(stack{2}(:,x,:)), squeeze(stack{3}(:,x,:)));
    hyz = imagesc(frameRGB, 'Parent', ha(2), 'HitTest', 'off');
end
hold(ha(2), 'on');
% slice line in y,z view
if strcmpi(ip.Results.GuideMode, 'full')
    hg(3) = plot(ha(2), slice*[1 1], [0.5 ny+0.5], guideOpts{:});
    hg(4) = plot(ha(2), [0.5 nz+0.5], y*[1 1], guideOpts{:});
else
    hg(3) = plot(ha(2), slice+zeros(1,5), [0.5 gw NaN ny-gw+0.5 ny+0.5], guideOpts{:});
    hg(4) = plot(ha(2), [-0.5 gw/za NaN (nz-gw+1)/za nz/za+0.5], y+zeros(1,5), guideOpts{:});
end

% x,z view
if ~isRGB
    if ~ip.Results.Project
        hxz = imagesc(squeeze(stack(y,:,:))', 'Parent', ha(3), 'HitTest', 'off');
    else
        hxz = imagesc(squeeze(max(stack,[],1))', 'Parent', ha(3), 'HitTest', 'off');
    end
else
    frameRGB = cat(3, squeeze(stack{1}(y,:,:))', squeeze(stack{2}(y,:,:))', squeeze(stack{3}(y,:,:))');
    hxz = imagesc(frameRGB, 'Parent', ha(3), 'HitTest', 'off');
end
hold(ha(3), 'on');
% slice line in x,z view
if strcmpi(ip.Results.GuideMode, 'full')
    hg(5) = plot(ha(3), [0.5 nx+0.5], slice*[1 1], guideOpts{:});
    hg(6) = plot(ha(3), x*[1 1], [0.5 nz+0.5], guideOpts{:});
else
    hg(5) = plot(ha(3), [0.5 gw NaN nx-gw+0.5 nx+0.5], slice+zeros(1,5), guideOpts{:});
    hg(6) = plot(ha(3), x+zeros(1,5), [0.5 gw/za NaN (nz-gw+1)/za nz/za+0.5], guideOpts{:});
end

set(ha, 'XTick', [], 'YTick', []);
if ~isRGB && diff(dRange)~=0
    arrayfun(@(i) caxis(i, dRange), ha, 'unif', 0);
end

% initialize plot handles in all windows
if ~isempty(X) % display detections
    opts = {'MarkerSize', 6, 'HitTest', 'off', 'LineWidth', 1};
    x0 = X(ismember(round(X(:,3)), slice+frameRange), [1 2]);
    if ~isempty(x0)
        hp(1) = plot(ha(1), x0(:,1), x0(:,2), 'go', opts{:});
    else
        hp(1) = plot(ha(1), NaN, NaN, 'go', opts{:});
    end
    
    x0 = X(ismember(round(X(:,2)), y+frameRange), [1 3]);
    if ~isempty(x0)
        hp(3) = plot(ha(3), x0(:,1), x0(:,2), 'go', opts{:});
    else
        hp(3) = plot(ha(3), NaN, NaN, 'go', opts{:});
    end
    
    x0 = X(ismember(round(X(:,1)), x+frameRange), [2 3]);
    if ~isempty(x0)
        hp(2) = plot(ha(2), x0(:,2), x0(:,1), 'go', opts{:});
    else
        hp(2) = plot(ha(2), NaN, NaN, 'go', opts{:});
    end
end
hpx = [];

axis(ha(1), 'equal');
% this fixes a bug with axis 'equal' that allows panning beyond boundaries
set(ha(1), 'XLim', [0.5 nx+0.5], 'YLim', [0.5 ny+0.5]);


% Set callbacks
set(ha, 'ButtonDownFcn', @click_Callback); % after 'hold on'
set(gcf, 'WindowScrollWheelFcn', @scroll_Callback);
set(gcf, 'KeyPressFcn', @key_Callback);

hpan = pan;
set(hpan,'ActionPreCallback',@panstart);
set(hpan,'ActionPostCallback',@panstop);

% settings
showGuides = ip.Results.ShowGuides;
if ~showGuides
    set(hg, 'Visible', 'off');
end
showLegend = ip.Results.ShowLegend;
if ~showLegend
    set(hl, 'Visible', 'off');
end
invertContrast = false;

if nargout>0
    varargout{1} = ha;
end

    % in this callback, set click-and-drag callbacks
    function click_Callback(varargin)
        switch gca
            case ha(1)
                a = get(gca, 'CurrentPoint');
                x = round(a(1,1));
                y = round(a(1,2));
                updateProj(); % required for clicking w/o dragging
                set(gcf, 'WindowButtonMotionFcn', @dragProj, 'WindowButtonUpFcn', @stopDragging);
            case ha(2)
                a = get(gca,'CurrentPoint');
                slice = round(a(1,1));
                updateSlice();
                set(gcf, 'WindowButtonMotionFcn', @dragSlice, 'WindowButtonUpFcn', @stopDragging);
            case ha(3)
                a = get(gca,'CurrentPoint');
                slice = round(a(1,2));
                updateSlice();
                set(gcf, 'WindowButtonMotionFcn', @dragSlice, 'WindowButtonUpFcn', @stopDragging);
        end
    end

    function dragProj(varargin)
        a = get(gca, 'CurrentPoint');
        x = min(max(1, round(a(1,1))), nx);
        y = min(max(1, round(a(1,2))), ny);
        updateProj();
    end

    function dragSlice(varargin)
        a = get(gca, 'CurrentPoint');
        switch gca
            case ha(2)
                slice = min(max(1,round(a(1,1))),round(nz/ip.Results.ZAnisotropy));
            case ha(3)
                slice = min(max(1,round(a(1,2))),round(nz/ip.Results.ZAnisotropy));
        end
        updateSlice();
    end

    function stopDragging(varargin)
        set(gcf, 'WindowButtonMotionFcn', '');
    end

    % scroll through stack slices
    function scroll_Callback(~, eventdata)
        if eventdata.VerticalScrollCount < 0
            if slice < floor(nz/ip.Results.ZAnisotropy)
                slice = slice + 1;
                updateSlice();
            end
        elseif eventdata.VerticalScrollCount > 0
            if slice > 1
                slice = slice - 1;
                updateSlice();
            end
        end
    end

    function key_Callback(~, eventdata)
        switch eventdata.Key
            case 'uparrow'
                if slice > 1
                    slice = slice - 1;
                    updateSlice();
                end
            case 'downarrow'
                if slice < nz/ip.Results.ZAnisotropy
                    slice = slice + 1;
                    updateSlice();
                end
            case 'comma'
                yshift = yshift-1;
                updateSlice();
            case 'period'
                yshift = yshift+1;
                updateSlice();
            case 'slash'
                yshift = 0;
                updateSlice;
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
            case 'g'
                showGuides = ~showGuides;
                if showGuides
                    set(hg, 'Visible', 'on');
                else
                    set(hg, 'Visible', 'off');
                end
            case 'l'
                showLegend = ~showLegend;
                if showLegend
                    set(hl, 'Visible', 'on');
                else
                    set(hl, 'Visible', 'off');
                end
        end
    end

    function setColormap()
        if ~invertContrast
            cmap = gray(256).^p;
        else
            cmap = gray(256).^(1/p);
            cmap = cmap(end:-1:1,:,:);
        end
        colormap(cmap);
    end

    function updateSlice()
        if ~isRGB
            set(hxy, 'CData', stack(:,:,slice));
        else
            frameRGB = cat(3, [zeros(max(yshift,0),nx); stack{1}(1+max(-yshift,0):end-max(yshift,0),:,slice); zeros(max(-yshift,0),nx)],...
                stack{2}(:,:,slice), stack{3}(:,:,slice));
            set(hxy, 'CData', frameRGB);
        end
        % x,y cross-sections
        set(hg(3), 'XData', slice+zeros(size(get(hg(3), 'XData'))));
        set(hg(5), 'YData', slice+zeros(size(get(hg(5), 'YData'))));
        
        set(hl, 'String', ['Slice ' num2str(slice) ' (' num2str(x) ',' num2str(y) ')']);
        
        if ~isempty(tstruct)
            vidx = ~isnan(tstruct.X(slice,:));
            if slice~=1
                delete(hpx);
                set(ha(1), 'ColorOrder', trackColormap(tstruct.idx(vidx),:));
                hpx = plot(ha(1), tstruct.X(1:slice,vidx), tstruct.Y(1:slice,vidx), 'HitTest', 'off');
            else
                delete(hpx);
                hpx = [];
            end
        end

        % points in current slice
        if ~isempty(X)
            x0 = X(ismember(round(X(:,3)), slice+frameRange), [1 2]);
            set(hp(1), 'XData', x0(:,1), 'YData', x0(:,2));
        end
    end

    function updateProj()
        % plot lines
        set(hg(1), 'XData', x+zeros(size(get(hg(1), 'XData'))));
        set(hg(2), 'YData', y+zeros(size(get(hg(1), 'XData'))));
        
        set(hg(4), 'YData', y+zeros(size(get(hg(4), 'YData'))));
        set(hg(6), 'XData', x+zeros(size(get(hg(6), 'XData'))));
        
        % update data
        xi = min(max(x,1),nx);
        yi = min(max(y,1),ny);
        if ~isRGB
            set(hyz, 'CData', squeeze(stack(:,xi,:)));
            set(hxz, 'CData', squeeze(stack(yi,:,:))');
        else
            frameRGB = cat(3, squeeze(stack{1}(:,xi,:)), squeeze(stack{2}(:,xi,:)),...
                squeeze(stack{3}(:,xi,:)));
            set(hyz, 'CData', frameRGB);
            frameRGB = cat(3, squeeze(stack{1}(yi,:,:))',...
                squeeze(stack{2}(yi,:,:))', squeeze(stack{3}(yi,:,:))');
            set(hxz, 'CData', frameRGB);
        end
        
        % points in current slice
        if ~isempty(X)
            % x,z view
            x0 = X(ismember(round(X(:,2)), yi+frameRange), [1 3]);
            set(hp(3), 'XData', x0(:,1), 'YData', x0(:,2));
            
            x0 = X(ismember(round(X(:,1)), xi+frameRange), [2 3]);
            set(hp(2), 'XData', x0(:,2), 'YData', x0(:,1));
        end
        
        set(hl, 'String', ['Slice ' num2str(slice) ' (' num2str(x) ',' num2str(y) ')']);
    end

    function czoom(~, eventdata)
        switch eventdata.Axes
            case ha(1)
                XLim = get(ha(1), 'XLim');
                YLim = get(ha(1), 'YLim');
                set(ha(2), 'YLim', YLim);
                set(ha(3), 'XLim', XLim);
        end
    end

    % Pan functions
    function panstart(~, eventdata)
        set(gcf, 'WindowButtonMotionFcn', {@dopan, eventdata});
    end

    function panstop(varargin)
        set(gcf, 'WindowButtonMotionFcn', '');
    end

    function dopan(~,~,eventdata)
        switch eventdata.Axes
            case ha(1)
                XLim = get(ha(1), 'XLim');
                YLim = get(ha(1), 'YLim');
                set(ha(2), 'YLim', YLim);
                set(ha(3), 'XLim', XLim);
            case ha(2)
                YLim = get(ha(2), 'YLim');
                set(ha(1), 'YLim', YLim);
            case ha(3)
                XLim = get(ha(3), 'XLim');
                set(ha(1), 'XLim', XLim);
        end
    end

end
