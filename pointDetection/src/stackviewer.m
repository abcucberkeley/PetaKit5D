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

% Francois Aguet (last update: 08/25/2013)

function varargout = stackviewer(stack, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('stack');
ip.addOptional('X', [], @(x) isstruct(x) || isnumeric(x));
ip.addParamValue('ZAnisotropy', 2.5);
ip.addParamValue('DynamicRange', []);
ip.addParamValue('Colormap', [], @(x) size(x,1)==nt && size(x,2)==3);
ip.addParamValue('Parent', []);
ip.addParamValue('Project', false, @islogical);
ip.addParamValue('EnhanceContrast', false, @islogical);
ip.addParamValue('ShowGuides', true, @islogical);
ip.addParamValue('SetCoord', [], @isvector);
ip.parse(stack, varargin{:});

% colormap for tracks
trackColormap = ip.Results.Colormap;

if ~iscell(stack)
    nc = 1;
    [ny,nx,nz] = size(stack);
else
    nc = numel(stack);
    [ny,nx,nz] = size(stack{1});
end
    
% dynamic range of the stack
dRange = ip.Results.DynamicRange;
if isempty(dRange)
    if nc==1
        dRange = [min(stack(:)) max(stack(:))];
    else
        dRange = cellfun(@(i) [min(i(:)) max(i(:))], stack, 'unif', 0);
    end
    %dRange = prctile(stack(:), [0.01 99.99]);
end

% select optional input: matrix of detection coordinates or 'tracks' structure
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

% spacing between axes/projections
dx = 6; % pixels

% normalized dimensions
fx = (w-dx)/(nx+nz);
fy = (h-dx)/(ny+nz);
f = min(fx,fy);
h = (ny+nz)*f+dx; % [pixels]
w = (nx+nz)*f+dx;

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
invertContrast = false;
cmap = gray(256).^p;
colormap(cmap);

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
        if ~isempty(ipos)
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

showLegend = true;

% x,y view
if nc==1
    if ~ip.Results.Project
        hxy = imagesc(stack(:,:,slice), 'Parent', ha(1), 'HitTest', 'off');
    else
        hxy = imagesc(max(stack,[],3), 'Parent', ha(1), 'HitTest', 'off');
    end
else
    frameRGB = cat(3, uint8(scaleContrast([zeros(max(yshift,0),nx); stack{1}(1+max(-yshift,0):end-max(yshift,0),:,slice); zeros(max(-yshift,0),nx)],dRange{1})),...
        uint8(scaleContrast(stack{2}(:,:,slice),dRange{2})), zeros(ny,nx,'uint8'));
    hxy = imagesc(frameRGB, 'Parent', ha(1), 'HitTest', 'off');
end
hold(ha(1), 'on');
% projection lines in x-y view
hl(1) = plot(ha(1), [x x], [0.5 ny+0.5], 'r', 'HitTest', 'off');
hl(2) = plot(ha(1), [0.5 nx+0.5], [y y], 'r', 'HitTest', 'off');
if showLegend
    ht = text(0.01, 0.01, ['Slice ' num2str(slice)], 'Color', 'r', 'Parent', ha(1),...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Units', 'normalized');
end

% y,z view
if nc==1
    if ~ip.Results.Project
        hyz = imagesc(squeeze(stack(:,x,:)), 'Parent', ha(2), 'HitTest', 'off');
    else
        hyz = imagesc(squeeze(max(stack,[],2)), 'Parent', ha(2), 'HitTest', 'off');
    end
else
    frameRGB = cat(3, uint8(scaleContrast(squeeze(stack{1}(:,x,:)),dRange{1})),...
        uint8(scaleContrast(squeeze(stack{2}(:,x,:)),dRange{2})), zeros(ny,size(stack{1},3),'uint8'));
    hyz = imagesc(frameRGB, 'Parent', ha(2), 'HitTest', 'off');
end
hold(ha(2), 'on');
% slice line in y,z view
hl(3) = plot(ha(2), slice*[1 1], [0.5 ny+0.5], 'r', 'HitTest', 'off');

% x,z view
if nc==1
    if ~ip.Results.Project
        hxz = imagesc(squeeze(stack(y,:,:))', 'Parent', ha(3), 'HitTest', 'off');
    else
        hxz = imagesc(squeeze(max(stack,[],1))', 'Parent', ha(3), 'HitTest', 'off');
    end
else
    frameRGB = cat(3, uint8(scaleContrast(squeeze(stack{1}(y,:,:))',dRange{1})),...
        uint8(scaleContrast(squeeze(stack{2}(y,:,:))',dRange{2})), zeros(size(stack{1},3),nx,'uint8'));
    hxz = imagesc(frameRGB, 'Parent', ha(3), 'HitTest', 'off');
end
hold(ha(3), 'on');
% slice line in x,z view
hl(4) = plot(ha(3), [0.5 nx+0.5], slice*[1 1], 'r', 'HitTest', 'off');

set(ha, 'XTick', [], 'YTick', []);
if nc==1 && diff(dRange)~=0
    arrayfun(@(i) caxis(i, dRange), ha, 'unif', 0);
end

% initialize plot handles in all windows
if ~isempty(X) % display detections
    hp(1) = plot(ha(1), NaN, NaN, 'go', 'MarkerSize', 10);
    hp(2) = plot(ha(2), NaN, NaN, 'go', 'MarkerSize', 10);
    hp(3) = plot(ha(3), NaN, NaN, 'go', 'MarkerSize', 10);
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

guides = ip.Results.ShowGuides;
if ~guides
    set([hl ht], 'Visible', 'off');
end

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
        x = round(a(1,1));
        y = round(a(1,2));
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
            if slice < nz/ip.Results.ZAnisotropy
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
                guides = ~guides;
                if guides
                    set([hl ht], 'Visible', 'on');
                else
                    set([hl ht], 'Visible', 'off');
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
        if nc==1
            set(hxy, 'CData', stack(:,:,slice));
        else
            frameRGB = cat(3, uint8(scaleContrast([zeros(max(yshift,0),nx); stack{1}(1+max(-yshift,0):end-max(yshift,0),:,slice); zeros(max(-yshift,0),nx)],dRange{1})),...
                uint8(scaleContrast(stack{2}(:,:,slice),dRange{2})), zeros(ny,nx,'uint8'));
            set(hxy, 'CData', frameRGB);
        end
        % x,y cross-sections
        set(hl(3), 'XData', slice*[1 1]);
        set(hl(4), 'YData', slice*[1 1]);
        set(ht, 'String', ['Slice ' num2str(slice)]);
        
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
            x0 = X(ismember(round(X(:,3)), slice+(-1:1)), [1 2]);
            set(hp(1), 'XData', x0(:,1), 'YData', x0(:,2));
        end
    end

    function updateProj()
        % plot lines
        set(hl(1), 'XData', x*[1 1]);
        set(hl(2), 'YData', y*[1 1]);
        
        % update data
        xi = min(max(x,1),nx);
        yi = min(max(y,1),ny);
        if nc==1
            set(hyz, 'CData', squeeze(stack(:,xi,:)));
            set(hxz, 'CData', squeeze(stack(yi,:,:))');
        else
            frameRGB = cat(3, uint8(scaleContrast(squeeze(stack{1}(:,x,:)),dRange{1})),...
                uint8(scaleContrast(squeeze(stack{2}(:,x,:)),dRange{2})), zeros(ny,size(stack{1},3),'uint8'));
            set(hyz, 'CData', frameRGB);
            frameRGB = cat(3, uint8(scaleContrast(squeeze(stack{1}(y,:,:))',dRange{1})),...
                uint8(scaleContrast(squeeze(stack{2}(y,:,:))',dRange{2})), zeros(size(stack{1},3),nx,'uint8'));
            set(hxz, 'CData', frameRGB);
        end
        
        
        % points in current slice
        if ~isempty(X)
            % x,z view
            x0 = X(ismember(round(X(:,2)), yi+(-1:1)), [1 3]);
            set(hp(3), 'XData', x0(:,1), 'YData', x0(:,2));
            
            x0 = X(ismember(round(X(:,1)), xi+(-1:1)), [2 3]);
            set(hp(2), 'XData', x0(:,2), 'YData', x0(:,1));
        end
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
