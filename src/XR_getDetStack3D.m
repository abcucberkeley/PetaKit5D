function [stack, xa, ya, za] = XR_getDetStack3D(data, det, t, varargin)
% get the cropped stack for each xyz detected coordinate @ specified t
% 
% copied from Gokul's code

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('det', @isstruct);
ip.addRequired('t', @isnumeric);
ip.addParameter('WindowWidth', [8 8], @isnumeric);
ip.addParameter('Mode', 'max', @(x) any(strcmpi(x, {'max', 'center'})));
ip.addParameter('Filepath', 'framePathsDS', @(x) any(strcmpi(x, {'framePaths', 'framePathsDS', 'framePathsDSR'})));
ip.addParameter('figureSize', [3 3], @isnumeric);
ip.parse(data, det, t, varargin{:});
% # of t need to match # of objects
% the res output assumes a background of 112 (sCMOS DC + 3*SD)

if strcmpi(ip.Results.Filepath, 'framePathsDS')
    filePath =data.framePathsDS;
elseif strcmpi(ip.Results.Filepath, 'framePaths')
    filePath =data.framePaths;
elseif strcmpi(ip.Results.Filepath, 'framePathsDSR')
    filePath =data.framePathsDSR;
end
nc = length(data.channels);

w = ceil(ip.Results.WindowWidth);

[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;

% coordinate matrices
xv = det.x;
yv = det.y;
zv = det.z;
ii = det.IntegratedIntensity;
A = det.A;

info = imfinfo(filePath{1}{1});
nz = numel(info);
nx = info(1).Width;
ny = info(1).Height;

% start and end buffer sizes
sb = 0;
eb = 0;

% frame index

nf = length(t);

% reverse z index
% zv = data.imagesize(3)+1-zv;


xi = round(xv(1,:));
yi = round(yv(1,:));
zi = round(zv(1,:));
% ensure that window falls within frame bounds
x0 = xi - min([xi-1 w(1)]);
x1 = xi + min([nx-xi w(1)]);
y0 = yi - min([yi-1 w(1)]);
y1 = yi + min([ny-yi w(1)]);
z0 = zi - min([zi-1 w(2)]);
z1 = zi + min([nz-zi w(2)]);
% axes for each frame
xa = arrayfun(@(i) x0(i):x1(i), 1:nf, 'unif', 0);
ya = arrayfun(@(i) y0(i):y1(i), 1:nf, 'unif', 0);
za = arrayfun(@(i) z0(i):z1(i), 1:nf, 'unif', 0);


% load all visible frames of this track and store
stack = cell(nc,nf);
for c = 1:nc
    for k = 1:numel(t)
        frame = readtiff(filePath{c}{t});
        stack{c,k} = frame(ya{k}, xa{k}, za{k});
        
        % Also not generate image by rounding the mean position, but use
        % subpixel Gaussian patterns
        % generate a model gaussian
        [sy, sx, sz] = size(stack{c,k});
        gim = zeros(sy, sx, sz);
        % gim(round(sy/2),round(sx/2),round(sz/2)) = ii;
        % gim = filterGauss3D(gim, [det.s(1), det.s(2)]);
        pts_loc = [xv(1,:), yv(1, :), zv(1, :)]; 
        data_range = [x0, x1; y0, y1; z0, z1];
        sigma = det.s;
        I_gauss = generate_gaussian_pattern_for_subpixel_location(pts_loc, sigma, data_range);
        gim = I_gauss * A;
        
        % imModel{c,k} = gim + det.c;
        imModel{c,k} = gim + det.c;
        
        % difference
        imModelDiff{c,k} = double(stack{c,k})-imModel{c,k};
        
%         % residual
%         exp = stack{c,k};
%         exp = exp - 108;
%         exp(exp<0) = 0;
%         
%         resExp = imModelDiff{c,k};
%         resExp = resExp - 108;
%         resExp(resExp<0) = 0;
%         res{c,k} = sum(resExp(:))/sum(exp(:))*100;
    end
end


% plot the image
if strcmpi(ip.Results.Mode, 'center')
    xzsec = cellfun(@(i) squeeze(i(w(1)+1,:,:)), stack, 'unif', 0);
    yzsec = cellfun(@(i) squeeze(i(:,w(1)+1,:))', stack, 'unif', 0);
    xysec = cellfun(@(i) i(:,:,w(2)+1), stack, 'unif', 0);
    
    gxzsec = cellfun(@(i) squeeze(i(w(1)+1,:,:)), imModel, 'unif', 0);
    gyzsec = cellfun(@(i) squeeze(i(:,w(1)+1,:))', imModel, 'unif', 0);
    gxysec = cellfun(@(i) i(:,:,w(2)+1), imModel, 'unif', 0);
    
    gdxzsec = cellfun(@(i) squeeze(i(w(1)+1,:,:)), imModelDiff, 'unif', 0);
    gdyzsec = cellfun(@(i) squeeze(i(:,w(1)+1,:))', imModelDiff, 'unif', 0);
    gdxysec = cellfun(@(i) i(:,:,w(2)+1), imModelDiff, 'unif', 0);
    
else
    % max proj
    xzsec = cellfun(@(i) squeeze(max(i,[],1)), stack, 'unif', 0);
    yzsec = cellfun(@(i) squeeze(max(i,[],2))', stack, 'unif', 0);
    xysec = cellfun(@(i) max(i,[],3), stack, 'unif', 0);
    
    gxzsec = cellfun(@(i) squeeze(max(i,[],1)), imModel, 'unif', 0);
    gyzsec = cellfun(@(i) squeeze(max(i,[],2))', imModel, 'unif', 0);
    gxysec = cellfun(@(i) max(i,[],3), imModel, 'unif', 0);
    
    gdxzsec = cellfun(@(i) squeeze(max(i,[],1)), imModelDiff, 'unif', 0);
    gdyzsec = cellfun(@(i) squeeze(max(i,[],2))', imModelDiff, 'unif', 0);
    gdxysec = cellfun(@(i) max(i,[],3), imModelDiff, 'unif', 0);
end


% tmp = cat(1,xysec{:});
% dRange = [min(tmp(:)) max(tmp(:))];
% dRange = prctile(double(tmp(:)), [0.1 99.9]);
dRange = prctile(double(xysec{:}(:)), [0.1 99.9]);

gdRange = prctile(double(gxysec{:}(:)), [0.1 99.9]);

gddRange = prctile(double(gdxysec{:}(:)), [0.1 99.9]);

figureSize = ip.Results.figureSize;
ha = setupFigure(3, 3, 'SameAxes', false, 'AxesWidth', figureSize(1) / 3, 'AxesHeight', figureSize(2) / 3,...
    'XSpace', [1 0.2 0.1], 'YSpace', [0.2 0.1 0.1]);
i = 1;
ha = reshape(ha, [3, 3]);
cmap = hot(256);
colormap(cmap);

set(ha(:), 'XTick', [], 'YTick', []);
imagesc(xa{i}, ya{i}, xysec{i}, 'Parent', ha(1,1)); caxis(ha(1,1), dRange);
hold on, plot(xi, yi, 'o', 'MarkerSize', 15, 'MarkerEdgeColor', ceG, 'Linewidth', 1.5, 'Parent', ha(1,1)), hold off;
imagesc(xa{i}, za{i}, xzsec{i}', 'Parent', ha(1,2)); caxis(ha(1,2), dRange);
hold on, plot(xi, zi, 'o', 'MarkerSize', 15, 'MarkerEdgeColor', ceG, 'Linewidth', 1.5, 'Parent', ha(1,2)), hold off;
imagesc(ya{i}, za{i}, yzsec{i}, 'Parent', ha(1,3)); caxis(ha(1,3), dRange);
hold on, plot(yi, zi, 'o', 'MarkerSize', 15, 'MarkerEdgeColor', ceG, 'Linewidth', 1.5, 'Parent', ha(1,3)), hold off;

imagesc(xa{i}, ya{i}, gxysec{i}, 'Parent', ha(2,1)); caxis(ha(2,1), dRange);
imagesc(xa{i}, za{i}, gxzsec{i}', 'Parent', ha(2,2)); caxis(ha(2,2), dRange);
imagesc(ya{i}, za{i}, gyzsec{i}, 'Parent', ha(2,3)); caxis(ha(2,3), dRange);

imagesc(xa{i}, ya{i}, gdxysec{i}, 'Parent', ha(3,1)); caxis(ha(3,1), gddRange);
hold on, plot(xi, yi, 'o', 'MarkerSize', 15, 'MarkerEdgeColor', ceG, 'Linewidth', 1.5, 'Parent', ha(3,1)), hold off;
imagesc(xa{i}, za{i}, gdxzsec{i}', 'Parent', ha(3,2)); caxis(ha(3,2), gddRange);
hold on, plot(xi, zi, 'o', 'MarkerSize', 15, 'MarkerEdgeColor', ceG, 'Linewidth', 1.5, 'Parent', ha(3,2)), hold off;
imagesc(ya{i}, za{i}, gdyzsec{i}, 'Parent', ha(3,3)); caxis(ha(3,3), gddRange); set(ha(3, 3), 'yDir', 'normal');
hold on, plot(yi, zi, 'o', 'MarkerSize', 15, 'MarkerEdgeColor', ceG, 'Linewidth', 1.5, 'Parent', ha(3,3)), hold off;
axis(ha(:), 'image');
axis(ha(:), 'off');

