% Author: Francois Aguet

function cropMovie3D(framePaths, outputDir, varargin)

nf = numel(framePaths);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('framePaths', @iscell);
ip.addRequired('outputDir', @ischar);
ip.addParamValue('FrameRange', [1 nf], @isnumeric);
ip.parse(framePaths, outputDir, varargin{:});

% max int projection of frame range
range = ip.Results.FrameRange;
proj = readtiff(framePaths{range(1)});

for i = 2:numel(range)
    proj = max(proj, readtiff(framePaths{range(i)}));
end
proj = max(proj, [], 3);

% get ROI (opens UI)
hf = figure;
imshow(proj,[]); axis image;
[ny,nx] = size(proj);
b = max(nx,ny)/10;
hr = imrect(gca, [b b nx-2*b ny-2*b]);
fcn = makeConstrainToRectFcn('imrect', get(gca,'XLim'), get(gca,'YLim'));
setPositionConstraintFcn(hr, fcn);
roi = round(wait(hr));
close(hf);
drawnow; % forces window to close

% slice vectors
xa = roi(1) + (0:roi(3)-1);
ya = roi(2) + (0:roi(4)-1);

% loop through frames, crop, save in output dir
[~,~] = mkdir(outputDir);
for i = 1:nf
    frame = readtiff(framePaths{i});
    frame = frame(ya,xa,:);
    [~,name,ext] = fileparts(framePaths{i});
    writetiff(frame, [outputDir filesep name '_crop' ext]);
end
