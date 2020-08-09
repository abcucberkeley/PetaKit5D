

function X = rotateDetections3D(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data');
ip.addOptional('X', []);
ip.addParamValue('Crop', true, @islogical);
ip.parse(data, varargin{:});

X = ip.Results.X;

vol = double(readtiff(data.framePathsDS{1}{1}));
[ny,nx,nz] = size(vol);

theta = data.angle*pi/180;

doCrop = false;
if ~data.objectiveScan
    % calculate height of rotated volume
    % project and find left/right edges of deskewed rhomboid, take median width
    proj = squeeze(max(vol,[],1))'~=0;
    proj = diff(proj,1,2);
    
    startIndex = ones(nz,1);
    endIndex = nx*ones(nz,1);
    [srow,scol] = find(proj==1);
    [erow,ecol] = find(proj==-1);
    startIndex(srow) = scol;
    endIndex(erow) = ecol;
    
    doCrop = any(startIndex>1 & endIndex<nx);
    
    % calculate height; first & last 2 frames have interpolation artifacts
    w = median(diff([startIndex endIndex],[],2));
    h = round(abs(w*tan(theta)*cos(theta)))-4;
end


center = ([ny nx nz]+1)/2;
T1 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    -center([2 1 3]) 1];

S = [1 0 0 0
    0 1 0 0
    0 0 data.zAniso 0
    0 0 0 1];

% Rotate x,z
R = [cos(theta) 0 -sin(theta) 0; % order for imwarp is x,y,z
    0 1 0 0;
    sin(theta) 0 cos(theta) 0;
    0 0 0 1];

if ip.Results.Crop && doCrop
    outSize = round([ny nx/cos(theta) h]);
else
    % exact proportions of rotated box
    outSize = round([ny nx*cos(theta)+nz*ip.Results.zxRatio*sin(abs(theta)) nz*ip.Results.zxRatio*cos(theta)+nx*sin(abs(theta))]);
end

T2 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    (outSize([2 1 3])+1)/2 1];

tform = affine3d(T1*S*R*T2);

[x,y,z] = transformPointsForward(tform, X(:,1), X(:,2), X(:,3));
X = [x(:) y(:) z(:)];
