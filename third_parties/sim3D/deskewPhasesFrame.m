function [] = deskewPhasesFrame(dataFile,xyPixelSize,dz,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataFile');
ip.addRequired('xyPixelSize'); % typical value: 0.1
ip.addRequired('dz'); % typical value: 0.2-0.5

ip.addParameter('Rotate',false,@islogical); % Rotate after deskew

ip.addOptional('SkewAngle', 32.45, @isscalar);
ip.addOptional('Reverse', true, @islogical);
ip.addParameter('nphases', 5, @isnumeric);
ip.addParameter('Save16bit', false, @islogical);

ip.parse(dataFile, xyPixelSize, dz, varargin{:});

pr = ip.Results;

Rotate = pr.Rotate;

SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
nphases = pr.nphases;
Save16bit = pr.Save16bit;

if iscell(dataFile)
    dataFile = char(dataFile);
end

[filepath,name,ext] = fileparts(dataFile);
fol = [filepath '/'];

if(Rotate)
    folStr = 'DSR';
else
    folStr = 'DS';
end

if ~exist([fol folStr],'dir')
    mkdir([fol folStr]);
    fileattrib([fol folStr], '+w', 'g');
end


tic
im = readtiff(dataFile);
[~,~,sz] = size(im);
for p = 1:nphases
    vol = im(:,:,p:nphases:end);
    volout = deskewFrame3D(vol, SkewAngle, dz, xyPixelSize, 'reverse', Reverse);
    [sy, sx, ~] = size(volout);
    if p == 1
        dsim = zeros([sy,sx,sz], 'single');
    end
    dsim(:,:,p:nphases:end) = volout;
end

if(Rotate)
    theta = SkewAngle * pi / 180;
    dz0 = sin(theta) * dz; % ~0.25 for dz0 = 0.45
    zAniso = dz0 / xyPixelSize;
    dsim = rotateFrame3D(dsim, SkewAngle, zAniso, Reverse, 'Crop', true, 'ObjectiveScan', true);
end

if Save16bit
    dsim = uint16(dsim);
else
    dsim = single(dsim);
end

writetiff(dsim, [fol folStr '/' name ext]);

toc


end

