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

ip.parse(dataFile, xyPixelSize, dz, varargin{:});

pr = ip.Results;

Rotate = pr.Rotate;

SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
nphases = pr.nphases;

if iscell(dataFile)
    dataFile = char(dataFile);
end

[filepath,name,ext] = fileparts(dataFile);
fol = [filepath filesep];

if(Rotate)
    if ~exist([fol 'DSR'],'dir')
        mkdir([fol 'DSR']);
        fileattrib([fol 'DSR'], '+w', 'g');
    end
else
    if ~exist([fol 'DS'],'dir')
        mkdir([fol 'DS']);
        fileattrib([fol 'DS'], '+w', 'g');
    end
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
writetiff(dsim, [fol 'DS' filesep name ext])
toc


end

