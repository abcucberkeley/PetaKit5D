function XR_unmix_channels_frame(frameFullpaths, unmixFactors, varargin)
% Unmixing a channel by subtracting another channel with a factor
%
% Author: Xiongtao Ruan (09/17/2025)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('unmixFactors', @(x) isnumeric(x));
ip.addParameter('mode', 'linear', @ischar); % linear vs gaussian
ip.addParameter('unmixSigmas', [], @isnumeric); 
ip.addParameter('resultDirName', 'Unmixed', @ischar); 
ip.addParameter('channelInd', 1, @isnumeric); % unmix for which channel
% flat field parameters
ip.addParameter('FFCorrection', false, @islogical);
ip.addParameter('lowerLimit', 0.4, @isnumeric);
ip.addParameter('FFImagePaths', {'','',''}, @iscell);
ip.addParameter('backgroundPaths', {'','',''}, @iscell);
% background subtraction with constant background instead of background image
ip.addParameter('constBackground', [], @isnumeric);
% constant offset after unmixing
ip.addParameter('constOffset', [], @(x) isnumeric(x));
% input and output parameters
ip.addParameter('save16bit', true, @islogical); 
ip.addParameter('zarrFile', false, @islogical); 
ip.addParameter('saveZarr', false, @islogical); 
ip.addParameter('blockSize', [256, 256, 256] , @isvector); % in y, x, z
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(frameFullpaths, unmixFactors, varargin{:});

pr = ip.Results;
mode = pr.mode;
unmixSigmas = pr.unmixSigmas;
resultDirName = pr.resultDirName;
channelInd = pr.channelInd;
FFCorrection = pr.FFCorrection;
lowerLimit = pr.lowerLimit;
FFImagePaths = pr.FFImagePaths;
backgroundPaths = pr.backgroundPaths;
constBackground = pr.constBackground;
constOffset = pr.constOffset;
save16bit = pr.save16bit;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;

frameFullpaths = strip(frameFullpaths, 'right', filesep);
[dataPath, fsname, ext] = fileparts(frameFullpaths{channelInd});

unmixPath = [dataPath, '/', resultDirName, '/'];
if ~exist(unmixPath, 'dir')
    mkdir(unmixPath);
end

if saveZarr
    unmixFullpath = sprintf('%s/%s.zarr', unmixPath, fsname);
    unmixTmppath = sprintf('%s/%s_%s.zarr', unmixPath, fsname, uuid);    
else
    unmixFullpath = sprintf('%s/%s.tif', unmixPath, fsname);
    unmixTmppath = sprintf('%s/%s_%s.tif', unmixPath, fsname, uuid);
end

done_flag = exist(unmixFullpath, 'file');
if all(done_flag)
    disp('The output results exist, skip it!');
    return;
end

fprintf('Start channel unmixing for %s...\n', fsname);

tic

imSize = getImageSize(frameFullpaths{channelInd});

nc = numel(frameFullpaths);

im = zeros(imSize, 'single');
for c = 1 : nc
    if zarrFile
        im_c = readzarr(frameFullpaths{c});
    else
        im_c = readtiff(frameFullpaths{c});
    end
    im_c = single(im_c);

    % flat field correction
    if FFCorrection
        FFImage = FFImagePaths{c};
        backgroundImage = backgroundPaths{c};
        fprintf(['Flat-field correction for frame %s...\n', ...
            '  Flat-field image: %s\n  Background image: %s\n'], ...
            frameFullpaths{c}, FFImage, backgroundImage);
        LSIm = readtiff(FFImage);
        BKIm = readtiff(backgroundImage);
        % do not add background offset after ff correction.
        im_c = XR_LSFlatFieldCorrection(im_c,LSIm,BKIm,'lowerLimit', lowerLimit, ...
            'constOffset', 0);
    end
    if ~isempty(constBackground)
        im_c = im_c - constBackground(c);
    end

    switch mode
        case 'linear'
        case 'gaussian'
            if ~all(unmixSigmas(c, :) == 0)
                im_c = imgaussfilt3(im_c, unmixSigmas(c, :));
            end
    end
    im = im + im_c * unmixFactors(c);
end

im = max(0, im);

% add background offset after unmixing if it is nozero.
if ~isempty(constOffset) && constOffset ~= 0
    im = im + constOffset;
end

if save16bit
    im = cast(im, 'uint16');
end

% write the result
if saveZarr
    writezarr(im, unmixTmppath, blockSize=blockSize, create=true)
else
    writetiff(im, unmixTmppath)
end

if exist(unmixTmppath, 'file') || exist(unmixTmppath, 'dir')
    movefile(unmixTmppath, unmixFullpath);
end

end
