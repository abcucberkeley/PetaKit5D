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

% the channels that contribute, and those flat-field corrected or filtered
% (in single) before the sum
used = false(1, nc);
filtered = false(1, nc);
for c = 1 : nc
    if unmixFactors(c) == 0
        continue;
    end
    used(c) = true;
    filtered(c) = strcmp(mode, 'gaussian') && ~all(unmixSigmas(c, :) == 0);
end
processed = used & (FFCorrection | filtered);

dtype = 'single';
if save16bit
    dtype = 'uint16';
end

% the mex sums the channels as read in one pass when they and the result
% take at most 8 bytes per voxel (the original code held at least a single
% sum and a single channel); otherwise the channels are added one at a time
% to the sum, which holds no more than the original code
atOnce = any(used) && ~any(processed) && sum(used) + dataTypeToByteNumber(dtype) <= 8 ...
    && exist('unmix_channels_sum_mex', 'file') == 3;
if atOnce
    dtypes = cellfun(@(x) frame_data_type(x, zarrFile), frameFullpaths(used), 'unif', 0);
    atOnce = all(ismember(dtypes, {'uint8', 'uint16', 'single', 'double'})) ...
        && sum(cellfun(@dataTypeToByteNumber, dtypes)) + dataTypeToByteNumber(dtype) <= 8;
end

ims = cell(1, nc);
bgs = zeros(1, nc);
im = [];
for c = find(used)
    if zarrFile
        im_c = readzarr(frameFullpaths{c});
    else
        im_c = readtiff(frameFullpaths{c});
    end

    % flat field correction, background and gaussian filtering in single
    % (the background of single channels in place); the background of the
    % other channels is subtracted in the sum
    bg_c = [];
    if processed(c) || isa(im_c, 'single')
        im_c = single(im_c);
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
        if filtered(c)
            im_c = imgaussfilt3(im_c, unmixSigmas(c, :));
        end
    elseif ~isempty(constBackground)
        bg_c = constBackground(c);
    end
    if atOnce
        ims{c} = im_c;
        if ~isempty(bg_c)
            bgs(c) = bg_c;
        end
    else
        im = unmix_channels_sum({im_c}, unmixFactors(c), bg_c, [], '', im);
    end
end
clear im_c;

if atOnce
    im = unmix_channels_sum(ims(used), unmixFactors(used), bgs(used), constOffset, dtype);
    clear ims;
else
    if isempty(im)
        % no channel contributes
        im = zeros(imSize, 'single');
    end
    im = unmix_channels_sum({}, [], [], constOffset, dtype, im);
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


function [dtype] = frame_data_type(frameFullpath, zarrFile)
% data type of a frame: from the zarr metadata, or from the first slice of a
% tiff with several (faster than reading the tiff's metadata)

if zarrFile
    dtype = char(getImageDataType(frameFullpath, true));
    return;
end
imSize = getImageSize(frameFullpath);
if numel(imSize) > 2 && imSize(3) > 1
    dtype = class(readtiff(frameFullpath, 'range', [1, 1]));
else
    dtype = char(getImageDataType(frameFullpath));
end

end
