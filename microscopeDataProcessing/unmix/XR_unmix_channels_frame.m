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

nC = numel(frameFullpaths);
if nC ~= numel(unmixFactors)
    error('The number of images does not match that of the unmixing factors!');
end

im = zeros(imSize, 'single');
for c = 1 : nC
    if zarrFile
        im_c = readzarr(frameFullpaths{c});
    else
        im_c = readtiff(frameFullpaths{c});
    end
    im_c = single(im_c);
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
