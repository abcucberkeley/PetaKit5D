function [] = unmix_channels_gaussian_block(batchInds, zarrFullpaths, unmixFullpath, ...
    unmixFactors, unmixSigmas, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin)
% unmixing a given block for a channel by subtracting another channel with a factor


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('zarrFullpaths', @(x) iscell(x));
ip.addRequired('unmixFullpath', @(x) ischar(x));
ip.addRequired('unmixFactors', @(x) isnumeric(x));
ip.addRequired('unmixSigmas', @(x) isnumeric(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addRequired('RegionBBoxes', @isnumeric);
ip.addRequired('localBBoxes', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpaths, unmixFullpath, unmixFactors, unmixSigmas, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
uuid = pr.uuid;

% we assume the path exists, otherwise return error (in case of completion 
% of processing for all blocks).
flagPath = fileparts(flagFullname);
if ~exist(flagPath, 'dir')
    error('The block directory %s does not exist, skip the processing!', flagPath);
end

if exist(flagFullname, 'file')
    if Overwrite
        delete(flagFullname);
    else
        fprintf('The batch files (%d - %d) already exist, skip them!\n', batchInds(1), batchInds(end));
        return;
    end
end

if ~exist(zarrFullpaths{1}, 'dir')
    error('The input zarr file %s doesnot exist!', zarrFullpaths{1});
end
nC = numel(zarrFullpaths);
if numel(unmixFactors) ~= nC
    error('The number of unmix factors must equal to the number of channesl!');
end

if ~exist(unmixFullpath, 'dir')
    error('The output zarr file %s does not exist!', unmixFullpath);
end
dtype = getImageDataType(unmixFullpath);

% the mex sums the channels as read in one pass when none is filtered and
% they and the result take at most 8 bytes per voxel (the original code held
% at least a single sum and a single channel); otherwise the channels are
% added one at a time to the sum, which holds no more than the original code
filtered = arrayfun(@(c) ~all(unmixSigmas(c, :) == 0), 1 : nC);
dtypes = cellfun(@(x) char(getImageDataType(x)), zarrFullpaths(:)', 'unif', 0);
atOnce = ~any(filtered) && exist('unmix_channels_sum_mex', 'file') == 3 ...
    && all(ismember([dtypes, {dtype}], {'uint8', 'uint16', 'single', 'double'})) ...
    && sum(cellfun(@dataTypeToByteNumber, dtypes)) + dataTypeToByteNumber(dtype) <= 8;

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d ... ', bi);
    tic;
    
    ibStart = BatchBBoxes(i, 1 : 3);
    ibEnd = BatchBBoxes(i, 4 : 6);
    
    % load the region in input 
    if atOnce
        in_batch = [];
        ims = cell(1, nC);
        for c = 1 : nC
            ims{c} = readzarr(zarrFullpaths{c}, 'inputBbox', [ibStart, ibEnd]);
        end
        in_batch = unmix_channels_sum(ims, unmixFactors, [], [], dtype);
        clear ims;
    else
        in_batch = [];
        for c = 1 : nC
            if filtered(c)
                im_c = single(readzarr(zarrFullpaths{c}, 'inputBbox', [ibStart, ibEnd]));
                im_c = imgaussfilt3(im_c, unmixSigmas(c, :));
            else
                im_c = readzarr(zarrFullpaths{c}, 'inputBbox', [ibStart, ibEnd]);
            end
            in_batch = unmix_channels_sum({im_c}, unmixFactors(c), [], [], '', in_batch);
        end
        clear im_c;
        in_batch = unmix_channels_sum({}, [], [], [], dtype, in_batch);
    end
    
    out_batch = crop3d(in_batch, localBBoxes(i, :));
    % write the block
    obStart = RegionBBoxes(i, 1 : 3);
    obEnd = RegionBBoxes(i, 4 : 6);
    writezarr(out_batch, unmixFullpath, 'bbox', [obStart, obEnd])

    done_flag(i) = true;
    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end
