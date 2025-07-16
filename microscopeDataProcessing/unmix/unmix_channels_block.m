function [] = unmix_channels_block(batchInds, zarrFullpaths, unmixFullpath, unmixFactors, flagFullname, BatchBBoxes, varargin)
% unmixing a given block for a channel by subtracting another channel with a factor


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('zarrFullpaths', @(x) iscell(x));
ip.addRequired('unmixFullpath', @(x) ischar(x));
ip.addRequired('unmixFactors', @(x) isnumeric(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpaths, unmixFullpath, unmixFactors, flagFullname, BatchBBoxes, varargin{:});

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

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d ... ', bi);
    tic;
    
    ibStart = BatchBBoxes(i, 1 : 3);
    ibEnd = BatchBBoxes(i, 4 : 6);
    
    % load the region in input 
    in_batch = zeros(ibEnd - ibStart + 1, 'single');
    for c = 1 : nC
        in_batch = in_batch + single(readzarr(zarrFullpaths{c}, 'inputBbox', [ibStart, ibEnd])) * unmixFactors(c);
    end
    in_batch = max(0, in_batch);
    in_batch = cast(in_batch, dtype);
    
    % write the block
    obStart = BatchBBoxes(i, 1 : 3);
    obEnd = BatchBBoxes(i, 4 : 6);    
    writezarr(in_batch, unmixFullpath, 'bbox', [obStart, obEnd], create=false)

    done_flag(i) = true;
    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end


end
