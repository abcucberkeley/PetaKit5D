function [] = resampleZarrBlock(batchInds, zarrFullpath, dsFullpath, flagFullname, dsFactor, varargin)
% resample each block for given block indices for zarr.
% 
% 
% Author: Xiongtao Ruan (12/19/2020) change to the point of view of output,
% especially for batch size and blocksize;


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('dsFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('dsFactor', @isnumeric);
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addParameter('bbox', [], @isnumeric);
ip.addParameter('BatchSize', [], @isnumeric);
ip.addParameter('BlockSize', [], @isnumeric);
ip.addParameter('BorderSize', [], @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('Interp', 'linear', @ischar);
% ip.addParameter('imdistPath', '', @ischar); % blurred sigma for blurred blend

ip.parse(batchInds, zarrFullpath, dsFullpath, flagFullname, dsFactor, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
BatchSize = pr.BatchSize;
Interp = pr.Interp;
bbox = pr.bbox;
BlockSize = pr.BlockSize;
BorderSize = pr.BorderSize;

if isempty(BorderSize) 
    BorderSize = [0, 0, 0];
end

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

if ~exist(zarrFullpath, 'dir')
    error('The input zarr file %s doesnot exist!', zarrFullpath);
end

if ~exist(dsFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', dsFullpath);
end

sz = getImageSize(zarrFullpath);
oSz = getImageSize(dsFullpath);
oBlockSize = BlockSize;
if isempty(BatchSize)
    BatchSize = oBlockSize * 2;
end

baSubSz = ceil(oSz ./ BatchSize);

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d... ', bi);
    tic;
    
    [suby, subx, subz] = ind2sub(baSubSz, bi);
    batchSub = [suby, subx, subz];
    
    % get input coordinates
    ibStart_orig = round((batchSub - 1) .* BatchSize .* dsFactor) + 1;
    ibEnd_orig = min(round(batchSub .* BatchSize .* dsFactor), sz);
    obStart = round((ibStart_orig - 1) ./ dsFactor) + 1;
    obStart = round((obStart - 1) ./ BatchSize) .* BatchSize + 1;    

    if ~isempty(bbox)
        ibStart_orig = ibStart_orig + bbox(1 : 3) - 1;
        ibEnd_orig = ibEnd_orig + bbox(1 : 3) - 1;
    end

    ibStart = max(ibStart_orig - BorderSize, 1);
    ibEnd = min(ibEnd_orig + BorderSize, sz);
    
    % load the region in input 
    % in_batch = bim.Adapter.getIORegion(ibStart, ibEnd);
    in_batch = readzarr(zarrFullpath, 'bbox', [ibStart, ibEnd]);
    
    % find the coresponding coordinates in the output
    obEnd = min(obStart + BatchSize - 1, oSz);
    oBatchSize = obEnd - obStart + 1;
    
    % find the start in out batch
    baStart = round((ibStart_orig - ibStart) ./ dsFactor) + 1;
    baEnd = baStart + oBatchSize - 1;
    
    % handle edge blocks
    outSize = max(round(size(in_batch) ./ dsFactor), baEnd);

    % resize batch
    if all(dsFactor == 1)
        out_batch = in_batch;
    else
        out_batch = imresize3(in_batch, outSize, Interp);
    end
    clear in_batch;

    try
        out_batch = crop3d_mex(out_batch, [baStart, baEnd]);
    catch ME
        disp(ME);
        out_batch = out_batch(baStart(1) : baEnd(1), baStart(2) : baEnd(2), baStart(3) : baEnd(3));
    end

    writezarr(out_batch, dsFullpath, bbox=[obStart, obEnd]);

    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

