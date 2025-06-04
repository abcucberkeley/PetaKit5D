function [] = resampleZarrBlock(batchInds, zarrFullpath, dsFullpath, flagFullname, resampleFactor, varargin)
% resample each block for given block indices for zarr.
% 
% 
% Author: Xiongtao Ruan (12/19/2020) change to the point of view of output,
% especially for batch size and blockSize;


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('dsFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('resampleFactor', @isnumeric);
ip.addParameter('inputBbox', [], @isnumeric);
ip.addParameter('batchSize', [], @isnumeric);
ip.addParameter('blockSize', [], @isnumeric);
ip.addParameter('borderSize', [], @isnumeric);
ip.addParameter('overwrite', false, @islogical);
ip.addParameter('interpMethod', 'linear', @ischar);

ip.parse(batchInds, zarrFullpath, dsFullpath, flagFullname, resampleFactor, varargin{:});

pr = ip.Results;
overwrite = pr.overwrite;
batchSize = pr.batchSize;
interpMethod = pr.interpMethod;
inputBbox = pr.inputBbox;
blockSize = pr.blockSize;
borderSize = pr.borderSize;

if isempty(borderSize) 
    borderSize = [0, 0, 0];
end

% we assume the path exists, otherwise return error (in case of completion 
% of processing for all blocks).
flagPath = fileparts(flagFullname);
if ~exist(flagPath, 'dir')
    error('The block directory %s does not exist, skip the processing!', flagPath);
end

if exist(flagFullname, 'file')
    if overwrite
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
oblockSize = blockSize;
if isempty(batchSize)
    batchSize = oblockSize * 2;
end

baSubSz = ceil(oSz ./ batchSize);

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d... ', bi);
    tic;
    
    [suby, subx, subz] = ind2sub(baSubSz, bi);
    batchSub = [suby, subx, subz];
    
    % get input coordinates
    ibStart_orig = round((batchSub - 1) .* batchSize .* resampleFactor) + 1;
    ibEnd_orig = min(round(batchSub .* batchSize .* resampleFactor), sz);
    obStart = round((ibStart_orig - 1) ./ resampleFactor) + 1;
    obStart = round((obStart - 1) ./ batchSize) .* batchSize + 1;    

    if ~isempty(inputBbox)
        ibStart_orig = ibStart_orig + inputBbox(1 : 3) - 1;
        ibEnd_orig = ibEnd_orig + inputBbox(1 : 3) - 1;
    end

    ibStart = max(ibStart_orig - borderSize, 1);
    ibEnd = min(ibEnd_orig + borderSize, sz);
    
    % load the region in input 
    % in_batch = bim.Adapter.getIORegion(ibStart, ibEnd);
    in_batch = readzarr(zarrFullpath, 'inputBbox', [ibStart, ibEnd]);
    
    % find the coresponding coordinates in the output
    obEnd = min(obStart + batchSize - 1, oSz);
    obatchSize = obEnd - obStart + 1;
    
    % find the start in out batch
    baStart = round((ibStart_orig - ibStart) ./ resampleFactor) + 1;
    baEnd = baStart + obatchSize - 1;
    
    % handle edge blocks
    outSize = max(round(size(in_batch) ./ resampleFactor), baEnd);

    % resize batch
    if all(resampleFactor == 1)
        out_batch = in_batch;
    else
        switch interpMethod
            case {'nearest', 'linear', 'cubic'}
                out_batch = imresize3(in_batch, outSize, interpMethod);
            case 'max'
                out_batch = max_pooling_3d(in_batch, resampleFactor);
            case 'mean'
                out_batch = imresize3_average(in_batch, resampleFactor);
        end
    end
    clear in_batch;

    try
        out_batch = crop3d_mex(out_batch, [baStart, baEnd]);
    catch ME
        disp(ME);
        out_batch = out_batch(baStart(1) : baEnd(1), baStart(2) : baEnd(2), baStart(3) : baEnd(3));
    end

    writezarr(out_batch, dsFullpath, bbox=[obStart, obEnd], create=false);

    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

