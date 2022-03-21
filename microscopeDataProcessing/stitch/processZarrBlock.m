function [] = processZarrBlock(batchInds, zarrFullpath, dsFullpath, flagFullname, dsFactor, varargin)
% process each block for given block indices for zarr.
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
ip.addParameter('BatchSize', [], @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('Interp', 'linear', @ischar);
ip.addParameter('BlockSize', [], @isnumeric);
ip.addParameter('BorderSize', [], @isnumeric);
% ip.addParameter('imdistPath', '', @ischar); % blurred sigma for blurred blend

ip.parse(batchInds, zarrFullpath, dsFullpath, flagFullname, dsFactor, varargin{:});

Overwrite = ip.Results.Overwrite;
BatchSize = ip.Results.BatchSize;
Interp = ip.Results.Interp;
BlockSize = ip.Results.BlockSize;
BorderSize = ip.Results.BorderSize;

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

bim = blockedImage(zarrFullpath, 'Adapter', ZarrAdapter);
iSz = bim.Size;

if ~exist(dsFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', dsFullpath);
end
nv_bim = blockedImage(dsFullpath, 'Adapter', ZarrAdapter);

oBlockSize = nv_bim.BlockSize;
oBlockSize = BlockSize;
if isempty(BatchSize)
    BatchSize = oBlockSize * 2;
end

Mode = nv_bim.Mode;
baSubSz = ceil(iSz ./ round((BatchSize .* dsFactor)));
dtype = nv_bim.ClassUnderlying;

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d... ', bi);
    tic;
    
    [suby, subx, subz] = ind2sub(baSubSz, bi);
    batchSub = [suby, subx, subz];
    
    % get input coordinates 
    ibStart_orig = round((batchSub - 1) .* BatchSize .* dsFactor) + 1;
    ibEnd_orig = min(round(batchSub .* BatchSize .* dsFactor), iSz);
    
    if any(BorderSize > 0)
        ibStart = max(ibStart_orig - BorderSize, 1);
        ibEnd = min(ibEnd_orig + BorderSize, iSz);
    end
    
    % load the region in input 
    in_batch = bim.Adapter.getIORegion(ibStart, ibEnd);
    
    % find the coresponding coordinates in the output
    obStart = round((ibStart_orig - 1) ./ dsFactor) + 1;
    obEnd = round(ibEnd_orig ./ dsFactor);
    oBatchSize = obEnd - obStart + 1;
    
    % find the start in out batch
    baStart = round((ibStart_orig - ibStart) ./ dsFactor) + 1;
    baEnd = baStart + oBatchSize - 1;
    
    % handle edge blocks
    outSize = max(round(size(in_batch) ./ dsFactor), baEnd);

    % resize batch
    out_batch_padded = imresize3(in_batch, outSize, Interp);

    out_batch = out_batch_padded(baStart(1) : baEnd(1), baStart(2) : baEnd(2), baStart(3) : baEnd(3));
    
    % write out_batch (in the future, directly write the whole region)
    % find the block inds in the output file and write each block
%     bSub_s = ceil(obStart ./ oBlockSize);
%     bSub_t = ceil(obEnd ./ oBlockSize);
%     
%     [Y, X, Z] = ndgrid(bSub_s(1) : bSub_t(1), bSub_s(2) : bSub_t(2), bSub_s(3) : bSub_t(3));
%     bSubs = [Y(:), X(:), Z(:)];
%     
%     for j = 1 : size(bSubs)
%         bSub_j = bSubs(j, :);
%         regionStart = (bSub_j-1) .* oBlockSize + 1;
%         blStart = regionStart - obStart + 1;
%         blEnd = min(blStart + oBlockSize - 1, size(out_batch, [1, 2, 3]));
%         out_block = out_batch(blStart(1) : blEnd(1), blStart(2) : blEnd(2), blStart(3) : blEnd(3));
%         out_block = cast(out_block, dtype);
%         writeZarrBlock(nv_bim, bSub_j, out_block, level, Mode)
%     end

    nv_bim.Adapter.setRegion(obStart, obEnd, out_batch);

    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

