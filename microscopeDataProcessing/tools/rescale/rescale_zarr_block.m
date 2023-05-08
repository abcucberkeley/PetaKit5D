function [] = rescale_zarr_block(batchInds, zarrFullpath, rsFullpath, flagFullname, rsFactor, rsRange, BatchBBoxes, RegionBBoxes, varargin)
% RL decon for given zarr block/batches and write to the certain output
% location. 
% 
% 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('rsFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('rsFactor', @isnumeric);
ip.addRequired('rsRange', @isnumeric);
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addRequired('RegionBBoxes', @isnumeric);
% ip.addParameter('BlockSize', [], @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
% ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('fixIter', false, @islogical); % CPU Memory in Gb
ip.addParameter('useGPU', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, rsFullpath, flagFullname, rsFactor, rsRange, BatchBBoxes, RegionBBoxes, varargin{:});

pr = ip.Results;
useGPU = pr.useGPU;
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

if ~exist(zarrFullpath, 'dir')
    error('The input zarr file %s doesnot exist!', zarrFullpath);
end

bim = blockedImage(zarrFullpath, 'Adapter', ZarrAdapter);

if ~exist(rsFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', rsFullpath);
end
nv_bim = blockedImage(rsFullpath, 'Adapter', ZarrAdapter);

oBlockSize =   nv_bim.BlockSize;
Mode = nv_bim.Mode;
dtype = nv_bim.ClassUnderlying;
level = 1;

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d... ', bi);
    tic;
    
    ibStart = BatchBBoxes(i, 1 : 3);
    ibEnd = BatchBBoxes(i, 4 : 6);
    
    % load the region in input 
    % in_batch = bim.getRegion(ibStart, ibEnd);
    in_batch = bim.Adapter.getIORegion(ibStart, ibEnd);
    
    % rescale
    out_batch = (single(in_batch) .* rsFactor) .* (rsRange(2) - rsRange(1)) + rsRange(1);
    
    obStart = RegionBBoxes(i, 1 : 3);
    obEnd = RegionBBoxes(i, 4 : 6);
    baStart = obStart - ibStart + 1;
    baEnd = obEnd - ibStart + 1;

    out_batch = out_batch(baStart(1) : baEnd(1), baStart(2) : baEnd(2), baStart(3) : baEnd(3));
    
    % write out_batch (in the future, directly write the whole region)
    % find the block inds in the output file and write each block
%     bSub_s = ceil(obStart ./ oBlockSize);
%     bSub_t = ceil(obEnd ./ oBlockSize);
    
%     [Y, X, Z] = ndgrid(bSub_s(1) : bSub_t(1), bSub_s(2) : bSub_t(2), bSub_s(3) : bSub_t(3));
%     bSubs = [Y(:), X(:), Z(:)];
    
%     for j = 1 : size(bSubs)
%         bSub_j = bSubs(j, :);
%         regionStart = (bSub_j-1) .* oBlockSize + 1;
%         blStart = regionStart - obStart + 1;
%         blEnd = min(blStart + oBlockSize - 1, size(out_batch, [1, 2, 3]));
%         out_block = out_batch(blStart(1) : blEnd(1), blStart(2) : blEnd(2), blStart(3) : blEnd(3));
%         out_block = cast(out_block, dtype);
%         writeZarrBlock(nv_bim, bSub_j, out_block, level, Mode)
%     end
    nv_bim.Adapter.setRegion(obStart, obEnd, out_batch)

    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

