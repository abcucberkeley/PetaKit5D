function [] = processStitchBlock(blockInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchBlockInfo, zarrHeaders, nv_bim, varargin)
% process each block for given block indices for zarr.
% 
% 
% Author: Xiongtao Ruan (10/04/2020)
% xruan (11/18/2020) change imdistPath to imdistFullpath to enable using
% distance map for the primary channel.
% xruan (02/06/2021): for singleDistMap, directly use the first one (so don't need to repmat). 

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('BlockInfoFullname', @(x) ischar(x));
ip.addRequired('PerBlockInfoFullname', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addOptional('stitchBlockInfo', []);
ip.addOptional('zarrHeaders', []);
ip.addOptional('nv_bim', []);
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
% ip.addParameter('stitchInfoDir', 'stitchInfo', @ischar);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('BlendMethod', 'mean', @ischar);
ip.addParameter('BorderSize', [], @isnumeric);
ip.addParameter('BlurSigma', 5, @isnumeric); % blurred sigma for blurred blend
% ip.addParameter('imdistPath', '', @ischar); % blurred sigma for blurred blend
ip.addParameter('imdistFullpaths', {}, @iscell); % image distance paths
ip.addParameter('weightDegree', 10, @isnumeric); % weight degree for image distances

ip.parse(blockInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchBlockInfo, zarrHeaders, nv_bim, varargin{:});

Overwrite = ip.Results.Overwrite;
BlendMethod = ip.Results.BlendMethod;
BorderSize = ip.Results.BorderSize;
BlurSigma = ip.Results.BlurSigma;
imdistFullpaths = ip.Results.imdistFullpaths;
wd = ip.Results.weightDegree;

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
        fprintf('The block files (%d - %d) already exist, skip them!\n', blockInds(1), blockInds(end));
        return;
    end
end

if ~isempty(PerBlockInfoFullname)
    a = load(PerBlockInfoFullname, 'stitchBlockInfo_t');
    stitchBlockInfo = a.stitchBlockInfo_t;
end

if ~isempty(BlockInfoFullname)
    a = load(BlockInfoFullname, 'nv_bim', 'zarrHeaders');
    nv_bim = a.nv_bim;
    nv_bim.Adapter.getInfo;
    
    if isempty(PerBlockInfoFullname)
        a = load(BlockInfoFullname, 'stitchBlockInfo');
        stitchBlockInfo = a.stitchBlockInfo(blockInds);
    end

    zarrHeaders = a.zarrHeaders;
end

Mode = nv_bim.Mode;
imSize = nv_bim.Size;
bSubSz = nv_bim.SizeInBlocks;
blockSize = nv_bim.BlockSize;
dtype = nv_bim.ClassUnderlying;
level = 1;

done_flag = false(numel(blockInds), 1);
for i = 1 : numel(blockInds)
    bi = blockInds(i);
    fprintf('Process block %d... ', bi);
    tic;
    [suby, subx, subz] = ind2sub(bSubSz, bi);
    blockSub = [suby, subx, subz];
    obStart = (blockSub - 1) .* blockSize + 1;
    obEnd = min(obStart + blockSize - 1, imSize);
    nbsz = obEnd - obStart + 1;

    % stchBlockInfo_i = stitchBlockInfo{bi};
    stchBlockInfo_i = stitchBlockInfo{i};
    numTiles = numel(stchBlockInfo_i);
    
    if numTiles == 0
        nv_block = zeros(nbsz, dtype);
        % writeBlock(nv_bim, blockSub, nv_block, level, Mode);
        nv_bim.Adapter.setRegion(obStart, obEnd, nv_block)
        done_flag(i) = true;
        toc;
        continue;
    end
    
    if any(BorderSize > 0)
        bsz = blockSize + BorderSize .* ((blockSub ~= 1 & blockSub ~= bSubSz) + 1);
    else
        bsz = blockSize;
    end
    
    tim_block = zeros([bsz, numTiles], dtype);
    tim_f_block = zeros([bsz, numTiles], dtype);
    
    if strcmpi(BlendMethod, 'feather')
        tim_d_block = zeros([bsz, numTiles], 'single');
    end

    % get the pixels for tile in the block
    for j = 1 : numTiles
        tileInd = stchBlockInfo_i{j}.tileInd;
        bim_j = zarrHeaders{tileInd};
        bCoords = stchBlockInfo_i{j}.bCoords; 
        bboxCoords = stchBlockInfo_i{j}.bboxCoords; 
        block_j = bim_j.Adapter.getIORegion(bboxCoords(1 : 3), bboxCoords(4 : 6));
        if ~isa(block_j, dtype)
            block_j = cast(block_j, dtype);
        end
        block_j_mregion = block_j;
        
        % remove overlap regions
        m_blockInfo = stchBlockInfo_i{j}.mblockInfo;
        for k = 1 : numel(m_blockInfo)
            m_bboxCoords = m_blockInfo{k}.bboxCoords;
            m_bbox = m_bboxCoords - repmat(bboxCoords(1 : 3), 1, 2) + 1;
            block_j_mregion(m_bbox(1) : m_bbox(4), m_bbox(2) : m_bbox(5), m_bbox(3) : m_bbox(6)) = 0;
        end
        
        try
            indexing4d_mex(tim_block, [bCoords(1 : 3), j, bCoords(4 : 6), j], block_j_mregion);
            indexing4d_mex(tim_f_block, [bCoords(1 : 3), j, bCoords(4 : 6), j], block_j);
        catch ME
            disp(ME);
            disp(ME.stack);            
            tim_block(bCoords(1) : bCoords(4), bCoords(2) : bCoords(5), bCoords(3) : bCoords(6), j) = block_j_mregion;
            tim_f_block(bCoords(1) : bCoords(4), bCoords(2) : bCoords(5), bCoords(3) : bCoords(6), j) = block_j;
        end

        if numTiles > 1 && strcmpi(BlendMethod, 'feather')
            % [~, fsname] = fileparts(bim_j.Source);
            if numel(imdistFullpaths) == 1
                imdistFullpath = imdistFullpaths{1};                
            else
                imdistFullpath = imdistFullpaths{tileInd};
            end
            % bim_d_j = blockedImage(imdistFullpath, 'Adapter', ZarrAdapter);
            % block_d_j = bim_d_j.Adapter.getIORegion(bboxCoords(1 : 3), bboxCoords(4 : 6));
            try
                indexing4d_mex(tim_d_block, [bCoords(1 : 3), j, bCoords(4 : 6), j], readzarr(imdistFullpath, 'bbox', bboxCoords));
            catch ME
                disp(ME);
                disp(ME.stack);
                tim_d_block(bCoords(1) : bCoords(4), bCoords(2) : bCoords(5), bCoords(3) : bCoords(6), j) = readzarr(imdistFullpath, 'bbox', bboxCoords);
            end
        end
    end
    
    if numTiles == 1
        nv_block = tim_block;
        if any(BorderSize > 0)
            s = (blockSub ~= 1) .* BorderSize;
            nv_block = nv_block(s(1) + 1 : s(1) + blockSize(1), s(2) + 1 : s(2) + blockSize(2), s(3) + 1 : s(3) + blockSize(3));
        end
        nv_block = cast(nv_block, dtype);
        % writeBlock(nv_bim, blockSub, nv_block, level, Mode);
        if any(nbsz ~= size(nv_block, [1, 2, 3]))
            nv_block = nv_block(1 : nbsz(1), 1 : nbsz(2), 1 : nbsz(3));
        end
        nv_bim.Adapter.setRegion(obStart, obEnd, nv_block);
        done_flag(i) = true;
        toc;
        continue;
    end
    
    % convert to single for processing
    tim_block = single(tim_block);
    tim_f_block = single(tim_f_block);
    
    if ~strcmp(BlendMethod, 'none') && ~strcmp(BlendMethod, 'blurred') && ~strcmp(BlendMethod, 'feather')
        tim_block(tim_block == 0) = nan; 
        tim_f_block(tim_f_block == 0) = nan;
    end
    
    % merge voxel info for the tiles
    switch BlendMethod
        case 'none'
            nv_block = zeros(bsz, 'single');
            nv_f_block = zeros(bsz, 'single');
            for j = 1 : numTiles
                if j == 1
                    nv_block = tim_block(:, :, :, j);
                    nv_f_block = tim_f_block(:, :, :, j);
                else
                    nv_block = nv_block + tim_block(:, :, :, j) .* (nv_block == 0);
                    nv_f_block = nv_f_block + tim_f_block(:, :, :, j) .* (nv_f_block == 0);
                end
            end
        case 'mean'
            nv_block = nanmean(tim_block, 4);
            nv_f_block = nanmean(tim_f_block, 4);
        case 'max'
            nv_block = max(tim_block, [], 4);
            nv_f_block = max(tim_f_block, [], 4);
        case 'median'
            nv_block = nanmedian(tim_block, 4);
            nv_f_block = nanmedian(tim_f_block, 4);
        case 'blurred'
            if numTiles == 1
                nv_block = tim_block(:, :, :, 1);
                nv_f_block = tim_f_block(:, :, :, 1);  
            else
                nv_block = zeros(bsz, 'single');
                nv_f_block = zeros(bsz, 'single');                
                nv_ind_block = zeros(bsz, 'single');
                nv_ind_f_block = zeros(bsz, 'single');
                for j = 1 : numTiles
                    nv_ind_block(nv_block == 0 & tim_block(:, :, :, j) ~= 0) = j; 
                    nv_block = nv_block + tim_block(:, :, :, j) .* (nv_block == 0);
                    nv_ind_f_block(nv_f_block == 0 & tim_f_block(:, :, :, j) ~= 0) = j;
                    nv_f_block = nv_f_block + tim_f_block(:, :, :, j) .* (nv_f_block == 0);
                end
            end
        case 'feather'
            if numTiles == 1
                nv_block = tim_block(:, :, :, 1);
                nv_f_block = tim_f_block(:, :, :, 1);  
            else
                % tim_d_block = (tim_d_block / 10) .^ wd;
                % tim_d_block = tim_d_block .^ wd;
                % tim_w_block = tim_d_block .* (tim_block ~= 0);
                % tim_w_block = tim_w_block ./ sum(tim_w_block, 4);
                % nv_block = sum(tim_block .* tim_w_block, 4);
                
                tim_w_block = tim_d_block .* (tim_f_block ~= 0);
                % tim_w_block = tim_w_block ./ sum(tim_w_block, 4);
                % nv_f_block = sum(tim_f_block .* tim_w_block, 4);
                nv_f_block = sum(tim_f_block .* tim_w_block, 4) ./ sum(tim_w_block, 4); 
                nv_block = nv_f_block;
            end
    end
    clear tim_block tim_f_block;
    
    if ~strcmp(BlendMethod, 'none') && ~strcmp(BlendMethod, 'blurred')
        try 
            nv_block = replace_nan_with_value(nv_block, 0);
            nv_f_block = replace_nan_with_value(nv_f_block, 0);
        catch ME
            disp(ME);
            nv_block(isnan(nv_block)) = 0;
            nv_f_block(isnan(nv_f_block)) = 0;
        end
    end

    nv_zero_inds = nv_block == 0 & nv_f_block ~= 0;
    nv_block(nv_zero_inds) = nv_f_block(nv_zero_inds);
    
    if strcmp(BlendMethod, 'blurred') && numTiles > 1
        nv_ind_block(nv_zero_inds) = nv_ind_f_block(nv_zero_inds);
        nv_ol_regions = false(bsz);
        nv_ind_block_bw = nv_ind_block > 0;
        ksz = max(BorderSize) * 2 + 1;        
        for j = 1 : numTiles
            nv_ol_j = nv_ind_block == j;
            nv_ol_j = imdilate(nv_ol_j, strel('cube', ksz));
            nv_ol_regions = nv_ol_regions | (nv_ol_j & nv_ind_block_bw & (nv_ind_block ~= j));
        end
        nv_block_blur = imgaussfilt3(nv_block, BlurSigma);
        nv_block(nv_ol_regions) = nv_block_blur(nv_ol_regions);
    end
    
    if any(BorderSize > 0)
        s = (blockSub ~= 1) .* BorderSize;
        nv_block = nv_block(s(1) + 1 : s(1) + blockSize(1), s(2) + 1 : s(2) + blockSize(2), s(3) + 1 : s(3) + blockSize(3));
    end
    
    nv_block = cast(nv_block, dtype);
    
    % write the block to zarr file
    % writeBlock(nv_bim, blockSub, nv_block, level, Mode);
    if any(nbsz ~= size(nv_block, [1, 2, 3]))
        nv_block = nv_block(1 : nbsz(1), 1 : nbsz(2), 1 : nbsz(3));
    end    
    nv_bim.Adapter.setRegion(obStart, obEnd, nv_block)
    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

