function [] = processStitchBlock(batchInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchFullname, stitchBlockInfo, tileFns, varargin)
% process each block for given block indices for zarr.
% 
% 
% Author: Xiongtao Ruan (10/04/2020)
% xruan (11/18/2020) change imdistPath to imdistFullpath to enable using
% distance map for the primary channel.
% xruan (02/06/2021): for singleDistMap, directly use the first one (so don't need to repmat). 
% xruan (06/10/2023): change blockSize to batchSize to allow larger batches with smaller chunks 

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('BlockInfoFullname', @(x) ischar(x));
ip.addRequired('PerBlockInfoFullname', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('stitchFullnanme', @(x) ischar(x));
ip.addOptional('stitchBlockInfo', []);
ip.addOptional('tileFns', []);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('imSize', [], @isnumeric);
ip.addParameter('batchSize', [], @isnumeric);
ip.addParameter('dtype', 'uint16', @ischar);
ip.addParameter('BlendMethod', 'mean', @ischar);
ip.addParameter('BorderSize', [], @isnumeric);
ip.addParameter('BlurSigma', 5, @isnumeric); % blurred sigma for blurred blend
ip.addParameter('imdistFullpaths', {}, @iscell); % image distance paths
ip.addParameter('imdistFileIdx', [], @isnumeric); % image distance paths indices
ip.addParameter('poolSize', [], @isnumeric); % distance matrix with max pooling factors
ip.addParameter('weightDegree', 10, @isnumeric); % weight degree for image distances

ip.parse(batchInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchFullname, stitchBlockInfo, tileFns, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
imSize = pr.imSize;
batchSize = pr.batchSize;
dtype = pr.dtype;
BlendMethod = pr.BlendMethod;
BorderSize = pr.BorderSize;
BlurSigma = pr.BlurSigma;
imdistFullpaths = pr.imdistFullpaths;
imdistFileIdx = pr.imdistFileIdx;
poolSize = pr.poolSize;
wd = pr.weightDegree;

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
        fprintf('The block files (%d - %d) already exist, skip them!\n', batchInds(1), batchInds(end));
        return;
    end
end

if ~isempty(PerBlockInfoFullname)
    [~, ~, ext] = fileparts(PerBlockInfoFullname);
    switch ext
        case '.mat'
            a = load(PerBlockInfoFullname, 'stitchBlockInfo_t');
            stitchBlockInfo = a.stitchBlockInfo_t;
        case '.json'
            fid = fopen(PerBlockInfoFullname);
            raw = fread(fid);
            fclose(fid);
            str = char(raw');
            stitchBlockInfo = jsondecode(str);
    end
end
if isstruct(stitchBlockInfo)
    stitchBlockInfo = arrayfun(@(x) stitchBlockInfo(x, :), 1 : size(stitchBlockInfo, 1), 'unif', 0);
end

if ~isempty(BlockInfoFullname)
    a = load(BlockInfoFullname, 'tileFns');
    tileFns = a.tileFns;
end

if isempty(BorderSize) 
    BorderSize = [0, 0, 0];
end

if strcmpi(BlendMethod, 'feather')
    switch numel(poolSize)
        case 3
            psz = [1, 1, poolSize(3)];
        case {6, 9}
            psz = poolSize([4, 5, 3]);
    end
end

bSubSz = ceil(imSize ./ batchSize);

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process batch %d... ', bi);
    tic;
    [suby, subx, subz] = ind2sub(bSubSz, bi);
    batchSub = [suby, subx, subz];
    obStart = (batchSub - 1) .* batchSize + 1;
    obEnd = min(obStart + batchSize - 1, imSize);
    nbsz = obEnd - obStart + 1;

    % stchBlockInfo_i = stitchBlockInfo{bi};
    stchBlockInfo_i = stitchBlockInfo{i};
    numTiles = numel(stchBlockInfo_i);
    
    if numTiles == 0
        nv_block = zeros(nbsz, dtype);
        % writeBlock(nv_bim, blockSub, nv_block, level, Mode);
        % nv_bim.Adapter.setRegion(obStart, obEnd, nv_block);
        writezarr(nv_block, stitchFullname, bbox=[obStart, obEnd]);
        done_flag(i) = true;
        toc;
        continue;
    end
    
    if any(BorderSize > 0)
        bsz = min(batchSize, nbsz) + BorderSize .* ((batchSub ~= 1 & batchSub ~= bSubSz) + 1);
    else
        bsz = min(batchSize, nbsz);
    end
    
    tim_f_block = zeros([bsz, numTiles], dtype);
    if strcmpi(BlendMethod, 'feather')
        if numTiles > 1
            tim_d_block = zeros([bsz, numTiles], 'single');
        end
    else
        tim_block = zeros([bsz, numTiles], dtype);
    end

    % get the pixels for tile in the block
    for j = 1 : numTiles
        tileInd = stchBlockInfo_i(j).tileInd;
        % bim_j = zarrHeaders{tileInd};
        bCoords = stchBlockInfo_i(j).bCoords; 
        bCoords = bCoords(:)';
        bboxCoords = stchBlockInfo_i(j).bboxCoords;
        bboxCoords = bboxCoords(:)';
        % block_j = bim_j.Adapter.getIORegion(bboxCoords(1 : 3), bboxCoords(4 : 6));
        block_j = readzarr(tileFns{tileInd}, bbox=bboxCoords);
        if ~isa(block_j, dtype)
            block_j = cast(block_j, dtype);
        end
        if ~strcmpi(BlendMethod, 'feather')
            block_j_mregion = block_j;
            
            % remove overlap regions
            m_blockInfo = stchBlockInfo_i(j).mblockInfo;
            for k = 1 : numel(m_blockInfo)
                m_bboxCoords = m_blockInfo(k).bboxCoords;
                m_bboxCoords = m_bboxCoords(:)';
                m_bbox = m_bboxCoords - repmat(bboxCoords(1 : 3), 1, 2) + 1;
                block_j_mregion(m_bbox(1) : m_bbox(4), m_bbox(2) : m_bbox(5), m_bbox(3) : m_bbox(6)) = 0;
            end
            
            try
                indexing4d_mex(tim_block, [bCoords(1 : 3), j, bCoords(4 : 6), j], block_j_mregion);
            catch ME
                disp(ME);
                disp(ME.stack);            
                tim_block(bCoords(1) : bCoords(4), bCoords(2) : bCoords(5), bCoords(3) : bCoords(6), j) = block_j_mregion;
            end
        end

        try
            indexing4d_mex(tim_f_block, [bCoords(1 : 3), j, bCoords(4 : 6), j], block_j);
        catch ME
            disp(ME);
            disp(ME.stack);            
            tim_f_block(bCoords(1) : bCoords(4), bCoords(2) : bCoords(5), bCoords(3) : bCoords(6), j) = block_j;
        end

        if numTiles > 1 && strcmpi(BlendMethod, 'feather')
            if numel(imdistFullpaths) == 1
                imdistFullpath = imdistFullpaths{1};                
            else
                imdistFullpath = imdistFullpaths{imdistFileIdx(tileInd)};
            end
            if isempty(poolSize)
                im_d_j = readzarr(imdistFullpath, 'bbox', bboxCoords);
            else
                % for now only consider pooling in z
                dsz = getImageSize(imdistFullpath);
                p_bboxCoords = bboxCoords;
                p_bboxCoords(1 : 3) = max(1, floor(bboxCoords(1 : 3) ./ psz));
                p_bboxCoords(4 : 6) = min(dsz, ceil(bboxCoords(4 : 6) ./ psz));
                
                im_d_j =  readzarr(imdistFullpath, 'bbox', p_bboxCoords);
                try
                    im_d_j = feather_distance_map_resize_3d_mex(im_d_j, bboxCoords(4 : 6) - bboxCoords(1 : 3) + 1, wd);
                catch ME
                    disp(ME);
                    im_d_j =  im_d_j .^ (1 / wd);
                    if bboxCoords(3) == bboxCoords(6)
                        im_d_j = imresize(im_d_j, bboxCoords(4 : 5) - bboxCoords(1 : 2) + 1, 'bilinear');                    
                    else
                        if size(im_d_j, 3) == 1 
                            im_d_j = repmat(im_d_j, 1, 1, 2);
                        end
                        im_d_j = imresize3(im_d_j, bboxCoords(4 : 6) - bboxCoords(1 : 3) + 1, 'linear');
                    end
                    % im_d_j = im_d_j .^ wd;
                    im_d_j = fastPower(im_d_j, wd);
                end
            end
            try
                indexing4d_mex(tim_d_block, [bCoords(1 : 3), j, bCoords(4 : 6), j], im_d_j);
            catch ME
                disp(ME);
                disp(ME.stack);
                tim_d_block(bCoords(1) : bCoords(4), bCoords(2) : bCoords(5), bCoords(3) : bCoords(6), j) = im_d_j;
            end
        end
    end
    
    if numTiles == 1
        % nv_block = tim_block;
        nv_block = tim_f_block;
        if any(BorderSize > 0)
            s = (batchSub ~= 1) .* BorderSize;
            try
                nv_block = crop3d_mex(nv_block, [s + 1, s + batchSize]);
            catch ME
                disp(ME)
                disp(ME.stack);
                nv_block = nv_block(s(1) + 1 : s(1) + batchSize(1), s(2) + 1 : s(2) + batchSize(2), s(3) + 1 : s(3) + batchSize(3));            
            end
        end
        nv_block = cast(nv_block, dtype);
        % writeBlock(nv_bim, blockSub, nv_block, level, Mode);
        if any(nbsz ~= size(nv_block, [1, 2, 3]))
            try
                nv_block = crop3d_mex(nv_block, [1, 1, 1, nbsz]);
            catch ME
                disp(ME)
                disp(ME.stack);
                nv_block = nv_block(1 : nbsz(1), 1 : nbsz(2), 1 : nbsz(3));
            end
        end
        % nv_bim.Adapter.setRegion(obStart, obEnd, nv_block);
        writezarr(nv_block, stitchFullname, bbox=[obStart, obEnd]);
        done_flag(i) = true;
        toc;
        continue;
    end
    
    % convert to single for processing
    if ~strcmp(BlendMethod, 'feather')
        tim_block = single(tim_block);
        tim_f_block = single(tim_f_block);        
    end

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
            % tim_d_block = (tim_d_block / 10) .^ wd;
            % tim_d_block = tim_d_block .^ wd;
            % tim_w_block = tim_d_block .* (tim_block ~= 0);
            % tim_w_block = tim_w_block ./ sum(tim_w_block, 4);
            % nv_block = sum(tim_block .* tim_w_block, 4);
            
            mex_compute = true;
            try
                nv_block = feather_blending_3d_mex(tim_f_block, tim_d_block);  
            catch ME
                disp(ME);
                mex_compute = false;
                tim_w_block = tim_d_block .* (tim_f_block ~= 0); 
                nv_block = sum(single(tim_f_block) .* tim_w_block, 4) ./ sum(tim_w_block, 4);                 
            end
    end
    % clear tim_f_block;
    
    if ~strcmp(BlendMethod, 'none') && ~strcmp(BlendMethod, 'blurred') && ~(strcmp(BlendMethod, 'feather') && mex_compute)
        try 
            nv_block = replace_nan_with_value(nv_block, 0);
        catch ME
            disp(ME);
            nv_block(isnan(nv_block)) = 0;
        end
    end
    if ~strcmp(BlendMethod, 'feather') 
        try 
            nv_f_block = replace_nan_with_value(nv_f_block, 0);
        catch ME
            disp(ME);
            nv_f_block(isnan(nv_f_block)) = 0;
        end        
        nv_zero_inds = nv_block == 0 & nv_f_block ~= 0;
        % nv_block(nv_zero_inds) = nv_f_block(nv_zero_inds);
        nv_block = nv_block .* nv_zero_inds + nv_f_block .* (1 - nv_zero_inds);
    end
        
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
        s = (batchSub ~= 1) .* BorderSize;
        try
            nv_block = crop3d_mex(nv_block, [s + 1, s + batchSize]);
        catch ME
            disp(ME)
            disp(ME.stack);
            nv_block = nv_block(s(1) + 1 : s(1) + batchSize(1), s(2) + 1 : s(2) + batchSize(2), s(3) + 1 : s(3) + batchSize(3));            
        end
    end
    
    % nv_block = cast(nv_block, dtype);
    
    % write the block to zarr file
    % writeBlock(nv_bim, blockSub, nv_block, level, Mode);
    if any(nbsz ~= size(nv_block, [1, 2, 3]))
        try
            nv_block = crop3d_mex(nv_block, [1, 1, 1, nbsz]);
        catch ME
            disp(ME)
            disp(ME.stack);
            nv_block = nv_block(1 : nbsz(1), 1 : nbsz(2), 1 : nbsz(3));
        end
    end    
    % nv_bim.Adapter.setRegion(obStart, obEnd, nv_block);
    writezarr(nv_block, stitchFullname, bbox=[obStart, obEnd]);    
    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fprintf('Processing of batches %d - %d are done!\n', batchInds(1), batchInds(end));
    fclose(fopen(flagFullname, 'w'));
end

end

