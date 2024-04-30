function [] = processStitchBlock(batchInds, BlockInfoFullname, PerBlockInfoFullname, flagFullname, stitchFullname, stitchBlockInfo, tileFns, varargin)
% process each block for given block indices for zarr.
% 
% 
% Author: Xiongtao Ruan (10/04/2020)
% xruan (11/18/2020) change imdistPath to imdistFullpath to enable using
% distance map for the primary channel.
% xruan (02/06/2021): for singleDistMap, directly use the first one (so don't need to repmat). 
% xruan (06/10/2023): change blockSize to batchSize to allow larger batches with smaller chunks
% xruan (04/18/2024): Optimize feather blending, directly save a region from a tile if the weight is dominant, update mex functions. 
% xruan (04/24/2024): remove blurred blending option. 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('BlockInfoFullname', @(x) ischar(x));
ip.addRequired('PerBlockInfoFullname', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('stitchFullname', @(x) ischar(x));
ip.addOptional('stitchBlockInfo', [], @(x) isempty(x) || isstruct(x));
ip.addOptional('tileFns', [], @(x) isempty(x) || iscell(x));
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('imSize', [], @isnumeric);
ip.addParameter('batchSize', [], @isnumeric);
ip.addParameter('dtype', 'uint16', @ischar);
ip.addParameter('BlendMethod', 'feather', @ischar);
ip.addParameter('BorderSize', [], @isnumeric);
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
imdistFullpaths = pr.imdistFullpaths;
imdistFileIdx = pr.imdistFileIdx;
poolSize = pr.poolSize;
wd = pr.weightDegree;

% we assume the path exists, otherwise return error (in case of completion 
% of processing for all blocks).
flagPath = fileparts(flagFullname);
if ~exist(flagPath, 'dir')
    fprintf('The block directory %s does not exist, skip the processing!\n', flagPath);
    return;
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
    dsz_mat = zeros(numel(imdistFullpaths), 3);
    for i = 1 : numel(imdistFullpaths)
        dsz_mat(i, :) = getImageSize(imdistFullpaths{i});
    end
end

bSubSz = ceil(imSize ./ batchSize);
% feather blending dominant factor
dmtFactor = 1000;

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

    stchBlockInfo_i = stitchBlockInfo{i};
    numTiles = numel(stchBlockInfo_i);
    
    if numTiles == 0
        nv_block = zeros(nbsz, dtype);
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
    
    tileInd_mat = zeros(numTiles, 1);
    bCoords_mat = zeros(numTiles, 6);
    bboxCoords_mat = zeros(numTiles, 6);

    for j = 1 : numTiles
        tileInd = stchBlockInfo_i(j).tileInd;
        tileInd_mat(j) = tileInd;
        bCoords = stchBlockInfo_i(j).bCoords; 
        bCoords = bCoords(:)';
        bCoords_mat(j, :) = bCoords;
        bboxCoords = stchBlockInfo_i(j).bboxCoords;
        bboxCoords = bboxCoords(:)';
        bboxCoords_mat(j, :) = bboxCoords;
    end

    if strcmpi(BlendMethod, 'feather')
        dbsz = bsz;
        if numTiles > 1
            % determine common overlap region
            bCoords_c = [1, 1, 1, bsz];
            if numTiles == 2
                bCoords_c = [max(bCoords_mat(:, 1 : 3)), min(bCoords_mat(:, 4 : 6))];
            elseif numTiles > 2
                [s, t] = find(triu(ones(numTiles), 1));
                [~, cuboid_overlap_mat] = check_cuboids_overlaps(bCoords_mat(s, :), bCoords_mat(t, :), false);
                if any(cuboid_overlap_mat == 0, 'all')
                    cuboid_overlap_mat = cuboid_overlap_mat(~any(cuboid_overlap_mat == 0, 2), :);
                end
                if isempty(cuboid_overlap_mat)
                    cuboid_overlap_mat = bCoords_mat;
                end
                bCoords_c = [min(cuboid_overlap_mat(:, 1 : 3), [], 1), max(cuboid_overlap_mat(:, 4 : 6), [], 1)];
                % include any tiles that does not have the overlap with other tiles (rare but may happen)
                if any(bCoords_mat(:, 1 : 3) > bCoords_c(4 : 6), 'all') || any(bCoords_mat(:, 4 : 6) < bCoords_c(1 : 3), 'all')
                    inds = any(bCoords_mat(:, 1 : 3) > bCoords_c(4 : 6) | bCoords_mat(:, 4 : 6) < bCoords_c(1 : 3), 2);
                    bCoords_c = [min([bCoords_mat(inds, 1 : 3); bCoords_c(1 : 3)], [], 1), max([bCoords_mat(inds, 4 : 6); bCoords_c(4 : 6)], [], 1)];
                end
            end
            dbsz = bCoords_c(4 : 6) - bCoords_c(1 : 3) + 1;
            % define the batch as full if the overlap region is over 80% of the total volume
            is_full_batch = prod(dbsz ./ bsz) >= 0.8;

            if is_full_batch
                bCoords_c = [1, 1, 1, bsz];
                dbsz = bsz;
            end

            if isempty(poolSize) || (~isempty(poolSize) && numTiles == 2)
                tim_d_block = zeros([dbsz, numTiles], 'single');
            else
                tim_d_block_cell = cell(numTiles, 1);
                d_ranges_mat = zeros(numTiles, 2);
                dbsz_mat = zeros(numTiles, 3);
            end
            
            for j = 1 : numTiles
                tileInd = tileInd_mat(j);
                bCoords = bCoords_mat(j, :);
                bboxCoords = bboxCoords_mat(j, :);
                dbCoords = [max(bCoords_c(1 : 3), bCoords(1 : 3)), min(bCoords_c(4 : 6), bCoords(4 : 6))];
                dbboxCoords = bboxCoords - bCoords + dbCoords;
                dbsz_j = dbCoords(4 : 6) - dbCoords(1 : 3) + 1;
                dbsz_mat(j, :) = dbsz_j;

                if numel(imdistFullpaths) == 1
                    imdistFullpath = imdistFullpaths{1};
                    dsz = dsz_mat;
                else
                    imdistFullpath = imdistFullpaths{imdistFileIdx(tileInd)};
                    dsz = dsz_mat(imdistFileIdx(tileInd), :);
                end
                if isempty(poolSize)
                    im_d_j = readzarr(imdistFullpath, 'bbox', dbboxCoords);
                    tim_d_block = indexing4d(tim_d_block, im_d_j, [dbCoords(1 : 3) - bCoords_c(1 : 3) + 1, j, dbCoords(4 : 6) - bCoords_c(1 : 3) + 1, j]);                    
                else
                    % for now only consider pooling in z
                    p_bboxCoords = bboxCoords;
                    p_bboxCoords(1 : 3) = max(1, round(dbboxCoords(1 : 3) ./ psz));
                    p_bboxCoords(4 : 6) = min(dsz, p_bboxCoords(1 : 3) + ceil(dbsz_j ./ psz) - 1);
                    if bboxCoords(6) == dbboxCoords(3)
                        p_bboxCoords(6) = p_bboxCoords(3);
                    end
                    
                    im_d_j =  readzarr(imdistFullpath, 'bbox', p_bboxCoords);
                    if numTiles == 2
                        im_d_j = feather_distance_map_resize_3d(im_d_j, [1, 1, 1, dbsz_j], wd);
                        tim_d_block = indexing4d(tim_d_block, im_d_j, [dbCoords(1 : 3) - bCoords_c(1 : 3) + 1, j, dbCoords(4 : 6) - bCoords_c(1 : 3) + 1, j]);                        
                    else
                        max_d_j = max(im_d_j, [], 'all');
                        min_d_j = min(im_d_j, [], 'all');
                        tim_d_block_cell{j} = im_d_j;
                        d_ranges_mat(j, :) = [max_d_j, min_d_j];
                    end
                end
            end
            if ~isempty(poolSize) && numTiles > 2 
                tim_d_block = zeros([dbsz, numTiles], 'single');
                
                minor_inds = false(numTiles, 1);
                if any(all(dbsz_mat == dbsz, 2))
                    minor_inds = max(d_ranges_mat(all(dbsz_mat == dbsz, 2), 2)) > dmtFactor .^ (1/10) * d_ranges_mat(:, 1);
                end
                % if there is only one major tile, do not load distance map
                if sum(~minor_inds) > 1
                    for j = 1 : numTiles
                        if minor_inds(j)
                            continue;
                        end
                        bCoords = bCoords_mat(j, :);
                        dbCoords = [max(bCoords_c(1 : 3), bCoords(1 : 3)), min(bCoords_c(4 : 6), bCoords(4 : 6))];
                        dbsz_j = dbsz_mat(j, :);
                        im_d_j = tim_d_block_cell{j};
                        im_d_j = feather_distance_map_resize_3d(im_d_j, [1, 1, 1, dbsz_j], wd);
                        tim_d_block = indexing4d(tim_d_block, im_d_j, [dbCoords(1 : 3) - bCoords_c(1 : 3) + 1, j, dbCoords(4 : 6) - bCoords_c(1 : 3) + 1, j]);
                    end
                end
            end
        end

        tim_f_block = zeros([dbsz, numTiles], dtype);
        if numTiles > 1 && ~is_full_batch
            tim_block = zeros([bsz, 1], dtype);
        end
    else
        tim_block = zeros([bsz, numTiles], dtype);
        tim_f_block = zeros([bsz, numTiles], dtype);
    end

    % get the pixels for tile in the block
    for j = 1 : numTiles
        tileInd = tileInd_mat(j);
        bCoords = bCoords_mat(j, :);
        bboxCoords = bboxCoords_mat(j, :);
        
        if strcmpi(BlendMethod, 'feather') && numTiles > 2 && ~isempty(poolSize) && minor_inds(j)
            continue;
        end
        block_j = readzarr(tileFns{tileInd}, bbox=bboxCoords);
        if ~isa(block_j, dtype)
            block_j = cast(block_j, dtype);
        end
        if strcmpi(BlendMethod, 'feather')
            if numTiles > 1
                dbCoords = [max(bCoords_c(1 : 3), bCoords(1 : 3)), min(bCoords_c(4 : 6), bCoords(4 : 6))];
                if ~is_full_batch
                    tim_block = indexing4d(tim_block, block_j, [bCoords(1 : 3), 1, bCoords(4 : 6), 1]);
                    tim_f_block = indexing4d(tim_f_block, block_j, [dbCoords(1 : 3) - bCoords_c(1 : 3) + 1, j, dbCoords(4 : 6) - bCoords_c(1 : 3) + 1, j], ...
                        [dbCoords(1 : 3) - bCoords(1 : 3) + 1, 1, dbCoords(4 : 6) - bCoords(1 : 3) + 1, 1]);
                else
                    tim_f_block = indexing4d(tim_f_block, block_j, [bCoords(1 : 3), j, bCoords(4 : 6), j]);
                end
                continue;
            end
        else
            block_j_mregion = block_j;
            
            % remove overlap regions
            m_blockInfo = stchBlockInfo_i(j).mblockInfo;
            for k = 1 : numel(m_blockInfo)
                m_bboxCoords = m_blockInfo(k).bboxCoords;
                m_bboxCoords = m_bboxCoords(:)';
                m_bbox = m_bboxCoords - repmat(bboxCoords(1 : 3), 1, 2) + 1;
                block_j_mregion(m_bbox(1) : m_bbox(4), m_bbox(2) : m_bbox(5), m_bbox(3) : m_bbox(6)) = 0;
            end
            tim_block = indexing4d(tim_block, block_j_mregion, [bCoords(1 : 3), j, bCoords(4 : 6), j]);
        end
        
        if numTiles == 1 && all(bCoords(4 : 6) - bCoords(1 : 3) + 1 == bsz)
            tim_f_block = block_j;
        else
            tim_f_block = indexing4d(tim_f_block, block_j, [bCoords(1 : 3), j, bCoords(4 : 6), j]);                
        end
    end

    % for feather blending large scale processing (numTile > 2), if any element in the
    % major tiles are zero, load the region of minor tiles
    if strcmpi(BlendMethod, 'feather') && numTiles > 2 && ~isempty(poolSize) && any(minor_inds)
        [major_valid_mat, unvalid_bbox] = check_major_tile_valid(tim_f_block, ~minor_inds);
        % also ensure the major tiles cover regions minor tiles covered
        [major_cover, uncover_bbox_mat] = check_major_tile_cover(bCoords_mat, ~minor_inds);
        orig_major_inds = ~minor_inds;

        if ~major_valid_mat || ~any(major_cover(~minor_inds))
            for j = 1 : numTiles
                if orig_major_inds(j) && sum(orig_major_inds) > 1 && major_cover(j)
                    continue;
                end
                if ~major_cover(j)
                    if ~major_valid_mat
                        uCoords = [min(unvalid_bbox(1 : 3), uncover_bbox_mat(j, 1 : 3)), max(unvalid_bbox(4 : 6), uncover_bbox_mat(j, 4 : 6))];
                    else
                        uCoords = uncover_bbox_mat(j, :);
                    end
                else
                    if ~major_valid_mat
                        uCoords = unvalid_bbox;
                    else
                        uCoords = [];
                    end
                end
                
                if isempty(uCoords)
                    continue;
                end

                tileInd = tileInd_mat(j);
                bCoords = bCoords_mat(j, :);
                bboxCoords = bboxCoords_mat(j, :);

                if any(bCoords(1 : 3) > uCoords(4 : 6) | bCoords(4 : 6) < uCoords(1 : 3))
                    continue;
                end
                
                uCoords = [max(bCoords(1 : 3), uCoords(1 : 3)), min(bCoords(4 : 6), uCoords(4 : 6))];
                uBboxCoords = [bboxCoords(1 :3) - bCoords(1 : 3) + uCoords(1 : 3), bboxCoords(1 :3) - bCoords(1 : 3) + uCoords(4 : 6)];

                dbCoords = [max(bCoords_c(1 : 3), bCoords(1 : 3)), min(bCoords_c(4 : 6), bCoords(4 : 6))];
                dbsz_j = dbCoords(4 : 6) - dbCoords(1 : 3) + 1;
                udbCoords = [max(bCoords_c(1 : 3), uCoords(1 : 3)), min(bCoords_c(4 : 6), uCoords(4 : 6))];
                udbsz_j = udbCoords(4 : 6) - udbCoords(1 : 3) + 1;

                if minor_inds(j)
                    % read the region
                    block_j = readzarr(tileFns{tileInd}, bbox=uBboxCoords);
                    if ~isa(block_j, dtype)
                        block_j = cast(block_j, dtype);
                    end
                    
                    if ~is_full_batch
                        tim_block = indexing4d(tim_block, block_j, [uCoords(1 : 3), 1, uCoords(4 : 6), 1]);
                        tim_f_block = indexing4d(tim_f_block, block_j, [udbCoords(1 : 3) - bCoords_c(1 : 3) + 1, j, udbCoords(4 : 6) - bCoords_c(1 : 3) + 1, j], ...
                            [udbCoords(1 : 3) - uCoords(1 : 3) + 1, 1, udbCoords(4 : 6) - uCoords(1 : 3) + 1, 1]);
                    else
                        tim_f_block = indexing4d(tim_f_block, block_j, [uCoords(1 : 3), j, uCoords(4 : 6), j]);
                    end
    
                    if sum(minor_inds) == 1
                        im_d_j = ones(udbsz_j, 'single') * d_ranges_mat(j, 1)^wd;
                        tim_d_block = indexing4d(tim_d_block, im_d_j, [udbCoords(1 : 3) - bCoords_c(1 : 3) + 1, j, udbCoords(4 : 6) - bCoords_c(1 : 3) + 1, j]);
                        minor_inds(j) = false;
                        continue;
                    end
                end

                % read the distance map and resize
                im_d_j = tim_d_block_cell{j};
                % there will only be one major tile in this case, so direct assign the weight
                if orig_major_inds(j)
                    im_d_j = ones(dbsz_j, 'single') * d_ranges_mat(j, 1)^wd;
                    tim_d_block = indexing4d(tim_d_block, im_d_j, [dbCoords(1 : 3) - bCoords_c(1 : 3) + 1, j, dbCoords(4 : 6) - bCoords_c(1 : 3) + 1, j]);
                else
                    p_udbCoords = udbCoords;
                    p_udbCoords(1 : 3) = max(1, round((udbCoords(1 : 3) - dbCoords(1 : 3) + 1) ./ psz));
                    p_udbCoords(4 : 6) = min(size(im_d_j, 1 : 3), p_udbCoords(1 : 3) + ceil(udbsz_j ./ psz) - 1);
                    if udbCoords(6) == udbCoords(3)
                        p_udbCoords(6) = p_udbCoords(3);
                    end

                    im_d_j = crop3d(im_d_j, p_udbCoords);
                    im_d_j = feather_distance_map_resize_3d(im_d_j, [1, 1, 1, udbsz_j], wd);
                    tim_d_block = indexing4d(tim_d_block, im_d_j, [udbCoords(1 : 3) - bCoords_c(1 : 3) + 1, j, udbCoords(4 : 6) - bCoords_c(1 : 3) + 1, j]);
                end
                minor_inds(j) = false;
            end
        end

        % if there is only one major tile, set numTiles as 1 and directly save it.
        if sum(~minor_inds) == 1 && all(major_cover)
            if ~is_full_batch
                tim_f_block = tim_block;
            else
                major_ind = find(~minor_inds);
                tim_f_block = crop4d(tim_f_block, [1, 1, 1, major_ind, dbsz, major_ind]);
            end
            numTiles = 1;
        end
    end
    
    if numTiles == 1
        if any(BorderSize > 0)
            s = (batchSub ~= 1) .* BorderSize;
            tim_f_block = crop3d(tim_f_block, [s + 1, s + batchSize]);
        end
        tim_f_block = cast(tim_f_block, dtype);
        if any(nbsz ~= size(tim_f_block, [1, 2, 3]))
            tim_f_block = crop3d(tim_f_block, [1, 1, 1, nbsz]);            
        end
        writezarr(tim_f_block, stitchFullname, bbox=[obStart, obEnd]);
        done_flag(i) = true;
        toc;
        continue;
    end
    
    % convert to single for processing
    if ~strcmp(BlendMethod, 'feather')
        tim_block = single(tim_block);
        tim_f_block = single(tim_f_block);        
    end

    if ~strcmp(BlendMethod, 'none') && ~strcmp(BlendMethod, 'feather')
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
        case 'feather'
            if is_full_batch
                [nv_block, mex_compute] = feather_blending_3d(tim_f_block, tim_d_block);
            else
                [nv_block, mex_compute] = feather_blending_3d(tim_f_block, tim_d_block, tim_block, bCoords_c);
            end
        otherwise
            error('Unsupported blending method: %s.', BlendMethod);
    end
    
    if ~strcmp(BlendMethod, 'none') && ~(strcmp(BlendMethod, 'feather') && mex_compute)
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
        nv_block = nv_block .* nv_zero_inds + nv_f_block .* (1 - nv_zero_inds);
    end
            
    if any(BorderSize > 0)
        s = (batchSub ~= 1) .* BorderSize;
        nv_block = crop3d(nv_block, [s + 1, s + batchSize]);
    end
        
    % write the block to zarr file
    if any(nbsz ~= size(nv_block, [1, 2, 3]))
        nv_block = crop3d(nv_block, [1, 1, 1, nbsz]);        
    end    
    writezarr(nv_block, stitchFullname, bbox=[obStart, obEnd]);    
    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fprintf('Processing of batches %d - %d are done!\n', batchInds(1), batchInds(end));
    fclose(fopen(flagFullname, 'w'));
end

end

