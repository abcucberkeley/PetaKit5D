function [] = stitch_organize_block_info(blockInds, taskSize, BlockInfoFullname, PerBlockInfoPath, flagFullname, varargin)
% organize block infor w.r.t. the blocks and directly save them to json files


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('taskSize', @isnumeric);
ip.addRequired('BlockInfoFullname', @(x) ischar(x));
ip.addRequired('PerBlockInfoPath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('BorderSize', [], @isnumeric);
ip.addParameter('uuid', '', @ischar);

ip.parse(blockInds, taskSize, BlockInfoFullname, PerBlockInfoPath, flagFullname, varargin{:});

t0 = tic();

pr = ip.Results;
Overwrite = pr.Overwrite;
BorderSize = pr.BorderSize;
uuid = pr.uuid;
if isempty(uuid)
    uuid = get_uuid();
end

% we assume the path exists, otherwise return error (in case of completion 
% of processing for all blocks).
if exist(flagFullname, 'file')
    if Overwrite
        delete(flagFullname);
    else
        fprintf('The block files (%d - %d) already exist, skip them!\n', blockInds(1), blockInds(end));
        return;
    end
end

a = load(BlockInfoFullname, 'overlap_matrix', 'half_ol_region_cell', 'ol_region_cell', ...
    'overlap_map_mat', 'bboxes', 'bboxStart_mat', 'bboxEnd_mat', 'nvSize', 'blockSize', 'bSubSz');
overlap_matrix = a.overlap_matrix;
half_ol_region_cell = a.half_ol_region_cell;
% ol_region_cell = a.ol_region_cell;
overlap_map_mat = a.overlap_map_mat;
bboxes = a.bboxes;
bboxStart_mat = a.bboxStart_mat;
bboxEnd_mat = a.bboxEnd_mat;
nvSize = a.nvSize;
blockSize = a.blockSize;
bSubSz = a.bSubSz;

nF = size(bboxes, 1);

% compute blockInfo for the included block indices
numBlocks_t = blockInds(end) - blockInds(1) + 1;
[bys, bxs, bzs] = ind2sub(bSubSz, (blockInds(1) : blockInds(end))');

blk_bboxes = [([bys, bxs, bzs] - 1)  .* blockSize + 1, min(nvSize, [bys, bxs, bzs] .* blockSize)];

% map each block with tile ind
cuboid_11 = bboxes(:, 1 : 3);
cuboid_12 = bboxes(:, 4 : 6);

include_tile_mat = false(nF, 1);
blockTileInds = cell(numBlocks_t, 1);
for bidx = 1 : numBlocks_t
    % blockInd = blockInds(1) + bidx - 1;
    cuboid_21 = blk_bboxes(bidx, 1 : 3);
    cuboid_22 = blk_bboxes(bidx, 4 : 6);
    
    ol_1 = cuboid_11 <= cuboid_21 & cuboid_21 <= cuboid_12;
    ol_2 = cuboid_11 <= cuboid_22 & cuboid_22 <= cuboid_12;
    ol_3 = cuboid_21 <= cuboid_12 & cuboid_12 <= cuboid_22;
    ol_4 = cuboid_21 <= cuboid_12 & cuboid_12 <= cuboid_22;

    is_overlap = all(ol_1 | ol_2 | ol_3 | ol_4, 2);
    include_tile_mat = include_tile_mat | is_overlap;
    blockTileInds{bidx} = find(is_overlap);
end

tileBlockInfo = cell(nF, 1);
for i = 1 : nF
    if ~include_tile_mat(i)
        continue;
    end
        
    bbox_i = bboxes(i, :);
    bboxStart = bboxStart_mat(i, :);
    bboxEnd = bboxEnd_mat(i, :);
    if all(bboxStart == 0)
        continue;
    end

    blockInfo = bboxToBlocks(bbox_i, bboxStart, blockSize, nvSize, BorderSize);
    
    % block info for regions not included for the tile
    numMblocks = sum(overlap_matrix(i, :)) + sum(overlap_matrix(:, i));
    mblockInfo_cell = cell(numMblocks, 1);
    mbboxes = zeros(numMblocks, 6);
    counter = 1;
    for j = 1 : nF
        if ~overlap_matrix(i, j) && ~overlap_matrix(j, i)
            continue;
        end

        if i < j
            ind = (j - 1) * nF + i;
            cind = overlap_map_mat(overlap_map_mat(:, 2) == ind, 1);
            mregion = half_ol_region_cell{cind}{1};
        else
            ind = (i - 1) * nF + j;
            cind = overlap_map_mat(overlap_map_mat(:, 2) == ind, 1);            
            mregion = half_ol_region_cell{cind}{2};
        end
        
        midx = mregion;
        % dsr_ol = dsr(mregion_f(2, 1) : mregion_f(2, 2), mregion_f(1, 1) : mregion_f(1, 2), mregion_f(3, 1) : mregion_f(3, 2));
        % dsr(midx(2, 1) : midx(2, 2), midx(1, 1) : midx(1, 2), midx(3, 1) : midx(3, 2)) = 0;

        mbbox = [midx(2, 1), midx(1, 1), midx(3, 1), midx(2, 2), midx(1, 2), midx(3, 2)];
        mbbox = [max(mbbox(1 : 3), bboxStart), min(mbbox(4 : 6), bboxEnd)];
        if any(mbbox(4 : 6) < mbbox(1 : 3))
            continue;
        end
        mbboxStart = mbbox(1 : 3);
        mbbox = mbbox - repmat(bboxStart - bbox_i(1 : 3), 1, 2);
        mblockInfo = bboxToBlocks(mbbox, mbboxStart, blockSize, nvSize, BorderSize);

        mblockInfo_cell{counter} = mblockInfo;
        mbboxes(counter, :) = mbbox;
        counter = counter + 1;
    end
    mblockInfo_cell(counter : end) = [];
    mbboxes(counter : end, :) = [];
    tileBlockInfo{i} = {blockInfo, mblockInfo_cell, mbboxes};
end

stitchBlockInfo = cell(numBlocks_t, 1);

for bidx = 1 : numBlocks_t
    blockInd = blockInds(1) + bidx - 1;
    tileIdx = blockTileInds{bidx};
    blk_bbox = blk_bboxes(bidx, :);
    cuboid_21 = blk_bbox(1 : 3);
    cuboid_22 = blk_bbox(4 : 6);

    for i = 1 : numel(tileIdx)
        ti = tileIdx(i);
        blockInfo = tileBlockInfo{ti}{1};
        mblockInfo_cell = tileBlockInfo{ti}{2};
        mbboxes = tileBlockInfo{ti}{3};

        ib_ind = blockInfo.blockInd == blockInd;

        bCoords = blockInfo.bCoords(ib_ind, :);
        wCoords = blockInfo.wCoords(ib_ind, :);
        bboxCoords = blockInfo.bboxCoords(ib_ind, :);
        
        cuboid_11 = mbboxes(:, 1 : 3);
        cuboid_12 = mbboxes(:, 4 : 6);
        
        ol_1 = cuboid_11 <= cuboid_21 & cuboid_21 <= cuboid_12;
        ol_2 = cuboid_11 <= cuboid_22 & cuboid_22 <= cuboid_12;
        ol_3 = cuboid_21 <= cuboid_12 & cuboid_12 <= cuboid_22;
        ol_4 = cuboid_21 <= cuboid_12 & cuboid_12 <= cuboid_22;
    
        is_overlap = all(ol_1 | ol_2 | ol_3 | ol_4, 2);
        
        mblock_inds = find(is_overlap);
        mblockInfo_j = struct([]);
        for k = 1 : numel(mblock_inds)
            mblockInfo_k = mblockInfo_cell{mblock_inds(k)};
            ind_k = mblockInfo_k.blockInd == blockInd;
            mblockInfo_j(k).bCoords = mblockInfo_k.bCoords(ind_k, :);
            mblockInfo_j(k).wCoords = mblockInfo_k.wCoords(ind_k, :);
            mblockInfo_j(k).bboxCoords = mblockInfo_k.bboxCoords(ind_k, :);
        end
        
        stitchBlockInfo{bidx}(i).tileInd = ti;
        stitchBlockInfo{bidx}(i).bCoords = bCoords;
        stitchBlockInfo{bidx}(i).wCoords = wCoords;
        stitchBlockInfo{bidx}(i).bboxCoords = bboxCoords;
        stitchBlockInfo{bidx}(i).mblockInfo = mblockInfo_j;
        stitchBlockInfo{bidx}(i).borderSize = BorderSize;
    end
end

% save stitchBlockInfo to json file
numTasks = ceil(numBlocks_t / taskSize);

for t = 1 : numTasks
    st = (t - 1) * taskSize + 1;
    tt = min(t * taskSize, numBlocks_t);
    bs = st + blockInds(1) - 1;
    bt = tt + blockInds(1) - 1;

    stitchBlockInfo_t = stitchBlockInfo(st : tt);
    
    PerBlockInfoFullname = sprintf('%s/stitch_block_info_blocks_%d_%d.json', PerBlockInfoPath, bs, bt);
    PerBlockInfoTmppath = sprintf('%s_%s.json', PerBlockInfoFullname(1 : end - 5), uuid);
    % save('-v7.3', PerBlockInfoTmppath, 'blockInds', 'stitchBlockInfo_t');
    s = jsonencode(stitchBlockInfo_t, 'PrettyPrint', true);
    fid = fopen(PerBlockInfoTmppath, 'w');
    fprintf(fid, s);
    fclose(fid);
    movefile(PerBlockInfoTmppath, PerBlockInfoFullname);
end

t1 = toc(t0);
flatTmpFn = sprintf('%s_%s.mat', flagFullname(1 : end - 4), uuid);
save(flatTmpFn, 'blockInds', 't1');
movefile(flatTmpFn, flagFullname);

end

