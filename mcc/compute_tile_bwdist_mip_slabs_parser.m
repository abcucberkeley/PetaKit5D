function [] = compute_tile_bwdist_mip_slabs_parser(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, singleDistMap, blockSize, shardSize, compressor, poolSize, distBbox, Overwrite)


if ischar(tileInd)
    tileInd = str2double(tileInd);
end
if ischar(weightDegree)
    weightDegree = str2double(weightDegree);
end
if ischar(singleDistMap)
    singleDistMap = strcmp(singleDistMap, 'true');
end
if ischar(blockSize)
    blockSize = str2num(blockSize);
end
if ischar(shardSize)
    shardSize = str2num(shardSize);
end
if ischar(poolSize)
    poolSize = str2num(poolSize);
end
if ischar(distBbox)
    distBbox = str2num(distBbox);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
end

compute_tile_bwdist_mip_slabs(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, ...
    singleDistMap, blockSize, shardSize, compressor, poolSize, distBbox, Overwrite)


end