function [] = compute_tile_bwdist_parser(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, singleDistMap, blockSize, shardSize, compressor, distBbox, Overwrite)


%#function compute_tile_bwdist

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
if ischar(distBbox)
    distBbox = str2num(distBbox);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
end

compute_tile_bwdist(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, singleDistMap, blockSize, shardSize, compressor, distBbox, Overwrite)


end