function [] = compute_tile_bwdist_parser(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, singleDistMap, blockSize, shardSize, compressor, distBbox, Overwrite)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInfoFullname', @ischar);
ip.addRequired('tileInd', @(x) isscalar(x) || ischar(x));
ip.addRequired('bwdistFullpath', @ischar);
ip.addRequired('weightDegree', @(x) isscalar(x) || ischar(x));
ip.addRequired('singleDistMap', @(x) islogical(x) || ischar(x));
ip.addRequired('blockSize', @(x) isvector(x) || ischar(x));
ip.addRequired('shardSize', @(x) isempty(x) || isvector(x) || ischar(x));
ip.addRequired('compressor', @ischar);
ip.addRequired('poolSize', @(x) isempty(x) || isvector(x) || ischar(x));
ip.addOptional('Overwrite', false, @(x) islogical(x) || ischar(x));

ip.parse(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, singleDistMap, blockSize, shardSize, compressor, distBbox, Overwrite);

pr = ip.Results;

if ischar(tileInd)
    tileInd = str2num(tileInd);
end
if ischar(weightDegree)
    weightDegree = str2num(weightDegree);
end
if ischar(singleDistMap)
    singleDistMap = str2num(singleDistMap);
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
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end

compute_tile_bwdist(blockInfoFullname, tileInd, bwdistFullpath, weightDegree, ...
    singleDistMap, blockSize, shardSize, compressor, poolSize, Overwrite);

end

