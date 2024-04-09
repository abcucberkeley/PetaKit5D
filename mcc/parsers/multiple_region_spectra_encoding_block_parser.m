function [] = multiple_region_spectra_encoding_block_parser(batchInds, zarrFullpath, maskFullpath, outFullpath, flagFullname, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('maskFullpath', @(x) ischar(x));
ip.addRequired('outFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addParameter('rescale', false, @(x) islogical(x) || ischar(x));
ip.addParameter('rsRanges', [1, 65535], @(x) isnumeric(x) || ischar(x)); % #region x 2
ip.addParameter('regionInds', [1], @(x) isnumeric(x) || ischar(x)); % #region
ip.addParameter('regionRanges', [1, 65535], @(x) isnumeric(x) || ischar(x)); % #region x 2
ip.addParameter('bbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BatchSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BlockSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BorderSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, maskFullpath, outFullpath, flagFullname, varargin{:});

pr = ip.Results;
rescale = pr.rescale;
rsRanges = pr.rsRanges;
regionInds = pr.regionInds;
regionRanges = pr.regionRanges;
bbox = pr.bbox;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
BorderSize = pr.BorderSize;
Overwrite = pr.Overwrite;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(rescale)
    rescale = str2num(rescale);
end
if ischar(rsRanges)
    rsRanges = str2num(rsRanges);
end
if ischar(regionInds)
    regionInds = str2num(regionInds);
end
if ischar(regionRanges)
    regionRanges = str2num(regionRanges);
end
if ischar(bbox)
    bbox = str2num(bbox);
end
if ischar(BatchSize)
    BatchSize = str2num(BatchSize);
end
if ischar(BlockSize)
    BlockSize = str2num(BlockSize);
end
if ischar(BorderSize)
    BorderSize = str2num(BorderSize);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end

multiple_region_spectra_encoding_block(batchInds, zarrFullpath, maskFullpath, ...
    outFullpath, flagFullname, rescale=rescale, rsRanges=rsRanges, regionInds=regionInds, ...
    regionRanges=regionRanges, bbox=bbox, BatchSize=BatchSize, BlockSize=BlockSize, ...
    BorderSize=BorderSize, Overwrite=Overwrite);

end

