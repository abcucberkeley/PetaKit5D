function [] = resampleZarrBlock_parser(batchInds, zarrFullpath, dsFullpath, flagFullname, dsFactor, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('dsFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('dsFactor', @(x) isnumeric(x) || ischar(x));
ip.addParameter('bbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BatchSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BlockSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BorderSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Interp', 'linear', @ischar);

ip.parse(batchInds, zarrFullpath, dsFullpath, flagFullname, dsFactor, varargin{:});

pr = ip.Results;
bbox = pr.bbox;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
BorderSize = pr.BorderSize;
Overwrite = pr.Overwrite;
Interp = pr.Interp;

if ischar(blockInds)
    blockInds = str2num(blockInds);
end
if ischar(dsFactor)
    dsFactor = str2num(dsFactor);
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

resampleZarrBlock(blockInds, zarrFullpath, dsFullpath, flagFullname, dsFactor, ...
    bbox=bbox, BatchSize=BatchSize, BlockSize=BlockSize, BorderSize=BorderSize, ...
    Overwrite=Overwrite, Interp=Interp);

end

