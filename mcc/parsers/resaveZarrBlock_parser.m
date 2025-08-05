function [] = resaveZarrBlock_parser(batchInds, zarrFullpath, resultFullpath, flagFullname, BatchBBoxes, RegionBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('resultFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('overwrite', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, resultFullpath, flagFullname, BatchBBoxes, RegionBBoxes, varargin{:});

pr = ip.Results;
overwrite = pr.overwrite;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(BatchBBoxes)
    BatchBBoxes = str2num(BatchBBoxes);
end
if ischar(RegionBBoxes)
    RegionBBoxes = str2num(RegionBBoxes);
end
if ischar(overwrite)
    overwrite = str2num(overwrite);
end

resaveZarrBlock(batchInds, zarrFullpath, resultFullpath, flagFullname, BatchBBoxes, ...
    RegionBBoxes, overwrite=overwrite);

end

