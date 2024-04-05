function [done_flag] = XR_crop_block_parser(batchInds, zarrFullpath, cropFullpath, flagFullname, BatchBBoxes, RegionBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('cropFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, cropFullpath, flagFullname, BatchBBoxes, RegionBBoxes, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
uuid = pr.uuid;
debug = pr.debug;

if ischar(blockInds)
    blockInds = str2num(blockInds);
end
if ischar(BatchBBoxes)
    BatchBBoxes = str2num(BatchBBoxes);
end
if ischar(RegionBBoxes)
    RegionBBoxes = str2num(RegionBBoxes);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_crop_block(blockInds, zarrFullpath, cropFullpath, flagFullname, BatchBBoxes, ...
    RegionBBoxes, Overwrite=Overwrite, uuid=uuid, debug=debug);

end

