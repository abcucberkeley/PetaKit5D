function [done_flag] = XR_crop_block_parser(batchInds, zarrFullpath, cropFullpath, flagFullname, BatchBBoxes, RegionBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('cropFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mipBatchSize', [], @(x) isnumeric(x) || ischar(x)); % save the MIPs of the blocks for MIP slabs with this batch size in the flag
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, cropFullpath, flagFullname, BatchBBoxes, RegionBBoxes, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
mipBatchSize = pr.mipBatchSize;
uuid = pr.uuid;
debug = pr.debug;

if ischar(batchInds)
    batchInds = str2num(batchInds);
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
if ischar(mipBatchSize)
    mipBatchSize = str2num(mipBatchSize);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_crop_block(batchInds, zarrFullpath, cropFullpath, flagFullname, BatchBBoxes, ...
    RegionBBoxes, Overwrite=Overwrite, mipBatchSize=mipBatchSize, uuid=uuid, ...
    debug=debug);

end

