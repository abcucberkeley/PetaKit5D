function [done_flag] = XR_ExM_Skeletonization_block_parser(batchInds, zarrFullpath, outFullpath, ...
    flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('outFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('localBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('minBranchLength', 50, @(x) isnumeric(x) || ischar(x));
ip.addParameter('areaThrsh', 5, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, outFullpath, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
minBranchLength = pr.minBranchLength;
areaThrsh = pr.areaThrsh;
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
if ischar(localBBoxes)
    localBBoxes = str2num(localBBoxes);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(minBranchLength)
    minBranchLength = str2num(minBranchLength);
end
if ischar(areaThrsh)
    areaThrsh = str2num(areaThrsh);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_ExM_Skeletonization_block(batchInds, zarrFullpath, outFullpath, flagFullname, ...
    BatchBBoxes, RegionBBoxes, localBBoxes, Overwrite=Overwrite, minBranchLength=minBranchLength, ...
    areaThrsh=areaThrsh, uuid=uuid, debug=debug);

end

