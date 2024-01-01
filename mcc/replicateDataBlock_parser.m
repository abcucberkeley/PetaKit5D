function [] = replicateDataBlock_parser(batchInds, zarrFullpath, repFullpath, flagFullname, RegionBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('repFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, repFullpath, flagFullname, RegionBBoxes, varargin{:});

pr = ip.Results;
uuid = pr.uuid;
debug = pr.debug;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(RegionBBoxes)
    RegionBBoxes = str2num(RegionBBoxes);
end
if ischar(debug)
    debug = strcmp(debug,'true');
end

replicateDataBlock(batchInds, zarrFullpath, repFullpath, flagFullname, RegionBBoxes, ...
    uuid=uuid, debug=debug);

end