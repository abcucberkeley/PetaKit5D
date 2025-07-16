function [] = unmix_channels_gaussian_block_parser(batchInds, zarrFullpaths, unmixFullpath, ...
    unmixFactors, unmixSigmas, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('unmixFullpath', @(x) ischar(x));
ip.addRequired('unmixFactors', @(x) isnumeric(x) || ischar(x));
ip.addRequired('unmixSigmas', @(x) isnumeric(x) || ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('localBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpaths, unmixFullpath, unmixFactors, unmixSigmas, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
uuid = pr.uuid;
debug = pr.debug;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(zarrFullpaths) && ~isempty(zarrFullpaths) && strcmp(zarrFullpaths(1), '{')
    zarrFullpaths = eval(zarrFullpaths);
end
if ischar(unmixFactors)
    unmixFactors = str2num(unmixFactors);
end
if ischar(unmixSigmas)
    unmixSigmas = str2num(unmixSigmas);
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
if ischar(debug)
    debug = str2num(debug);
end

unmix_channels_gaussian_block(batchInds, zarrFullpaths, unmixFullpath, unmixFactors, ...
    unmixSigmas, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, Overwrite=Overwrite, ...
    uuid=uuid, debug=debug);

end

