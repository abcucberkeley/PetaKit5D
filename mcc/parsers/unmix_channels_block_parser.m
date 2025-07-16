function [] = unmix_channels_block_parser(batchInds, zarrFullpaths, unmixFullpath, unmixFactors, flagFullname, BatchBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('unmixFullpath', @(x) ischar(x));
ip.addRequired('unmixFactors', @(x) isnumeric(x) || ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpaths, unmixFullpath, unmixFactors, flagFullname, BatchBBoxes, varargin{:});

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
if ischar(BatchBBoxes)
    BatchBBoxes = str2num(BatchBBoxes);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(debug)
    debug = str2num(debug);
end

unmix_channels_block(batchInds, zarrFullpaths, unmixFullpath, unmixFactors, ...
    flagFullname, BatchBBoxes, Overwrite=Overwrite, uuid=uuid, debug=debug);

end

