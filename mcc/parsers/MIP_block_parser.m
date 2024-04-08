function [done_flag] =  MIP_block_parser(batchInds, zarrFullpath, MIPFullpaths, flagFullname, BatchBBoxes, poolSize, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('MIPFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('poolSize', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, MIPFullpaths, flagFullname, BatchBBoxes, poolSize, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
uuid = pr.uuid;
debug = pr.debug;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(MIPFullpaths) && ~isempty(MIPFullpaths) && strcmp(MIPFullpaths(1), '{')
    MIPFullpaths = eval(MIPFullpaths);
end
if ischar(BatchBBoxes)
    BatchBBoxes = str2num(BatchBBoxes);
end
if ischar(poolSize)
    poolSize = str2num(poolSize);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(debug)
    debug = str2num(debug);
end

MIP_block(batchInds, zarrFullpath, MIPFullpaths, flagFullname, BatchBBoxes, ...
    poolSize, Overwrite=Overwrite, uuid=uuid, debug=debug);

end

