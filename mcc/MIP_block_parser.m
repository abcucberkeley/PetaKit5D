function MIP_block_parser(batchInds, zarrFullpath, MIPFullpaths, flagFullname, BatchBBoxes, bSubs, varargin)
% MIP for all axises for given blocks


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('MIPFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('bSubs', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, MIPFullpaths, flagFullname, BatchBBoxes, bSubs, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
uuid = pr.uuid;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(MIPFullpaths)
    MIPFullpaths = eval(MIPFullpaths);
end
if ischar(BatchBBoxes)
    BatchBBoxes = str2num(BatchBBoxes);
end
if ischar(bSubs)
    bSubs = str2num(bSubs);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite,'true');
end

MIP_block(batchInds,zarrFullpath,MIPFullpaths,flagFullname,BatchBBoxes,...
    bSubs,'Overwrite',Overwrite,'uuid',uuid);