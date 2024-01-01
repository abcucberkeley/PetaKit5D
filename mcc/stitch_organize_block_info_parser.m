function [] = stitch_organize_block_info_parser(blockInds, taskSize, BlockInfoFullname, PerBlockInfoPath, flagFullname, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('taskSize', @(x) isnumeric(x) || ischar(x));
ip.addRequired('BlockInfoFullname', @(x) ischar(x));
ip.addRequired('PerBlockInfoPath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('BorderSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);

ip.parse(blockInds, taskSize, BlockInfoFullname, PerBlockInfoPath, flagFullname, varargin{:});

Overwrite = ip.Results.Overwrite;
BorderSize = ip.Results.BorderSize;
uuid = ip.Results.uuid;

if ischar(blockInds)
    blockInds = str2num(blockInds);
end
if ischar(taskSize)
    taskSize = str2num(taskSize);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite, 'true');
end
if ischar(BorderSize)
    BorderSize = str2num(BorderSize);
end

stitch_organize_block_info(blockInds, taskSize, BlockInfoFullname, PerBlockInfoPath, ...
    flagFullname, Overwrite=Overwrite, BorderSize=BorderSize, uuid=uuid);

end

