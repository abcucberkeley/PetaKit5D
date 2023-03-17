function [] = resampleZarrBlock_parser(batchInds, zarrFullpath, dsFullpath, flagFullname, dsFactor, varargin)
% resample each block for given block indices for zarr.
% 
% 
% Author: Xiongtao Ruan (12/19/2020) change to the point of view of output,
% especially for batch size and blocksize;


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('dsFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('dsFactor', @(x) isnumeric(x) || ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addParameter('BatchSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Interp', 'linear', @ischar);
ip.addParameter('BlockSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('BorderSize', [], @(x) isnumeric(x) || ischar(x));
% ip.addParameter('imdistPath', '', @ischar); % blurred sigma for blurred blend

ip.parse(batchInds, zarrFullpath, dsFullpath, flagFullname, dsFactor, varargin{:});

Overwrite = ip.Results.Overwrite;
BatchSize = ip.Results.BatchSize;
Interp = ip.Results.Interp;
BlockSize = ip.Results.BlockSize;
BorderSize = ip.Results.BorderSize;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(dsFactor)
    dsFactor = str2num(dsFactor);
end
if ischar(BatchSize)
    BatchSize = str2num(BatchSize);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite,'true');
end
if ischar(BlockSize)
    BlockSize = str2num(BlockSize);
end
if ischar(BorderSize)
    BorderSize = str2num(BorderSize);
end

 resampleZarrBlock(batchInds, zarrFullpath, dsFullpath, flagFullname, dsFactor, ...
     'BatchSize', BatchSize, 'Overwrite', Overwrite, 'Interp', Interp, ...
     'BlockSize', BlockSize, 'BorderSize', BorderSize);