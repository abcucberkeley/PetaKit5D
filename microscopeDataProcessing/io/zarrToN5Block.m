function [done_flag] = zarrToN5Block(batchInds, zarrFullpath, outFullpath, flagFullname, inBatchBBoxes, outBatchBBoxes, varargin)
% convert zarr to n5 block by block
%
% Author: Xiongtao Ruan (09/05/2023)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('outFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('inBatchBBoxes', @isnumeric);
ip.addRequired('outBatchBBoxes', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);

ip.parse(batchInds, zarrFullpath, outFullpath, flagFullname, inBatchBBoxes, outBatchBBoxes, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;

% we assume the path exists, otherwise return error (in case of completion 
% of processing for all blocks).
flagPath = fileparts(flagFullname);
if ~exist(flagPath, 'dir')
    error('The block directory %s does not exist, skip the processing!', flagPath);
end

if exist(flagFullname, 'file')
    if Overwrite
        delete(flagFullname);
    else
        fprintf('The batch files (%d - %d) already exist, skip them!\n', batchInds(1), batchInds(end));
        return;
    end
end

if ~exist(zarrFullpath, 'dir')
    error('The input zarr file %s doesnot exist!', zarrFullpath);
end

if ~exist(outFullpath, 'dir')
    error('The output n5 file %s doesnot exist!', outFullpath);
end
bim = blockedImage(outFullpath, "Adapter", N5Adapter);

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d... ', bi);
    tic;
        
    inBbox = inBatchBBoxes(i, :);
    outBbox = outBatchBBoxes(i, :);

    % load the region in input 
    in_batch = readzarr(zarrFullpath, 'bbox', inBbox);
    bim.Adapter.setRegion(outBbox(1 : 3), outBbox(4 : 6), in_batch);

    done_flag(i) = true;

    toc;
end
bim.Adapter.close();

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

