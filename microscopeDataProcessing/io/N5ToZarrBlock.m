function [done_flag] = N5ToZarrBlock(batchInds, N5Fullpath, outFullpath, flagFullname, inBatchBBoxes, outBatchBBoxes, varargin)
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
ip.addParameter('flipEmptyValue', false, @islogical);

ip.parse(batchInds, N5Fullpath, outFullpath, flagFullname, inBatchBBoxes, outBatchBBoxes, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
flipEmptyValue = pr.flipEmptyValue;

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

if ~exist(N5Fullpath, 'dir')
    error('The input n5 file %s doesnot exist!', N5Fullpath);
end

if ~exist(outFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', outFullpath);
end
bim = blockedImage(N5Fullpath, "Adapter", N5Adapter);
dtype = bim.ClassUnderlying;

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d... ', bi);
    tic;
        
    inBbox = inBatchBBoxes(i, :);
    outBbox = outBatchBBoxes(i, :);

    % load the region in input 
    in_batch = bim.Adapter.getIORegion(inBbox(1 : 3), inBbox(4 : 6));
    if flipEmptyValue
        in_batch = in_batch .* cast(in_batch ~= 65535, dtype);
    end
    writezarr(in_batch, outFullpath, bbox=outBbox);

    done_flag(i) = true;

    toc;
end
bim.Adapter.close();

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

