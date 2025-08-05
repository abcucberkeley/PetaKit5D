function [] = resaveZarrBlock(batchInds, zarrFullpath, resultFullpath, flagFullname, BatchBBoxes, RegionBBoxes, varargin)
% resave each block for given block indices for zarr.
% 
% 
% Author: Xiongtao Ruan (12/19/2020) change to the point of view of output,
% especially for batch size and blockSize;


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('resultFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addRequired('RegionBBoxes', @isnumeric);
ip.addParameter('overwrite', false, @islogical);

ip.parse(batchInds, zarrFullpath, resultFullpath, flagFullname, BatchBBoxes, RegionBBoxes, varargin{:});

pr = ip.Results;
overwrite = pr.overwrite;

% we assume the path exists, otherwise return error (in case of completion 
% of processing for all blocks).
flagPath = fileparts(flagFullname);
if ~exist(flagPath, 'dir')
    error('The batch flag directory %s does not exist, skip the processing!', flagPath);
end

if exist(flagFullname, 'file')
    if overwrite
        delete(flagFullname);
    else
        fprintf('The batch files (%d - %d) already exist, skip them!\n', batchInds(1), batchInds(end));
        return;
    end
end

if ~exist(zarrFullpath, 'dir')
    error('The input zarr file %s doesnot exist!', zarrFullpath);
end

if ~exist(resultFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', resultFullpath);
end

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d... ', bi);
    tic;
    ibStart = BatchBBoxes(i, 1 : 3);
    ibEnd = BatchBBoxes(i, 4 : 6);

    obStart = RegionBBoxes(i, 1 : 3);
    obEnd = RegionBBoxes(i, 4 : 6);
    
    % load the region in input 
    out_batch = readzarr(zarrFullpath, 'inputBbox', [ibStart, ibEnd]);

    writezarr(out_batch, resultFullpath, bbox=[obStart, obEnd], create=false);

    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

