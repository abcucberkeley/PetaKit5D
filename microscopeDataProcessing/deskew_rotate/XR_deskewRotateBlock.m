function [done_flag] = XR_deskewRotateBlock(batchInds, zarrFullpath, dsrFullpath, flagFullname, ...
    inBatchBBoxes, outBatchBBoxes, outRegionBBoxes, outLocalBboxes, xyPixelSize, dz, varargin)
% Deskew and/or rotate data for given blocks
%
% xruan (08/02/2023): add support to directly crop the data during deskew/rotate

t0 = tic;

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('dsrFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('inBatchBBoxes', @isnumeric);
ip.addRequired('outBatchBBoxes', @isnumeric);
ip.addRequired('outRegionBBoxes', @isnumeric);
ip.addRequired('outLocalBboxes', @isnumeric);
ip.addRequired('xyPixelSize', @isnumeric); %in um
ip.addRequired('dz', @isnumeric); %in um
% ip.addParameter('blockSize', [], @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('SkewAngle', 32.45 , @isnumeric);
ip.addParameter('Reverse', true, @islogical); 
ip.addParameter('flipZstack', false, @islogical); 
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('interpMethod', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.addParameter('resampleFactor', [], @(x) isempty(x) || isnumeric(x)); % resampling after rotation 
ip.addParameter('mipBatchSize', [], @isnumeric); % save the MIPs of the blocks for MIP slabs with this batch size in the flag
ip.addParameter('inputRows', [], @isnumeric); % [first, last] row of the input region: the batches are resampled on the grid of its result
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, dsrFullpath, flagFullname, inBatchBBoxes, ...
    outBatchBBoxes, outRegionBBoxes, outLocalBboxes, xyPixelSize, dz, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
flipZstack = pr.flipZstack;
save16bit = pr.save16bit;
interpMethod = pr.interpMethod;
resampleFactor = pr.resampleFactor;
mipBatchSize = pr.mipBatchSize;
inputRows = pr.inputRows;
uuid = pr.uuid;

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

% a block of a task whose job ended (e.g., was cancelled): nothing would use
% it (see result_writers)
if result_writers('orphaned', flagFullname)
    return;
end

if ~exist(zarrFullpath, 'dir')
    error('The input zarr file %s doesnot exist!', zarrFullpath);
end

if ~exist(dsrFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', dsrFullpath);
end

dtype = getImageDataType(zarrFullpath);
if ~isempty(mipBatchSize)
    outDtype = getImageDataType(dsrFullpath, true);
    mips = cell(numel(batchInds), 3);
end

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d ... ', bi);
    tic;
    
    inBbox = inBatchBBoxes(i, :);
    outBbox = outRegionBBoxes(i, :);
    
    % load the region in input 
    in_batch = readzarr(zarrFullpath, 'inputBbox', inBbox);
    % (a batch has the whole z range of the input box: flipped after the crop,
    % as in XR_deskewRotateFrame)
    if flipZstack
        in_batch = flip(in_batch, 3);
    end
    if (save16bit && ~strcmp(dtype, 'uint16')) || ~save16bit
        in_batch = single(in_batch);
    end

    % deskew and rotate    
    objectiveScan = false;
    
    dsrBbox = outLocalBboxes(i, :) + [1, outBatchBBoxes(i, 2 : 3), 1, outBatchBBoxes(i, 2 : 3)] - 1;
    % with inputRows, a resampling along y on the grid of the whole result:
    % the rows of the box are the result's, and the batch's input rows are
    % placed in the input region's
    batchRows = [];
    if ~isempty(inputRows)
        dsrBbox([1, 4]) = dsrBbox([1, 4]) + outBatchBBoxes(i, 1) - 1;
        ny = inputRows(2) - inputRows(1) + 1;
        batchRows = [inBbox(1) - inputRows(1) + 1, ny];
        % a last row of the result that samples the region's last row exactly
        % came out as the fill value: imwarp takes a sample past the last row
        % as outside, and the batch's coordinate is past it by its last bit
        % where the whole volume's is not. The last batch gets two copies of
        % that row past it (they only reach samples past the last row)
        nyOut = round(ny / resampleFactor(2));
        if inBbox(4) == inputRows(2) && outRegionBBoxes(i, 4) == nyOut ...
                && abs((nyOut - (nyOut + 1) / 2) * resampleFactor(2) + (ny + 1) / 2 - ny) < 1e-6
            in_batch = cat(1, in_batch, repmat(in_batch(end, :, :), 2, 1, 1));
        end
    end
    
    out_batch = deskewRotateFrame3D(in_batch, SkewAngle, dz, xyPixelSize, ...
                'reverse', Reverse, 'bbox', dsrBbox, 'objectiveScan', objectiveScan, ...
                'resampleFactor', resampleFactor, 'save16bit', save16bit, 'interpMethod', interpMethod, ...
                'batchRows', batchRows);
    clear in_batch;
    if save16bit && ~strcmp(dtype, 'uint16')
        % uint16 rounds to the nearest integer
        out_batch = uint16(out_batch);
    end
    writezarr(out_batch, dsrFullpath, 'bbox', outBbox, create=false);
    if ~isempty(mipBatchSize)
        mips(i, :) = MIP_zarr_block_mips(out_batch, outBbox, mipBatchSize, outDtype);
    end
    clear out_batch;
    done_flag(i) = true;

    toc;
end

if all(done_flag)
    % fclose(fopen(flagFullname, 'w'));
    t = toc(t0);
    if isempty(mipBatchSize)
        save(flagFullname, 't')
    else
        MIP_zarr_save_block_mips(flagFullname, outRegionBBoxes, mips, mipBatchSize, t);
    end
end

end

