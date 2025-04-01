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

if ~exist(zarrFullpath, 'dir')
    error('The input zarr file %s doesnot exist!', zarrFullpath);
end

if ~exist(dsrFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', dsrFullpath);
end

dtype = getImageDataType(zarrFullpath);

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d ... ', bi);
    tic;
    
    inBbox = inBatchBBoxes(i, :);
    outBbox = outRegionBBoxes(i, :);
    
    % load the region in input 
    in_batch = readzarr(zarrFullpath, 'inputBbox', inBbox);
    if (save16bit && ~strcmp(dtype, 'uint16')) || ~save16bit
        in_batch = single(in_batch);
    end

    % deskew and rotate    
    objectiveScan = false;
    
    dsrBbox = outLocalBboxes(i, :) + [1, outBatchBBoxes(i, 2 : 3), 1, outBatchBBoxes(i, 2 : 3)] - 1;
    
    out_batch = deskewRotateFrame3D(in_batch, SkewAngle, dz, xyPixelSize, ...
                'reverse', Reverse, 'bbox', dsrBbox, 'objectiveScan', objectiveScan, ...
                'resampleFactor', resampleFactor, 'save16bit', save16bit, 'interpMethod', interpMethod);
    clear in_batch;
    if save16bit && ~strcmp(dtype, 'uint16')
        out_batch = uint16(out_batch + 0.5);
    end
    writezarr(out_batch, dsrFullpath, 'bbox', outBbox, create=false);
    clear out_batch;
    done_flag(i) = true;

    toc;
end

if all(done_flag)
    % fclose(fopen(flagFullname, 'w'));
    t = toc(t0);
    save(flagFullname, 't')
end

end

