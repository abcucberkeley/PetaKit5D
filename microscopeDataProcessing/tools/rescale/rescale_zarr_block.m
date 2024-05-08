function [] = rescale_zarr_block(batchInds, zarrFullpath, rsFullpath, flagFullname, rescaleFactor, rescaleRange, batchBboxes, regionBboxes, varargin)
% rescale zarr block
% 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('rsFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('rescaleFactor', @isnumeric);
ip.addRequired('rescaleRange', @isnumeric);
ip.addRequired('batchBboxes', @isnumeric);
ip.addRequired('regionBboxes', @isnumeric);
% ip.addParameter('BlockSize', [], @isnumeric);
ip.addParameter('overwrite', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, rsFullpath, flagFullname, rescaleFactor, rescaleRange, batchBboxes, regionBboxes, varargin{:});

% we assume the path exists, otherwise return error (in case of completion 
% of processing for all blocks).
flagPath = fileparts(flagFullname);
if ~exist(flagPath, 'dir')
    error('The block directory %s does not exist, skip the processing!', flagPath);
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

if ~exist(rsFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', rsFullpath);
end

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d... ', bi);
    tic;
        
    % load the region in input 
    % in_batch = bim.getRegion(ibStart, ibEnd);
    % in_batch = bim.Adapter.getIORegion(ibStart, ibEnd);
    in_batch = readzarr(zarrFullpath, inputBbox=batchBboxes(i, :));
    
    % rescale
    out_batch = (single(in_batch) .* rescaleFactor) .* (rescaleRange(2) - rescaleRange(1)) + rescaleRange(1);
    
    obStart = regionBboxes(i, 1 : 3);
    obEnd = regionBboxes(i, 4 : 6);
    baStart = obStart - ibStart + 1;
    baEnd = obEnd - ibStart + 1;

    out_batch = out_batch(baStart(1) : baEnd(1), baStart(2) : baEnd(2), baStart(3) : baEnd(3));
    
    % nv_bim.Adapter.setRegion(obStart, obEnd, out_batch)
    writezarr(out_batch, bbox=[obStart, obEnd]);
    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

