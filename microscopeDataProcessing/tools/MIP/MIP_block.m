function [done_flag] =  MIP_block(batchInds, zarrFullpath, MIPFullpaths, flagFullname, BatchBBoxes, startCoord, poolSize, varargin)
% MIP for all axises for given blocks
% 
% xruan (05/01/2023): add support for user defined max pooling size (only
% in one dimension for MIP slabs for now). 
% xruan (09/06/2023): add support for resampling before max pooling (reducing noise)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('MIPFullpaths', @(x) iscell(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addRequired('startCoord', @isnumeric);
ip.addRequired('poolSize', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, MIPFullpaths, flagFullname, BatchBBoxes, startCoord, poolSize, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
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

poolSize_orig = poolSize;
if numel(poolSize_orig) == 3
    poolSize_orig = [poolSize_orig, 1, 1, 1];
end

% if pool size is smaller than the batch size, use the max pooling routines.
max_pooling = any(BatchBBoxes(1, 4 : 6) - BatchBBoxes(1, 1 : 3) + 1 > poolSize(1 : 3)) || (numel(poolSize) == 6 && any(poolSize(4 : 6) > 1));
resampling = any(BatchBBoxes(1, 4 : 6) - BatchBBoxes(1, 1 : 3) + 1 > poolSize(1 : 3)) && (numel(poolSize) == 9 && any(poolSize(7 : 9) > 1));
poolSize_1 = [1, 1, 1];
dsfactor = [1, 1, 1];
if resampling && numel(poolSize) == 9
    dsfactor = poolSize(7 : 9);
    poolSize(1 : 3) = poolSize(1 : 3) ./ dsfactor;
    poolSize(4 : 6) = poolSize(4 : 6) ./ dsfactor;
end
preMIPCompute = false;
if max_pooling && numel(poolSize) >= 6
    poolSize_1 = poolSize(4 : 6);
    poolSize = poolSize(1 : 3);
    if all(rem(poolSize, poolSize_1) == 0, 'all') && any(poolSize_1 > 1)
        poolSize = poolSize ./ poolSize_1;
        preMIPCompute = true;
    end
end

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d ... ', bi);
    tic;
    
    ibStart = BatchBBoxes(i, 1 : 3);
    ibEnd = BatchBBoxes(i, 4 : 6);
    
    % load the region in input 
    % in_batch = bim.getRegion(ibStart, ibEnd);
    % in_batch = bim.Adapter.getIORegion(ibStart, ibEnd);
    in_batch = readzarr(zarrFullpath, 'inputBbox', [ibStart, ibEnd]);
    
    if resampling
        in_batch = imresize3(in_batch, round(size(in_batch) ./ dsfactor), 'linear');
    end

    if preMIPCompute
        in_batch = max_pooling_3d_mex(in_batch, poolSize_1);
    end

    % MIP for each axis
    for j = 1 : 3
        if preMIPCompute
            poolSize_j = [1, 1, 1];
        else
            poolSize_j = poolSize_1;
        end
        poolSize_j(j) = poolSize(j);
        if max_pooling
            try 
                out_batch = max_pooling_3d_mex(in_batch, poolSize_j);
            catch ME
                disp(ME)
                out_batch = max_pooling_3d(in_batch, poolSize_j);
            end
        else
            out_batch = max(in_batch, [], j);
        end
        
        poolSize_out = poolSize_orig(4 : 6);
        poolSize_out(j) = poolSize_orig(j);

        obStart = BatchBBoxes(i, 1 : 3) - startCoord + 1;
        obStart = floor((obStart - 1) ./ poolSize_out) + 1;
        obEnd = BatchBBoxes(i, 4 : 6) - startCoord + 1;
        obEnd = floor((obEnd - 1) ./ poolSize_out) + 1;

        if any(size(out_batch, 1 : 3) ~= obEnd - obStart + 1)
            warning('The MIP out size does not match that in the bounding box!')
            dsz = obEnd - obStart + 1 - size(out_batch);
            if any(dsz > 0)
                out_batch = padarray(out_batch, max(dsz, 0), 'post', 'replicate');
            end
            if any(dsz < 0)
                out_batch = crop3d(out_batch, [1, 1, 1, obEnd - obStart + 1 ]);
            end
        end
        
        % nv_bim_cell{j}.Adapter.setRegion(obStart, obEnd, out_batch);
        writezarr(out_batch, MIPFullpaths{j}, 'bbox', [obStart, obEnd], create=false)
    end

    done_flag(i) = true;
    
    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

