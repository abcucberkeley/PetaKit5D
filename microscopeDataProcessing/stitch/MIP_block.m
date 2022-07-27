function [done_flag] =  MIP_block(batchInds, zarrFullpath, MIPFullpaths, flagFullname, BatchBBoxes, bSubs, varargin)
% MIP for all axises for given blocks


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('MIPFullpaths', @(x) iscell(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addRequired('bSubs', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, MIPFullpaths, flagFullname, BatchBBoxes, bSubs, varargin{:});

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

% bim = blockedImage(zarrFullpath, 'Adapter', ZarrAdapter);

% nv_bim_cell = cell(3, 1);
% for i = 1 : 3
%     if ~exist(MIPFullpaths{i}, 'dir')
%         error('The output zarr file %s does not exist!', MIPFullpaths{i});
%     end
%     nv_bim_cell{i} = blockedImage(MIPFullpaths{i}, 'Adapter', ZarrAdapter);
% end

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
    in_batch = readzarr(zarrFullpath, 'bbox', [ibStart, ibEnd]);

    % MIP for each axis
    for j = 1 : 3
        out_batch = max(in_batch, [], j);
        
        obStart = BatchBBoxes(i, 1 : 3);
        obStart(j) = bSubs(i, j);
        obEnd = BatchBBoxes(i, 4 : 6);
        obEnd(j) = bSubs(i, j);
        
        % nv_bim_cell{j}.Adapter.setRegion(obStart, obEnd, out_batch);
        writezarr(out_batch, MIPFullpaths{j}, 'bbox', [obStart, obEnd])
    end

    done_flag(i) = true;
    
    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end


end

