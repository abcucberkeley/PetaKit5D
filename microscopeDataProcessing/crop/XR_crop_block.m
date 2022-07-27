function [done_flag] = XR_crop_block(batchInds, zarrFullpath, cropFullpath, flagFullname, BatchBBoxes, RegionBBoxes, varargin)
% crop for give zarr blocks
% 
% if the bbox is out of bound of the input image, it will automatically pad


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('cropFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addRequired('RegionBBoxes', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, cropFullpath, flagFullname, BatchBBoxes, RegionBBoxes, varargin{:});

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

imSize = getImageSize(zarrFullpath);
try 
    nv_bim = blockedImage(cropFullpath, 'Adapter', CZarrAdapter);
catch ME
    disp(ME);
    nv_bim = blockedImage(cropFullpath, 'Adapter', ZarrAdapter);
end    
dtype = nv_bim.ClassUnderlying;

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d ... ', bi);

    tic;
    ibStart = BatchBBoxes(i, 1 : 3);
    ibEnd = BatchBBoxes(i, 4 : 6);

    obStart = RegionBBoxes(i, 1 : 3);
    obEnd = RegionBBoxes(i, 4 : 6);
    
    % load the region in input 
    % in_batch = bim.getRegion(ibStart, ibEnd);
    % in_batch = bim.Adapter.getIORegion(ibStart, ibEnd);
    [is_overlap, cuboid_overlap] = cuboids_overlaps([[1, 1, 1]', imSize'], [ibStart', ibEnd'], false);
    cuboid_overlap = cuboid_overlap(:)';
    if ~is_overlap
        out_batch = zeros(obEnd - obStart, dtype);
    else
        out_batch = readzarr(zarrFullpath, 'bbox', [cuboid_overlap(1 : 3), cuboid_overlap(4 : 6)]);
        if ~all(cuboid_overlap == [ibStart, ibEnd])
            out_batch = padarray(out_batch, cuboid_overlap(1 : 3) - ibStart, 0, 'pre');
            out_batch = padarray(out_batch, ibEnd - cuboid_overlap(4 : 6), 0, 'post');
        end
    end

    writezarr(out_batch, cropFullpath, 'bbox', [obStart, obEnd])

    done_flag(i) = true;
    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end
