function [done_flag] = XR_ExM_PunctaRemoval_block(batchInds, zarrFullpath, outFullpath, ...
    flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin)
% remove puncta for a given block
%
% xruan (02/22/2022) add support for upper bound, and weighted threshold.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('outFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addRequired('RegionBBoxes', @isnumeric);
ip.addRequired('localBBoxes', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('Sigma', 2.5, @isnumeric); % for 3D Gauss filtering
ip.addParameter('OTSUMaxPer', 99.9, @isnumeric); % Max percentile of data for OTSU calculation
ip.addParameter('MinThreshold', 0, @isnumeric); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('MaxThreshold', [], @isnumeric); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('BaseThreshold', [], @isnumeric); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('volThrsh', 1000, @isvector); % set to zere if no "pre-cleaning" is necessary to remove high-freq noise; second value is for the DAN channel
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, outFullpath, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin{:});

pr = ip.Results;
Sigma = pr.Sigma;
OTSUMaxPer = pr.OTSUMaxPer;
MinThreshold = pr.MinThreshold;
MaxThreshold = pr.MaxThreshold;
BaseThreshold = pr.BaseThreshold;
volThrsh = pr.volThrsh;
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

bim = blockedImage(zarrFullpath, 'Adapter', ZarrAdapter);

if ~exist(outFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', outFullpath);
end
nv_bim = blockedImage(outFullpath, 'Adapter', ZarrAdapter);
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
    
    baStart = localBBoxes(i, 1 : 3);
    baEnd = localBBoxes(i, 4 : 6);

    % load the region in input 
    % in_batch = bim.getRegion(ibStart, ibEnd);
    batch_i = bim.Adapter.getIORegion(ibStart, ibEnd);
    if all(batch_i == 0, 'all')
        nv_bim.Adapter.setRegion(obStart, obEnd, zeros(obEnd - obStart + 1, dtype));
        done_flag(i) = true;        
        continue;
    end
    
    % puncta removal    
    im_g = imgaussfilt3(batch_i, Sigma);
    T = thresholdOtsu(batch_i(batch_i > 0 & batch_i < prctile(batch_i(:), OTSUMaxPer)));
    if ~isempty(BaseThreshold)
        % weighted average to consider both local and global effects
        if isempty(T)
            T = MinThreshold;
        end
        if T > MaxThreshold
            T = MaxThreshold;
        end
        T = 0.75 * BaseThreshold + 0.25 * T;
    else
        if isempty(T) || T < MinThreshold
            T = MinThreshold;
        end
    end
        
    im_g = bwareaopen(im_g > T, volThrsh, 26);
    
    batch_i = single(batch_i) .* im_g;
    batch_i = cast(batch_i, dtype);
    batch_i = batch_i(baStart(1) : baEnd(1), baStart(2) : baEnd(2), baStart(3) : baEnd(3));    
    % write out_batch (in the future, directly write the whole region)
    nv_bim.Adapter.setRegion(obStart, obEnd, batch_i)

    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end


