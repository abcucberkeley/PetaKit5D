function [done_flag] = XR_ExM_PunctaRemoval_background_estimation_block(batchInds, zarrFullpath, resultFullname, BatchBBoxes, varargin)
% estimate background for puncta removal for a given block
%

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('resultFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('Sigma', 2.5, @isnumeric); % for 3D Gauss filtering
ip.addParameter('OTSUMaxPer', 99.9, @isnumeric); % Max percentile of data for OTSU calculation
ip.addParameter('MinThreshold', 0, @isnumeric); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('MaxThreshold', [], @isnumeric); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('BaseThreshold', [], @isnumeric); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, resultFullname, BatchBBoxes, varargin{:});

pr = ip.Results;
Sigma = pr.Sigma;
OTSUMaxPer = pr.OTSUMaxPer;
MinThreshold = pr.MinThreshold;
MaxThreshold = pr.MaxThreshold;
BaseThreshold = pr.BaseThreshold;
Overwrite = pr.Overwrite;
uuid = pr.uuid;

if isempty(uuid)
    uuid = get_uuid();
end

% we assume the path exists, otherwise return error (in case of completion 
% of processing for all blocks).

if exist(resultFullname, 'file')
    if Overwrite
        delete(resultFullname);
    else
        fprintf('The batch files (%d - %d) already exist, skip them!\n', batchInds(1), batchInds(end));
        return;
    end
end

if ~exist(zarrFullpath, 'dir')
    error('The input zarr file %s doesnot exist!', zarrFullpath);
end

done_flag = false(numel(batchInds), 1);
bg_mat = zeros(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d ... ', bi);
    tic;
    
    ibStart = BatchBBoxes(i, 1 : 3);
    ibEnd = BatchBBoxes(i, 4 : 6);
    
    % load the region in input 
    batch_i = readzarr(zarrFullpath, 'inputBbox', [ibStart, ibEnd]);

    if all(batch_i == 0, 'all')
        bg_mat(i) = 0;
        done_flag(i) = true;        
        continue;
    end
    
    % puncta removal    
    im_g = imgaussfilt3(batch_i, Sigma);
    % use MIP to determine the upper bound, which is much faster
    mip_i = max(im_g, [], 3);
    im_g_t = im_g(im_g > 0 & im_g < prctile(mip_i(:), OTSUMaxPer));
    if isempty(im_g_t)
        im_g_t = im_g(im_g > 0 & im_g < prctile(mip_i(mip_i > 0), OTSUMaxPer));
    end
    
    if isempty(im_g_t)
        bg_mat(i) = 0;
        done_flag(i) = true;        
        continue;
    end

    T_g = thresholdOtsu(im_g_t);
    T = T_g;
    % use MIP mask to determine the actual base threshold
    T_min = prctile(mip_i(mip_i > 0), 1);    
    T_med = median(mip_i(mip_i > 0));    
    mip_bg = imgaussfilt(min(mip_i, MaxThreshold), Sigma * 20);
    mip_mask_i = mip_i > mip_bg;
    alpha = 1.1;
    beta = alpha;
    while mean(mip_mask_i, 'all') > 0.9
        mip_mask_i = mip_i > mip_bg * beta;
        beta = beta * alpha;
    end
    mip_mask_i = ~imdilate(mip_mask_i, strel('disk', 5)) & (mip_i > max(BaseThreshold, T_min) & mip_i < max(T_med, MaxThreshold));
    mip_mask_i = imopen(mip_mask_i, strel('disk', 5));
    mip_mask_i = imerode(mip_mask_i, strel('disk', 5));
    T_mask = mode(mip_i(mip_mask_i));
    if T_mask < prctile(mip_i(mip_mask_i), 25)
        T_mask = median(mip_i(mip_mask_i));
    end
    T_mask = max(BaseThreshold, min(T_mask, MaxThreshold));

    if ~isempty(BaseThreshold)
        % weighted average to consider both local and global effects
        if isempty(T)
            T = MinThreshold;
        end
        if T > MaxThreshold
            T = MaxThreshold;
        end
        T = 0.9 * max([BaseThreshold, T_mask]) + 0.1 * T;
    else
        if isempty(T) || T < MinThreshold
            T = MinThreshold;
        end
    end

    bg_mat(i) = T;
        
    done_flag(i) = true;

    toc;
end

if all(done_flag)
    [resultPath, fsn, ext] = fileparts(resultFullname);
    resultTempname = sprintf('%s/%s_%s%s', resultPath, fsn, uuid, ext);
    save('-v7.3', resultTempname, 'bg_mat', 'batchInds');
    movefile(resultTempname, resultFullname);
end

end


