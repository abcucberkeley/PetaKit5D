function [done_flag] = XR_ExM_PunctaRemoval_updated_block(batchInds, zarrFullpath, outFullpath, ...
    flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin)
% remove puncta for a given block
%
% xruan (02/22/2022) add support for upper bound, and weighted threshold.
% xruan (02/23/2022) add support for point detection based puncta removal
% xruan (07/14/2023) change to use mip z to determine the upper bound for
% Otsu thresholding; update point detection code to remove bright puncta


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
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
ip.addParameter('offset', 0, @isnumeric); % offset to add to the cleaned image to make the background non-zero
ip.addParameter('localWinSize', [15, 15, 23], @isvector); % local window size for cropping local region
ip.addParameter('SigmaThrsh', 4, @isnumeric); % sigma threshold for point detection
ip.addParameter('intThrsh', 10000, @isnumeric); % intensity threshold for the peak to be removed
ip.addParameter('initDetect', '2d', @ischar); % initial detection method, 2d mip or 3d stack
ip.addParameter('detVolThrsh', 5000, @isnumeric); % volume threshold for the peak to be removed
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, outFullpath, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin{:});

pr = ip.Results;
Sigma = pr.Sigma;
OTSUMaxPer = pr.OTSUMaxPer;
MinThreshold = pr.MinThreshold;
MaxThreshold = pr.MaxThreshold;
BaseThreshold = pr.BaseThreshold;
volThrsh = pr.volThrsh;
offset = pr.offset;
localWinSize = pr.localWinSize;
SigmaThrsh = pr.SigmaThrsh;
intThrsh = pr.intThrsh;
initDetect = pr.initDetect;
detVolThrsh = pr.detVolThrsh;
Overwrite = pr.Overwrite;

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

if ~exist(outFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', outFullpath);
end

zarrFile = true;
dtype = getImageDataType(outFullpath, zarrFile);

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
    batch_i = readzarr(zarrFullpath, 'bbox', [ibStart, ibEnd]);
    if offset ~= 0
        mask_i = batch_i ~= 0;
    end

    if all(batch_i == 0, 'all')
        writezarr(zeros(obEnd - obStart + 1, dtype), outFullpath, 'bbox', [obStart, obEnd]);
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
        writezarr(zeros(obEnd - obStart + 1, dtype), outFullpath, 'bbox', [obStart, obEnd]);
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
        
    im_g = bwareaopen(im_g > T, volThrsh, 26);
    
    batch_i = single(batch_i) .* im_g;
    clear im_g;
    
    % puncta removal based on point detection
    if ~false && any(batch_i ~= 0, 'all')
        % remove bright points iteratively
        for j = 1 : numel(SigmaThrsh)
            SigmaThrsh_j = SigmaThrsh(j);
            intThrsh_j = intThrsh(j);
            detVolThrsh_j = detVolThrsh(j);
            batch_i = XR_ExM_PunctalRemovalPointDetection(batch_i, 'Sigma', Sigma, ...
                'OTSUMaxPer', OTSUMaxPer,  'localWinSize', localWinSize, 'SigmaThrsh', SigmaThrsh_j, ...
                'intThrsh', intThrsh_j, 'initDetect', initDetect, 'detVolThrsh', detVolThrsh_j);
        end
    end

    if offset ~= 0
        batch_i = batch_i + mask_i * offset;
    end
    batch_i = cast(batch_i, dtype);
    try 
        batch_i = crop3d_mex(batch_i, [baStart, baEnd]);
    catch ME
        disp(ME)
        batch_i = batch_i(baStart(1) : baEnd(1), baStart(2) : baEnd(2), baStart(3) : baEnd(3));
    end
    % write out_batch 
    writezarr(batch_i, outFullpath, 'bbox', [obStart, obEnd]);

    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end


