function [done_flag] = XR_ExM_PunctaRemoval_interpolated_background_block(batchInds, zarrFullpath, outFullpath, ...
    bgFullpath, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin)
% remove puncta for a given block by removing backgrounds interpolated from
% the estimated ones. 
%

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('outFullpath', @(x) ischar(x));
ip.addRequired('bgFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addRequired('RegionBBoxes', @isnumeric);
ip.addRequired('localBBoxes', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('bgFactor', 1.0, @isnumeric);
ip.addParameter('Sigma', 2.5, @isnumeric); % for 3D Gauss filtering
ip.addParameter('OTSUMaxPer', 99.9, @isnumeric); % Max percentile of data for OTSU calculation
ip.addParameter('volThrsh', 1000, @isvector); % set to zere if no "pre-cleaning" is necessary to remove high-freq noise; second value is for the DAN channel
ip.addParameter('offset', 0, @isnumeric); % offset to add to the cleaned image to make the background non-zero
ip.addParameter('localWinSize', [15, 15, 23], @isvector); % local window size for cropping local region
ip.addParameter('SigmaThrsh', 4, @isnumeric); % sigma threshold for point detection
ip.addParameter('intThrsh', 10000, @isnumeric); % intensity threshold for the peak to be removed
ip.addParameter('initDetect', '2d', @ischar); % initial detection method, 2d mip or 3d stack
ip.addParameter('detVolThrsh', 5000, @isnumeric); % volume threshold for the peak to be removed
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, outFullpath, bgFullpath, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin{:});

pr = ip.Results;
bgFactor = pr.bgFactor;
Sigma = pr.Sigma;
OTSUMaxPer = pr.OTSUMaxPer;
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
    error('The input zarr file %s does not exist!', zarrFullpath);
end

if ~exist(bgFullpath, 'dir')
    error('The background zarr file %s does not exist!', bgFullpath);
end

if ~exist(outFullpath, 'dir')
    error('The output zarr file %s does not exist!', outFullpath);
end

zarrFile = true;
dtype = getImageDataType(outFullpath, zarrFile);

% get the background batch size
tmp = regexp(bgFullpath, 'batch_(\d+)_(\d+)_(\d+)_background_values', 'tokens');
bgBatchSize = cellfun(@(x) str2num(x), tmp{1});
bgSize = getImageSize(bgFullpath);

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
    % batch_i = bim.Adapter.getIORegion(ibStart, ibEnd);
    batch_i = readzarr(zarrFullpath, 'inputBbox', [ibStart, ibEnd]);
    if offset ~= 0
        mask_i = batch_i ~= 0;
    end

    % batch_i = readzarr(zarrFullpath, 'bbox', BatchBBoxes(i, 4 : 6));
    if all(batch_i == 0, 'all')
        % nv_bim.Adapter.setRegion(obStart, obEnd, zeros(obEnd - obStart + 1, dtype));
        writezarr(zeros(obEnd - obStart + 1, dtype), outFullpath, 'bbox', [obStart, obEnd], create=false);
        done_flag(i) = true;        
        continue;
    end
    
    % puncta removal    
    im_g = imgaussfilt3(batch_i, Sigma);

    % load the corresponding background region
    bg_ibStart = floor((ibStart - 1) ./ bgBatchSize) + 1;
    bg_ibEnd = min(bgSize, bg_ibStart + ceil((ibEnd - ibStart + 1) ./ bgBatchSize) - 1);

    bg = readzarr(bgFullpath, inputBbox=[bg_ibStart, bg_ibEnd]);
    im_bg = imresize3(bg, (ibEnd - ibStart + 1), 'linear');
        
    im_g = bwareaopen(im_g > im_bg .* bgFactor, volThrsh, 26);
    
    batch_i = single(batch_i) .* im_g;
    clear im_g;
    
    % puncta removal based on point detection
    if ~false && any(batch_i ~= 0, 'all')
        % remove bright points iteratively
        for j = 1 : numel(SigmaThrsh)
            SigmaThrsh_j = SigmaThrsh(j);
            intThrsh_j = intThrsh(j);
            detVolThrsh_j = detVolThrsh(j);
            batch_i = XR_ExM_PunctaRemovalPointDetection(batch_i, 'Sigma', Sigma, ...
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
    writezarr(batch_i, outFullpath, 'bbox', [obStart, obEnd], create=false);

    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end


