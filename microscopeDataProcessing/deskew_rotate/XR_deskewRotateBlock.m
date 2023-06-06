function [done_flag] = XR_deskewRotateBlock(batchInds, zarrFullpath, dsrFullpath, flagFullname, ...
    inBatchBBoxes, outBatchBBoxes, outRegionBBoxes, xyPixelSize, dz, varargin)
% Deskew and/or rotate data for given blocks


t0 = tic;

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('dsrFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('inBatchBBoxes', @isnumeric);
ip.addRequired('outBatchBBoxes', @isnumeric);
ip.addRequired('outRegionBBoxes', @isnumeric);
ip.addRequired('pixelSize', @isnumeric); %in um
ip.addRequired('dz', @isnumeric); %in um
% ip.addParameter('BlockSize', [], @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('SkewAngle', 32.45 , @isnumeric);
ip.addParameter('Reverse', true, @islogical); 
ip.addParameter('flipZstack', false, @islogical); 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x)); % resampling after rotation 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, dsrFullpath, flagFullname, inBatchBBoxes, ...
    outBatchBBoxes, outRegionBBoxes, xyPixelSize, dz, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
flipZstack = pr.flipZstack;
Interp = pr.Interp;
resample = pr.resample;
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

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d ... ', bi);
    tic;
    
    ibStart = inBatchBBoxes(i, 1 : 3);
    ibEnd = inBatchBBoxes(i, 4 : 6);
    
    % load the region in input 
    in_batch = readzarr(zarrFullpath, 'bbox', [ibStart, ibEnd]);
    in_batch = single(in_batch);
    
    % deskew and rotate    
    ObjectiveScan = false;
    
    out_batch = deskewRotateFrame3D(in_batch, SkewAngle, dz, xyPixelSize, ...
                'reverse', Reverse, 'Crop', true, 'ObjectiveScan', ObjectiveScan, ...
                'resample', resample, 'Interp', Interp);
    clear in_batch;
    
    obStart = outRegionBBoxes(i, 1 : 3);
    obEnd = outRegionBBoxes(i, 4 : 6);

    borderSize = [obStart - outBatchBBoxes(i, 1 : 3), outBatchBBoxes(i, 4 : 6) - obEnd];

    if ~isempty(borderSize) && any(borderSize ~= 0)
        try 
            out_batch = crop3d_mex(out_batch, [1 + borderSize(1 : 3), size(out_batch) - borderSize(4 : 6)]);
        catch ME
            disp(ME)
            out_batch = out_batch(borderSize(1) + 1 : end - borderSize(4), borderSize(2) + 1 : end - borderSize(5), ...
                borderSize(3) + 1 : end - borderSize(6)); 
        end
    end
            
    writezarr(out_batch, dsrFullpath, 'bbox', [obStart, obEnd]);
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

