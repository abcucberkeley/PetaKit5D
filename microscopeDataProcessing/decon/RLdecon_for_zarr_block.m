function [] = RLdecon_for_zarr_block(batchInds, zarrFullpath, psfFullpath, deconFullpath, ...
    flagFullname, BatchBBoxes, RegionBBoxes, pixelSize, dz, varargin)
% RL decon for given zarr block/batches and write to the certain output
% location. 
% 
% 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('psfFullpath', @(x) ischar(x));
ip.addRequired('dsFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addRequired('RegionBBoxes', @isnumeric);
ip.addRequired('pixelSize', @isnumeric); %in um
ip.addRequired('dz', @isnumeric); %in um
% ip.addParameter('BlockSize', [], @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('SkewAngle', -32.45 , @isnumeric);
ip.addParameter('flipZstack', false, @islogical); 
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('dzPSF', 0.1 , @isnumeric); %in um
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('scaleFactor', 1e8 , @isnumeric); % scale factor for data
% ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('fixIter', false, @islogical); % CPU Memory in Gb
ip.addParameter('useGPU', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(batchInds, zarrFullpath, psfFullpath, deconFullpath, flagFullname, ...
    BatchBBoxes, RegionBBoxes, pixelSize, dz, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
SkewAngle = pr.SkewAngle;
flipZstack = pr.flipZstack;
Background = pr.Background;
dzPSF = pr.dzPSF;
DeconIter = pr.DeconIter;
scaleFactor = pr.scaleFactor;
useGPU = pr.useGPU;
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

if ~exist(deconFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', deconFullpath);
end
nv_bim = blockedImage(deconFullpath, 'Adapter', ZarrAdapter);

oBlockSize =   nv_bim.BlockSize;
Mode = nv_bim.Mode;
dtype = nv_bim.ClassUnderlying;
level = 1;

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d... ', bi);
    tic;
    
    ibStart = BatchBBoxes(i, 1 : 3);
    ibEnd = BatchBBoxes(i, 4 : 6);
    
    % load the region in input 
    % in_batch = bim.getRegion(ibStart, ibEnd);
    in_batch = bim.Adapter.getIORegion(ibStart, ibEnd);
    
    % deconvolution
    frameFullpath = '';
    deconTmpPath = '';
    Deskew = false;
    Rotate = false;
    Save16bit = false;
    Crop = [];
    zFlip = false;
    GenMaxZproj = [0, 0, 0];
    ResizeImages = false;
    RLMethod = 'simplified';
    fixIter = true;
    errThresh = 1e-12;
    debug = false;
    
    out_batch = RLdecon(frameFullpath, deconTmpPath, psfFullpath, Background, DeconIter, ...
        dzPSF, dz, Deskew, [], SkewAngle, pixelSize, Rotate, Save16bit, Crop, zFlip, ...
        GenMaxZproj, ResizeImages, [], RLMethod, fixIter, errThresh, flipZstack, debug, ...
        'rawdata', in_batch, 'scaleFactor', scaleFactor, 'useGPU', useGPU);
    
    obStart = RegionBBoxes(i, 1 : 3);
    obEnd = RegionBBoxes(i, 4 : 6);
    baStart = obStart - ibStart + 1;
    baEnd = obEnd - ibStart + 1;

    out_batch = out_batch(baStart(1) : baEnd(1), baStart(2) : baEnd(2), baStart(3) : baEnd(3));
    
    % write out_batch (in the future, directly write the whole region)
    % find the block inds in the output file and write each block
%     bSub_s = ceil(obStart ./ oBlockSize);
%     bSub_t = ceil(obEnd ./ oBlockSize);
    
%     [Y, X, Z] = ndgrid(bSub_s(1) : bSub_t(1), bSub_s(2) : bSub_t(2), bSub_s(3) : bSub_t(3));
%     bSubs = [Y(:), X(:), Z(:)];
    
%     for j = 1 : size(bSubs)
%         bSub_j = bSubs(j, :);
%         regionStart = (bSub_j-1) .* oBlockSize + 1;
%         blStart = regionStart - obStart + 1;
%         blEnd = min(blStart + oBlockSize - 1, size(out_batch, [1, 2, 3]));
%         out_block = out_batch(blStart(1) : blEnd(1), blStart(2) : blEnd(2), blStart(3) : blEnd(3));
%         out_block = cast(out_block, dtype);
%         writeZarrBlock(nv_bim, bSub_j, out_block, level, Mode)
%     end
    nv_bim.Adapter.setRegion(obStart, obEnd, out_batch)

    done_flag(i) = true;

    toc;
end

if all(done_flag)
    fclose(fopen(flagFullname, 'w'));
end

end

