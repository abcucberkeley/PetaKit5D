function [] = RLdecon_for_zarr_block(batchInds, zarrFullpath, psfFullpath, deconFullpath, ...
    flagFullname, BatchBBoxes, RegionBBoxes, xyPixelSize, dz, varargin)
% RL decon for given zarr block/batches and write to the certain output
% location. 
% 
% xruan (05/24/2022): add support for masked decon, directly skip empty
% regions if the region is empty in the masks. Masks are 2D projections.

t0 = tic;

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @isnumeric);
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('psfFullpath', @(x) ischar(x));
ip.addRequired('deconFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('BatchBBoxes', @isnumeric);
ip.addRequired('RegionBBoxes', @isnumeric);
ip.addRequired('pixelSize', @isnumeric); %in um
ip.addRequired('dz', @isnumeric); %in um
% ip.addParameter('BlockSize', [], @isnumeric);
ip.addParameter('Save16bit', false , @islogical);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('SkewAngle', -32.45 , @isnumeric);
ip.addParameter('flipZstack', false, @islogical); 
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('dzPSF', 0.1 , @isnumeric); %in um
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('scaleFactor', 1e8 , @isnumeric); % scale factor for data
ip.addParameter('deconMaskFns', {} , @iscell); % Full paths of 2D mask zarr files, in xy, xz, yz order
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @isnumeric);
ip.addParameter('OTFCumThresh', 0.9, @isnumeric); % OTF cumutative sum threshold
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x)); % decon in skewed space
ip.addParameter('fixIter', false, @islogical); % CPU Memory in Gb
ip.addParameter('useGPU', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('psfGen', true, @islogical); % psf generation

ip.parse(batchInds, zarrFullpath, psfFullpath, deconFullpath, flagFullname, ...
    BatchBBoxes, RegionBBoxes, xyPixelSize, dz, varargin{:});

pr = ip.Results;
Save16bit = pr.Save16bit;
Overwrite = pr.Overwrite;
SkewAngle = pr.SkewAngle;
flipZstack = pr.flipZstack;
Background = pr.Background;
dzPSF = pr.dzPSF;
DeconIter = pr.DeconIter;
scaleFactor = pr.scaleFactor;
deconMaskFns = pr.deconMaskFns;
RLMethod = pr.RLMethod;
skewed = pr.skewed;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
useGPU = pr.useGPU;
uuid = pr.uuid;
psfGen = pr.psfGen;

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

if ~exist(deconFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', deconFullpath);
end
try 
    nv_bim = blockedImage(deconFullpath, 'Adapter', CZarrAdapter);
catch 
    nv_bim = blockedImage(deconFullpath, 'Adapter', ZarrAdapter);
end
dtype = nv_bim.ClassUnderlying;
finds = [1, 2; 2, 3; 1, 3];

done_flag = false(numel(batchInds), 1);
for i = 1 : numel(batchInds)
    bi = batchInds(i);
    fprintf('Process Batch %d... ', bi);
    tic;
    
    ibStart = BatchBBoxes(i, 1 : 3);
    ibEnd = BatchBBoxes(i, 4 : 6);
    obStart = RegionBBoxes(i, 1 : 3);
    obEnd = RegionBBoxes(i, 4 : 6);

    % use masks to determine whether to run the decon or directly save an empty region
    if ~(isempty(deconMaskFns) || isempty(deconMaskFns{1}))
        skipDecon = false;
        for f = 1 : 3
            im_f = readzarr(deconMaskFns{f}, 'bbox', [obStart(finds(f, :)), 1, obEnd(finds(f, :)), 1]);
            if ~any(im_f(:))
                skipDecon = true;
                break;
            end
        end
        
        if skipDecon
            % nv_bim.Adapter.setRegion(obStart, obEnd, zeros(obEnd - obStart + 1, dtype));
            writezarr(zeros(obEnd - obStart + 1, dtype), deconFullpath, 'bbox', [obStart, obEnd]);
            done_flag(i) = true;
            toc;
            continue;
        end
    end

    % load the region in input 
    in_batch = readzarr(zarrFullpath, 'bbox', [ibStart, ibEnd]);

    % deconvolution
    inputFn = '';
    outputFn = deconFullpath;
    Deskew = false;
    Rotate = false;
    DSRCombined = true;
    Reverse = true;
    fixIter = true;
    debug = false;
    save3Dstack = [false, false, false];
    mipAxis = [0, 0, 0];

    baStart = obStart - ibStart + 1;
    baEnd = obEnd - ibStart + 1;
    deconBbox = [baStart, baEnd];

    out_batch = RLdecon(inputFn, outputFn, psfFullpath, xyPixelSize, dz, dzPSF, ...
        'rawdata', in_batch, 'Save16bit', Save16bit, 'SkewAngle', SkewAngle, ...
        'Deskew', Deskew, 'Rotate', Rotate, 'DSRCombined', DSRCombined, 'Reverse', Reverse, ...
        'Background', Background, 'DeconIter', DeconIter, 'RLMethod', RLMethod, ...
        'skewed', skewed, 'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, ...
        'fixIter', fixIter, 'scaleFactor', scaleFactor, 'deconBbox', deconBbox, ...
        'useGPU', useGPU, 'psfGen', psfGen, 'debug', debug, 'save3Dstack', save3Dstack, ...
        'mipAxis', mipAxis);
    
    clear in_batch;

    writezarr(out_batch, deconFullpath, 'bbox', [obStart, obEnd]);
    done_flag(i) = true;
    toc;
end

if all(done_flag)
    % fclose(fopen(flagFullname, 'w'));
    t = toc(t0);
    save(flagFullname, 't')    
end

end

