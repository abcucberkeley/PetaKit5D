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
ip.addRequired('xyPixelSize', @isnumeric); %in um
ip.addRequired('dz', @isnumeric); %in um
% ip.addParameter('BlockSize', [], @isnumeric);
ip.addParameter('save16bit', true , @islogical);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('SkewAngle', -32.45 , @isnumeric);
ip.addParameter('flipZstack', false, @islogical); 
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('dzPSF', 0.1 , @isnumeric); %in um
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('dampFactor', 1, @isnumeric); % damp factor for decon result
ip.addParameter('scaleFactor', 1.0, @isnumeric); % scale factor for decon result
ip.addParameter('deconOffset', 0, @isnumeric); % offset for decon result
ip.addParameter('EdgeErosion', 0, @isnumeric); % edge erosion for decon result
ip.addParameter('maskFullpaths', {} , @iscell); % Full paths of 2D mask zarr files, in xy, xz, yz order
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @isnumeric);
ip.addParameter('OTFCumThresh', 0.9, @isnumeric); % OTF cumutative sum threshold
ip.addParameter('hannWinBounds', [0.8, 1.0], @isnumeric); % apodization range for distance matrix
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x)); % decon in skewed space
ip.addParameter('fixIter', false, @islogical); % CPU Memory in Gb
ip.addParameter('useGPU', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('psfGen', true, @islogical); % psf generation

ip.parse(batchInds, zarrFullpath, psfFullpath, deconFullpath, flagFullname, ...
    BatchBBoxes, RegionBBoxes, xyPixelSize, dz, varargin{:});

pr = ip.Results;
save16bit = pr.save16bit;
Overwrite = pr.Overwrite;
SkewAngle = pr.SkewAngle;
flipZstack = pr.flipZstack;
Background = pr.Background;
dzPSF = pr.dzPSF;
DeconIter = pr.DeconIter;
dampFactor = pr.dampFactor;
scaleFactor = pr.scaleFactor;
deconOffset = pr.deconOffset;
EdgeErosion = pr.EdgeErosion;
maskFullpaths = pr.maskFullpaths;
RLMethod = pr.RLMethod;
skewed = pr.skewed;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
hannWinBounds = pr.hannWinBounds;
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

if ~exist(deconFullpath, 'dir')
    error('The output zarr file %s doesnot exist!', deconFullpath);
end

zarrFile = true;
dtype = getImageDataType(deconFullpath, zarrFile);
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
    if ~(isempty(maskFullpaths) || isempty(maskFullpaths{1}))
        skipDecon = false;
        for f = 1 : 3
            im_f = readzarr(maskFullpaths{f}, 'inputBbox', [obStart(finds(f, :)), 1, obEnd(finds(f, :)), 1]);
            if ~any(im_f(:))
                skipDecon = true;
                break;
            end
        end
        
        if skipDecon
            % nv_bim.Adapter.setRegion(obStart, obEnd, zeros(obEnd - obStart + 1, dtype));
            writezarr(zeros(obEnd - obStart + 1, dtype), deconFullpath, 'bbox', [obStart, obEnd], create=false);
            done_flag(i) = true;
            toc;
            continue;
        end
    end

    % load the region in input 
    in_batch = readzarr(zarrFullpath, 'inputBbox', [ibStart, ibEnd]);

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
        'rawdata', in_batch, 'save16bit', save16bit, 'SkewAngle', SkewAngle, ...
        'Deskew', Deskew, 'Rotate', Rotate, 'DSRCombined', DSRCombined, 'Reverse', Reverse, ...
        'Background', Background, 'DeconIter', DeconIter, 'RLMethod', RLMethod, ...
        'skewed', skewed, 'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, ...
        'hannWinBounds', hannWinBounds, 'fixIter', fixIter, 'dampFactor', dampFactor, ...
        'scaleFactor', scaleFactor, 'deconOffset', deconOffset, 'EdgeErosion', EdgeErosion, ...
        'deconBbox', deconBbox, 'useGPU', useGPU, 'psfGen', psfGen, 'debug', debug, ...
        'save3Dstack', save3Dstack, 'mipAxis', mipAxis);
    
    clear in_batch;

    writezarr(out_batch, deconFullpath, 'bbox', [obStart, obEnd], create=false);
    done_flag(i) = true;
    toc;

    % in case the job is finished by other jobs, check the flag file every
    % five loops if not using GPU
    if ~useGPU && rem(i, 5) == 0 && exist(flagFullname, 'file')
        break;
    end
end

if all(done_flag)
    % fclose(fopen(flagFullname, 'w'));
    t = toc(t0);
    save(flagFullname, 't')    
end

end

