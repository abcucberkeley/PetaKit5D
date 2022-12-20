function [] = RLdecon_large_in_memory(frameFullpath, psfFullpath, deconFullpath, xyPixelSize, dz, varargin)
% RL decon for big data that can be loaded to the memory. 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('deconFullpath', @(x) ischar(x));
ip.addRequired('psfFullpath', @(x) ischar(x));
ip.addRequired('xyPixelSize', @isnumeric);
ip.addRequired('dz', @isnumeric);
ip.addParameter('Save16bit', false , @islogical);
ip.addParameter('Rotate', false , @islogical);
ip.addParameter('Deskew', false , @islogical);
ip.addParameter('SkewAngle', -32.45 , @isnumeric);
ip.addParameter('flipZstack', false, @islogical); 
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('EdgeSoften', 5, @isnumeric); % # ofxy px to soften
ip.addParameter('zEdgeSoften', 2, @isnumeric); % # ofxy px to soften
ip.addParameter('Crop', [], @isnumeric); % requires lower and higher values for cropping
ip.addParameter('zFlip', false, @islogical);
ip.addParameter('GenMaxZproj', [0,0,1] , @isnumeric);
ip.addParameter('dzPSF', 0.1 , @isnumeric); %in um
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('EdgeErosion', 8 , @isnumeric); % erode edges for certain size.
ip.addParameter('ErodeMaskfile', '', @ischar); % erode edges file
ip.addParameter('SaveMaskfile', false, @islogical); % save mask file for common eroded mask
ip.addParameter('scaleFactor', 1e8 , @isnumeric); % scale factor for data
ip.addParameter('deconMaskFns', {} , @iscell); % Full paths of 2D mask zarr files, in xy, xz, yz order
ip.addParameter('saveZarr', false, @islogical); % save as zarr
% ip.addParameter('DoNotAdjustResForFFT', true , @islogical); % not crop chunks for deconvolution
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @isnumeric);
ip.addParameter('OTFCumThresh', 0.9, @isnumeric); % OTF cumutative sum threshold
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x)); % decon in skewed space
ip.addParameter('fixIter', true, @islogical); % 
ip.addParameter('errThresh', [], @isnumeric); % error threshold for simplified code
ip.addParameter('BatchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('Overlap', 200, @isnumeric); % block overlap
ip.addParameter('CPUMaxMem', 500, @isnumeric); % CPU Memory in Gb
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('psfGen', true, @islogical); % psf generation
ip.addParameter('useGPU', true, @islogical); % use gpu if it is available

ip.parse(frameFullpath, psfFullpath, deconFullpath, xyPixelSize, dz, varargin{:});

pr = ip.Results;
% Overwrite = pr.Overwrite;
Save16bit = pr.Save16bit;
Deskew = pr.Deskew;
Rotate = pr.Rotate;
SkewAngle = pr.SkewAngle;
flipZstack = pr.flipZstack;
dzPSF = pr.dzPSF;
DeconIter = pr.DeconIter;
% scaleFactor = pr.scaleFactor;
deconMaskFns = pr.deconMaskFns;
saveZarr = pr.saveZarr;
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
skewed = pr.skewed;
fixIter = pr.fixIter;
errThresh = pr.errThresh;
BatchSize = pr.BatchSize;
psfGen = pr.psfGen;
useGPU = pr.useGPU;
EdgeErosion = pr.EdgeErosion;
ErodeMaskfile = pr.ErodeMaskfile;
if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
    EdgeErosion = 0; % set EdgeErosion length as 0.
end
SaveMaskfile = pr.SaveMaskfile;

% check if background information available, if not, estimate background
% info. Currently use 99. 
Background = pr.Background;
if isempty(Background)
    Background = 99;
end
uuid = pr.uuid;

% simplified version related options
tic
[dataPath, fsname, ext] = fileparts(frameFullpath);
[deconPath, fsn] = fileparts(deconFullpath);
fprintf('Start Large-file RL Decon for %s...\n', fsname);

% generate psf and use the cropped psf to decide the overlap region
pp = readtiff(psfFullpath);
if psfGen
    medFactor = 1.5;
    PSFGenMethod = 'masked';
    if Deskew
        dz_data_ratio = sind(SkewAngle);
    else
        dz_data_ratio = 1;
    end
    psf = psf_gen_new(pp, dzPSF, dz*dz_data_ratio, medFactor, PSFGenMethod);
    py = find(squeeze(sum(psf, [2, 3])));
    px = find(squeeze(sum(psf, [1, 3])));
    pz = find(squeeze(sum(psf, [1, 2])));
    cropSz = [min(py(1) - 1, size(psf, 1) - py(end)), min(px(1) - 1, size(psf, 2) - px(end)), min(pz(1) - 1, size(psf, 3) - pz(end))] - 1;
    cropSz = max(0, cropSz);
    bbox = [cropSz + 1, size(psf) - cropSz];
    psf = psf(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
    
    % write generated psf to disk in case it will be deleted in the c
    psfgen_folder = sprintf('%s/psfgen/', deconPath);
    mkdir(psfgen_folder);
    [~, psf_fsn] = fileparts(psfFullpath);
    psfgen_filename = sprintf('%s/%s.tif', psfgen_folder, psf_fsn);
    if ~exist(psfgen_filename, 'file')
        writetiff(psf, psfgen_filename);
    end
else
    psf = pp;
end

tic
fprintf(['reading ' fsname '...\n'])
switch ext
    case {'.tif', '.tiff'}
        im = readtiff(frameFullpath);
    case '.zarr'
        im = readzarr(frameFullpath);
end
if flipZstack
    im = flip(im, 3);
end
toc

if EdgeErosion > 0
    if SaveMaskfile
        maskPath = [deconPath, '/', 'Masks'];
        if ~exist(maskPath, 'dir')
            mkdir(maskPath);
        end
        maskFullPath = sprintf('%s/%s_eroded.tif', maskPath, fsname);
        maskTmpPath = sprintf('%s/%s_eroded_%s.tif', maskPath, fsname, uuid);
    end
    if ~(SaveMaskfile && exist(maskFullPath, 'file'))
        im_bw_erode = im > 0;
        im_bw_erode([1, end], :, :) = false;
        im_bw_erode(:, [1, end], :) = false;
        im_bw_erode(:, :, [1, end]) = false;
        im_bw_erode = bwdist(~im_bw_erode) > EdgeErosion - 1;
    else
        im_bw_erode = readtiff(maskFullPath) > 0;
    end

    % save mask file as common one for other time points/channels
    if SaveMaskfile
        writetiff(uint8(im_bw_erode), maskTmpPath);
        movefile(maskTmpPath, maskFullPath);
    end
end

% calculate number of chunks to break the image file
imSize = size(im);
% dtype = class(im);
if Save16bit
    dtype = 'uint16';
else
    dtype = 'single';
end

SameBatchSize = ~true;
BorderSize = round((size(psf) + 10) / 2);
BlockSize = BatchSize;
[BatchBBoxes, RegionBBoxes] = XR_zarrChunkCoordinatesExtraction(imSize, 'BatchSize', BatchSize, ...
    'BlockSize', BlockSize, 'SameBatchSize', SameBatchSize, 'BorderSize', BorderSize);
% scaleFactor = prod(min(imSize, BatchSize));
scaleFactor = 1.0;

t0 = tic;
fprintf('Deconvolving chunks...\n')
batchInds = 1 : size(BatchBBoxes, 1);
finds = [1, 2; 2, 3; 1, 3];
done_flag = false(numel(batchInds), 1);
% imout = 0 * im;
imout = zeros(imSize, dtype);
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
            imout(obStart(1) : obEnd(1), obStart(2) : obEnd(2), obStart(3) : obEnd(3)) = zeros(obEnd - obStart + 1, dtype);
            done_flag(i) = true;
            toc;
            continue;
        end
    end

    % load the region in input 
    try
        in_batch = crop3d_mex(im, [ibStart, ibEnd]);
    catch ME
        disp(ME);
        in_batch = im(ibStart(1) : ibEnd(1), ibStart(2) : ibEnd(2), ibStart(3) : ibEnd(3));
    end

    % deconvolution
    %{
    frameFullpath = '';
    deconTmpPath = '';
    Deskew = false;
    Rotate = false;
    Save16bit = ~false;
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
        'rawdata', in_batch, 'scaleFactor', scaleFactor, 'useGPU', useGPU, 'psfGen', psfGen);
    %}

    inputFn = '';
    outputFn = deconFullpath;
    DSRCombined = true;
    Reverse = true;
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
        'skewed', skewed, 'wienerAlpha', wienerAlpha, 'fixIter', fixIter, 'scaleFactor', scaleFactor, ...
        'deconBbox', deconBbox, 'useGPU', useGPU, 'psfGen', psfGen, 'debug', debug, ...
        'save3Dstack', save3Dstack, 'mipAxis', mipAxis);    
    %{
    try
        out_batch = crop3d_mex(out_batch, [baStart, baEnd]);
    catch ME
        disp(ME);
        out_batch = out_batch(baStart(1) : baEnd(1), baStart(2) : baEnd(2), baStart(3) : baEnd(3));
    end
    %}
    
    try 
        indexing3d_mex(imout, [obStart, obEnd], out_batch);
    catch ME
        disp(ME);
        imout(obStart(1) : obEnd(1), obStart(2) : obEnd(2), obStart(3) : obEnd(3)) = out_batch;
    end

    done_flag(i) = true;

    toc;
end

im = imout;
clear imout;

if ~all(done_flag)
    error('Deconvolution of chunk %s is missing!', mat2str(find(~done_flag)));
end

tic
fprintf('Saving combined deconvolved chunks...\n')
if EdgeErosion > 0
    fprintf('Erode edges of deconvolved data w.r.t. raw data...\n');        
    im = im .* cast(im_bw_erode, class(im));
end

if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
    fprintf('Erode edges of deconvolved data using a predefined mask...\n');                
    im_bw_erode = readtiff(ErodeMaskfile);   
    im = im .* cast(im_bw_erode, class(im));
end

if pr.Save16bit
    im = uint16(im);
else
    im = single(im);
end

if saveZarr
    deconTmpPath = sprintf('%s_%s.zarr', deconFullpath(1 : end - 5), uuid);
    writezarr(im, deconTmpPath);
else
    deconTmpPath = sprintf('%s_%s.tif', deconFullpath(1 : end - 4), uuid);
    writetiff(im, deconTmpPath);
end
movefile(deconTmpPath, deconFullpath);

tmp_xy = max(im,[],3);
deconMIPPath = sprintf('%s/MIPs/', deconPath);
mkdir(deconMIPPath);
deconMIPFullPath = sprintf('%s%s_MIP_z.tif', deconMIPPath, fsname);
writetiff(tmp_xy, deconMIPFullPath);
toc(t0);

end

