function [] = tiffToZarr(tifFilename, zarrFilename, frame, varargin)
% use image block function to convert tiff to zarr
% 
% Author: Xiongtao Ruan (10/02/2020)
% 
% xruan (10/11/2020): add function handle for processing before saving to zarr
% xruan (10/24/2020): add support for variable in memory
% xruan (11/10/2020): check whether the image is empty, and return error for empty image
% xruan (11/18/2020): use bioformat to open tif file when fail to read (python saved file)
% xruan (07/26/2021): add support for flipped tile
% xruan (07/27/2021): add support for resampling
% xruan (08/18/2021): change to write to all zarr blocks for a block of tiff z slices
% xruan (09/23/2021): add support for including partial files 
% xruan (10/13/2021): add support for cropping data after processing
% xruan (02/02/2022): add support for cropping the input data
% xruan (07/05/2022): change zarr writer to writezarr.m
% xruan (08/25/2022): change CropToSize to tileOutBbox (more generic)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tifFilename', @(x) ischar(x) || iscell(x));
ip.addRequired('zarrFilename', @ischar);
ip.addOptional('frame', [], @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('blockSize', [500, 500, 250], @isnumeric);
ip.addParameter('shardSize', [], @isnumeric);
ip.addParameter('expand2dDim', true, @islogical); % expand the z dimension for 2d data
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('resampleFactor', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x));
ip.addParameter('tileOutBbox', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('readWholeTiff', true, @islogical);
ip.addParameter('compressor', 'zstd', @ischar);
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x) || isstring(x));
ip.addParameter('uuid', '', @ischar);

ip.parse(tifFilename, zarrFilename, frame, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
blockSize = pr.blockSize;
shardSize = pr.shardSize;
expand2dDim = pr.expand2dDim;
flipZstack = pr.flipZstack;
resampleFactor = pr.resampleFactor;
inputBbox = pr.inputBbox;
tileOutBbox = pr.tileOutBbox;
readWholeTiff = pr.readWholeTiff;
compressor = pr.compressor;
usrFcn = pr.usrFcn;
uuid = pr.uuid;

% t0 = tic;

% remove the last slash
zarrFilename = strip(zarrFilename, 'right', '/');
if exist(zarrFilename, 'dir')
    if ~Overwrite
        fprintf('Zarr result %s already exists, skip it!\n', zarrFilename);
        return;
    else
        rmdir(zarrFilename, 's');
    end
end

if isempty(uuid)
    uuid = get_uuid();
end
% with a full random uuid made here (get_uuid makes shorter ones on Windows
% and as a fallback) no other process writes the temporary result, so its
% chunks are written without cpp-zarr's per-chunk temporary files
privateTmp = isempty(pr.uuid) && numel(uuid) >= 32;

if isempty(tifFilename)
    fprintf('Convert frame to Zarr for %s ...\n', zarrFilename);    
else
    if ischar(tifFilename)
        fprintf('Convert Tiff to Zarr for %s ...\n', tifFilename);
    else
        fprintf('Convert Tiff to Zarr for %s ...\n', tifFilename{1});        
    end
end

applyUserFcn = ~(isstring(usrFcn) || isempty(usrFcn)) || (isstring(usrFcn) && ~isempty(usrFcn{1}));

% if flipZstack is true, or resampleFactor is nonempty, or apply user
% function, set read whole tiff as true.
if flipZstack || ~isempty(resampleFactor) || applyUserFcn
    readWholeTiff = true;
end

% load the image; it is processed as an array, and as a blocked image only
% when the tiff is read by regions (readWholeTiff false)
I = [];
bim = [];
if ~isempty(frame)
    if ismatrix(frame)
        blockSize = blockSize(1 : 2);
    end
    blockSize = min(blockSize, size(frame));    
    I = frame;
else
    if ischar(tifFilename) || numel(tifFilename) == 1
        if iscell(tifFilename)
            tifFilename = tifFilename{1};
        end
        if readWholeTiff
            if ~isempty(inputBbox)
                I = readtiff(tifFilename, range=[inputBbox(3), inputBbox(6)]);
                I = crop3d(I, [inputBbox(1 : 2), 1, inputBbox(4 : 5), size(I, 3)]);
            else
                I = readtiff(tifFilename);
            end
            sz = size(I);
            if ismatrix(I)
                if expand2dDim
                    blockSize(3) = 1;
                    sz(3) = 1;
                else
                    blockSize =  blockSize(1 : 2);
                end
            end
            blockSize = min(blockSize, sz);
        else
            try
                bim = blockedImage(tifFilename, "Adapter", MPageTiffAdapter);
            catch ME
                disp(ME)
                fprintf('Use BioFormats to read the tiff file ...\n');
                im = bfOpen3DVolume(tifFilename);
                I = im{1,1}{1,1};
                clear im;
            end
            if ~isempty(inputBbox)
                if isempty(bim)
                    I = crop3d(I, inputBbox);
                else
                    bim = bim.crop(inputBbox(1 : 3), inputBbox(4 : 6));
                end
            end
        end
    else
        % with partial files
        nF = numel(tifFilename);
        I = cell(nF, 1);
        for i = 1 : nF
            I{i} = readtiff(tifFilename{i});
        end
        I = cat(3, I{:});
        if ~isempty(inputBbox)
            I = crop3d(I, inputBbox);
        end
        if ismatrix(I)
            blockSize = blockSize(1 : 2);
        end        
        blockSize = min(blockSize, size(I));        
    end
    if ~isempty(gcp('nocreate'))
        parfevalOnAll(@clearvars, 0)
    end
end

if isempty(bim)
    sz = size(I);
    dtype = class(I);
else
    sz = bim.Size;
    dtype = bim.ClassUnderlying;
end
if any(sz == 0)
    fprintf('The input image is empty, skip it!\n');
    return;
end    

if flipZstack
    I = flip(I, 3);
end    

% bounding box crop for output
if ~isempty(tileOutBbox)
    % halfCrop = floor((bim.Size - CropToSize) / 2);
    % bbox = [halfCrop + 1, halfCrop +  CropToSize];
    bbox = tileOutBbox;
    if isempty(bim)
        I = crop3d(I, bbox);
    else
        bim = bim.crop(bbox(1 : 3), bbox(4 : 6));
    end
end

if ~isempty(resampleFactor) && ~all(resampleFactor == 1)
    rs = resampleFactor(:)';
    % complete rs to 3d in case it is not
    rs = [ones(1, 4 - numel(rs)) * rs(1), rs(2:end)];    
    outSize = round(size(I) ./ rs);
    I = imresize3(I, outSize);
end
if isempty(bim)
    sz = size(I);
else
    sz = bim.Size;
end

if applyUserFcn
    disp(usrFcn)
    if ischar(usrFcn) || isstring(usrFcn)
        usrFcn = str2func(usrFcn);
    end
    I = usrFcn(I);
end
if numel(sz) == 2
    if expand2dDim 
        sz(3) = 1;
        blockSize(3) = 1;
    else
        if sz(1) <= sz(2)
            blockSize = min(sz, [sz(1), round(prod(blockSize) / sz(1))]);
        else
            blockSize = min(sz, [round(prod(blockSize) / sz(2)), sz(2)]);
        end  
    end
else
    blockSize = min(sz, blockSize);
end
tmpFilename = [zarrFilename '_' uuid];
if ~exist(tmpFilename, 'dir')
    dimSeparator = '.';
    if prod(ceil(sz ./ blockSize)) > 10000
        dimSeparator = '/';
    end
    createzarr(tmpFilename, dataSize=sz, blockSize=blockSize, shardSize=shardSize, ...
        dtype=dtype, expand2dDim=expand2dDim, compressor=compressor, dimSeparator=dimSeparator);
end

% write zarr
if isempty(bim)
    writezarr(cast(I, dtype), tmpFilename, 'blockSize', blockSize, 'uuid', ~privateTmp);
else
    for i = 1 : ceil(sz(3) / blockSize(3))
        zs = (i - 1) * blockSize(3) + 1;
        zt = min(i * blockSize(3), sz(3));
        bbox = [1, 1, zs, sz(1), sz(2), zt];
        if ~isempty(inputBbox) || ~isempty(tileOutBbox)
            bbox_in = bbox + [floor(bim.WorldStart), floor(bim.WorldStart)];
        else
            bbox_in = bbox;
        end
        im_i = bim.Adapter.getIORegion(bbox_in(1 : 3), bbox_in(4 : 6));
        writezarr(im_i, tmpFilename, bbox=bbox, uuid=~privateTmp);
    end
end

% mv tmp result folder to output folder
if exist(zarrFilename, 'dir')
    if ~Overwrite
        fprintf('Zarr result %s already exists, skip it!\n', zarrFilename);
        rmdir(tmpFilename, 's');
        return;
    else
        rmdir(zarrFilename, 's');
    end
end
movefile(tmpFilename, zarrFilename);
% toc(t0);
% fprintf('Done!\n');

end

