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


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tifFilename', @(x) ischar(x) || iscell(x));
ip.addRequired('zarrFilename', @ischar);
ip.addOptional('frame', [], @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('blockSize', [500, 500, 250], @isnumeric);
ip.addParameter('expand2dDim', true, @islogical); % expand the z dimension for 2d data
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('InputBbox', [], @(x) isnumeric(x));
ip.addParameter('CropToSize', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('readWholeTiff', true, @islogical);
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x));
ip.addParameter('uuid', '', @ischar);

ip.parse(tifFilename, zarrFilename, frame, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
blockSize = pr.blockSize;
expand2dDim = pr.expand2dDim;
flipZstack = pr.flipZstack;
resample = pr.resample;
InputBbox = pr.InputBbox;
CropToSize = pr.CropToSize;
readWholeTiff = pr.readWholeTiff;
usrFcn = pr.usrFcn;
uuid = pr.uuid;

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

if isempty(tifFilename)
    fprintf('Convert frame to Zarr for %s ...\n', zarrFilename);    
else
    if ischar(tifFilename)
        fprintf('Convert Tiff to Zarr for %s ...\n', tifFilename);
    else
        fprintf('Convert Tiff to Zarr for %s ...\n', tifFilename{1});        
    end
end

% if InputBbox is nonempty, read whole tiff
if ~isempty(InputBbox)
    readWholeTiff = true;
end

% load tiff file as image block file
% Direct conversion should be faster. 
if ~isempty(frame)
    if ismatrix(frame)
        blockSize = blockSize(1 : 2);
    end
    blockSize = min(blockSize, size(frame));    
    bim = blockedImage(frame, "BlockSize", blockSize);
else
    if ischar(tifFilename) || numel(tifFilename) == 1
        if iscell(tifFilename)
            tifFilename = tifFilename{1};
        end
        if readWholeTiff
            % I = readtiff(tifFilename);
            try 
                I = parallelReadTiff(tifFilename);
            catch 
                I = readtiff_parallel(tifFilename);
            end
            if ~isempty(InputBbox)
                I = I(InputBbox(1) : InputBbox(4), InputBbox(2) : InputBbox(5), InputBbox(3) : InputBbox(6));
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
            bim = blockedImage(I, "BlockSize", blockSize(1 : ndims(I)));
            clear I;
        else
            try
                bim = blockedImage(tifFilename, "Adapter", MPageTiffAdapter);
            catch ME
                disp(ME)
                fprintf('Use BioFormats to read the tif file ...\n');
                im = bfOpen3DVolume(tifFilename);
                bim = blockedImage(im{1,1}{1,1});
                clear im;
            end
        end
    else
        % with partial files
        nF = numel(tifFilename);
        I = cell(nF, 1);
        for i = 1 : nF
            % I{i} = readtiff(tifFilename{i});
            try 
                I{i} = parallelReadTiff(tifFilename{i});
            catch
                I{i} = readtiff_parallel(tifFilename{i});
            end
        end
        I = cat(3, I{:});
        if ~isempty(InputBbox)
            I = I(InputBbox(1) : InputBbox(4), InputBbox(2) : InputBbox(5), InputBbox(3) : InputBbox(6));
        end
        if ismatrix(I)
            blockSize = blockSize(1 : 2);
        end        
        blockSize = min(blockSize, size(I));        
        bim = blockedImage(I, "BlockSize", blockSize);
        clear I;
    end
    if ~isempty(gcp('nocreate'))
        parfevalOnAll(@clearvars, 0)
    end
end

sz = bim.Size;
dtype = bim.ClassUnderlying;
if any(sz == 0)
    fprintf('The input image is empty, skip it!\n');
    return;
end    

if flipZstack
    bim = apply(bim, @(bs) flip(bs.Data, 3), 'BlockSize', sz);
end    

% centerred crop
if ~isempty(CropToSize)
    halfCrop = floor((bim.Size - CropToSize) / 2);
    bbox = [halfCrop + 1, halfCrop +  CropToSize];
    if isa(bim.Adapter, 'images.blocked.InMemory')
        bim = gather(bim);
        bim = bim(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
        bim = blockedImage(bim, 'BlockSize', blockSize);
    else
        bim = bim.crop(bbox(1 : 3), bbox(4 : 6));
    end
end

if ~isempty(resample) && ~all(resample == 1)
    rs = resample(:)';
    % complete rs to 3d in case it is not
    rs = [ones(1, 4 - numel(rs)) * rs(1), rs(2:end)];    
    outSize = round(bim.Size ./ rs);
    bim = apply(bim, @(bs) imresize3(bs.Data, outSize), 'BlockSize', sz);
end
sz = bim.Size;

if ~isempty(usrFcn)
    if ischar(usrFcn)
        usrFcn = str2func(usrFcn);
    end
    bim = apply(bim, @(bs) usrFcn(bs.Data), 'BlockSize', sz);
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
% setBlock(bo, [1, 1, 1], blockData);
tmpFilename = [zarrFilename '_' uuid];
% write(bim, tmpFilename, "BlockSize", blockSize, "Adapter", ZarrAdapter);
if ~exist(tmpFilename, 'dir')
    init_val = zeros(1, dtype);
    nv_bim = blockedImage(tmpFilename, sz, blockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');
    nv_bim.Adapter.close();
end

% write zarr
writezarr(cast(gather(bim), dtype), tmpFilename, 'blockSize', blockSize);

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
fprintf('Done!\n');


end

