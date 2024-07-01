function [] = stitch_process_zarr_tile(inputFilename, zarrFilename, frame, varargin)
% process zarr tile: resample, crop, flip or applying user defined function handle
% 
% copied from tiffToZarr.m
%
% Author: Xiongtao Ruan (04/29/2023)



ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFilename', @(x) ischar(x) || iscell(x));
ip.addRequired('zarrFilename', @ischar);
ip.addOptional('frame', [], @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('blockSize', [500, 500, 250], @isnumeric);
ip.addParameter('shardSize', [], @isnumeric); 
ip.addParameter('expand2dDim', true, @islogical); % expand the z dimension for 2d data
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('inputBbox', [], @(x) isnumeric(x));
ip.addParameter('tileOutBbox', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('compressor', 'lz4', @ischar);
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x) || isstring(x));
ip.addParameter('uuid', '', @ischar);

ip.parse(inputFilename, zarrFilename, frame, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
blockSize = pr.blockSize;
shardSize = pr.shardSize;
expand2dDim = pr.expand2dDim;
flipZstack = pr.flipZstack;
resample = pr.resample;
inputBbox = pr.inputBbox;
tileOutBbox = pr.tileOutBbox;
compressor = pr.compressor;
usrFcn = pr.usrFcn;
uuid = pr.uuid;

t0 = tic;

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

if isempty(inputFilename)
    fprintf('Process zarr tile for %s ...\n', zarrFilename);    
else
    if ischar(inputFilename)
        fprintf('Process zarr tile for %s ...\n', inputFilename);
    else
        fprintf('Process zarr tile for %s ...\n', inputFilename{1});        
    end
end

% Direct conversion should be faster. 
if ~isempty(frame)
    if ismatrix(frame)
        blockSize = blockSize(1 : 2);
    end
else
    frame = readzarr(inputFilename, inputBbox=inputBbox);

    sz = size(frame);
    if ismatrix(frame)
        if expand2dDim
            blockSize(3) = 1;
            sz(3) = 1;
        else
            blockSize =  blockSize(1 : 2);                    
        end
    end
    blockSize = min(blockSize, sz);
end

dtype = class(frame);
if any(sz == 0)
    fprintf('The input image is empty, skip it!\n');
    return;
end    

if flipZstack
    frame = flip(frame, 3);
end    

% bounding box crop for output
frame = crop3d(frame, tileOutBbox);

if ~isempty(resample) && ~all(resample == 1)
    rs = resample(:)';
    % complete rs to 3d in case it is not
    rs = [ones(1, 4 - numel(rs)) * rs(1), rs(2:end)];    
    outSize = round(size(frame) ./ rs);
    frame = imresize3(frame, outSize);
end

if ~(isstring(usrFcn) || isempty(usrFcn)) || (isstring(usrFcn) && ~isempty(usrFcn{1}))
    disp(usrFcn)
    if ischar(usrFcn) || isstring(usrFcn)
        usrFcn = str2func(usrFcn);
    end
    frame = usrFcn(frame);
end
sz = size(frame);
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
writezarr(cast(frame, dtype), tmpFilename, 'blockSize', blockSize);

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
toc(t0);
fprintf('Done!\n');

end

