function [] = XR_resampleSingleZarr(zarrFullpath, dsFullpath, dsFactor, varargin)
% 
% Author: Xiongtao Ruan (12/14/2020)

% xruan (11/09/2021): enable arbitray blockSize
% xruan ( 06/22/2023): add support for bbox for input


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath'); 
ip.addRequired('dsFullpath'); 
ip.addRequired('dsFactors'); 
ip.addParameter('bbox', [], @isnumeric); % bbox for input
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @isnumeric); % size to process in one batch 
ip.addParameter('zarrSubSize', [20, 20, 20], @isnumeric);
ip.addParameter('BorderSize', [5, 5, 5], @isnumeric); % padded boarder for each batch
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear', 'nearest'})));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(zarrFullpath, dsFullpath, dsFactor, varargin{:});

pr = ip.Results;
bbox = pr.bbox;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
zarrSubSize = pr.zarrSubSize;
BorderSize = pr.BorderSize;
Interp = pr.Interp;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if isempty(uuid)
    uuid = get_uuid();
end

if numel(dsFactor) == 1
    dsFactor = ones(1, 3) * dsFactor;
elseif numel(dsFactor) == 2
    dsFactor = [dsFactor(1), dsFactor(1), dsFactor(2)];
end

% check if the group folder exist

if exist(dsFullpath, 'dir')
    disp('The output result exists, skip it!');
    return;
end

[dsPath, dsfsname] = fileparts(dsFullpath);
if strcmp(dsFullpath(end), '/')
    dsTmppath = [dsFullpath(1 : end - 1), '_', uuid];
else
    dsTmppath = [dsFullpath, '_', uuid];
end    

if ~exist(zarrFullpath, 'dir')
    error('The input file %s does not exist, pleast double check!', zarrFullpath);
end
try
    bim = blockedImage(zarrFullpath, "Adapter", CZarrAdapter);
catch ME
    disp(ME);
    bim = blockedImage(zarrFullpath, "Adapter", ZarrAdapter);
end    
dtype = bim.ClassUnderlying;
sz = bim.Size;

inSize = sz;
if ~isempty(bbox)
    inSize = bbox(4 : 6) - bbox(1 : 3) + 1;
end

ds_size = round(inSize ./ [dsFactor(1), dsFactor(2), dsFactor(3)]);

% framework
blockSize = min(ds_size, blockSize);
batchSize = min(ds_size, batchSize);
batchSize = ceil(batchSize ./ blockSize) .* blockSize;
bSubSz = ceil(ds_size ./ batchSize);
numBatch = prod(bSubSz);

% process for each block based on all BlockInfo use distributed computing
fprintf('Process blocks for resampled data...\n')

zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', dsPath, dsfsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end

% initialize zarr file
if exist(dsTmppath, 'dir')
    bim = blockedImage(dsTmppath, 'Adapter', CZarrAdapter);
    if any(bim.BlockSize ~= blockSize) || any(bim.Size ~= ds_size)
        rmdir(dsTmppath, 's');
        rmdir(zarrFlagPath, 's');
        mkdir(zarrFlagPath);
    end
end
if ~exist(dsTmppath, 'dir')
    createzarr(dsTmppath, dataSize=ds_size, blockSize=blockSize, dtype=dtype, zarrSubSize=zarrSubSize);
end

taskSize = max(5, ceil(numBatch / 5000)); % the number of batches a job should process
numTasks = ceil(numBatch / taskSize);

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    
    zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;
    
    funcStrs{i} = sprintf(['resampleZarrBlock([%s],''%s'',''%s'',''%s'',[%s],', ...
        '''Interp'',''%s'',''bbox'',[%s],''BorderSize'',[%s],''batchSize'',[%s],', ...
        '''blockSize'',%s)'], strrep(num2str(batchInds, '%d,'), ' ', ''), zarrFullpath, ...
        dsTmppath, zarrFlagFullpath, strrep(num2str(dsFactor(:)', '%d,'), ' ', ''), Interp, ...
        strrep(num2str(bbox(:)', '%d,'), ' ', ''), strrep(num2str(BorderSize(:)', '%d,'), ' ', ''), ...
        strrep(num2str(batchSize(:)', '%d,'), ' ', ''), strrep(mat2str(blockSize(:)'), ' ', ','));
end

inputFullpaths = repmat({zarrFullpath}, numTasks, 1);

dtype = getImageDataType(zarrFullpath);
byteNum = dataTypeToByteNumber(dtype);
memAllocate = prod(batchSize) * byteNum / 2^30 * (prod(dsFactor) + 1) * 3;
is_done_flag= generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
    funcStrs, cpusPerTask=cpusPerTask, memAllocate=memAllocate, parseCluster=parseCluster, ...
    mccMode=mccMode, ConfigFile=ConfigFile);

if ~all(is_done_flag)
    is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, cpusPerTask=cpusPerTask, memAllocate=memAllocate, parseCluster=parseCluster, ...
        mccMode=mccMode, ConfigFile=ConfigFile);
end

if ~all(is_done_flag)
    error('%d / %d tasks failed!', sum(~is_done_flag), numel(is_done_flag));
end

if exist(dsFullpath, 'dir')
    rmdir(dsFullpath, 's');
end
movefile(dsTmppath, dsFullpath);
rmdir(zarrFlagPath, 's');
fprintf('Done!\n')

end

