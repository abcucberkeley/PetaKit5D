function [] = XR_resampleSingleZarr(zarrFullpath, rsFullpath, resampleFactor, varargin)
% 
% Author: Xiongtao Ruan (12/14/2020)

% xruan (11/09/2021): enable arbitray blockSize
% xruan ( 06/22/2023): add support for inputBbox for input


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath', @ischar); 
ip.addRequired('rsFullpath', @ischar); 
ip.addRequired('resampleFactor', @isnumeric); 
ip.addParameter('inputBbox', [], @isnumeric); % inputBbox for input
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @isnumeric); % size to process in one batch 
ip.addParameter('borderSize', [5, 5, 5], @isnumeric); % padded boarder for each batch
ip.addParameter('interpMethod', 'linear', @(x) ischar(x) && any(strcmpi(x, {'cubic', 'linear', 'nearest'})));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('cpusPerTask', 1, @isscalar);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(zarrFullpath, rsFullpath, resampleFactor, varargin{:});

pr = ip.Results;
inputBbox = pr.inputBbox;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
borderSize = pr.borderSize;
interpMethod = pr.interpMethod;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if isempty(uuid)
    uuid = get_uuid();
end

if numel(resampleFactor) == 1
    resampleFactor = ones(1, 3) * resampleFactor;
elseif numel(resampleFactor) == 2
    resampleFactor = [resampleFactor(1), resampleFactor(1), resampleFactor(2)];
end

% check if the group folder exist

if exist(rsFullpath, 'dir')
    disp('The output result exists, skip it!');
    return;
end

[rsPath, rsfsname] = fileparts(rsFullpath);
if strcmp(rsFullpath(end), '/')
    rsTmpPath = [rsFullpath(1 : end - 1), '_', uuid];
else
    rsTmpPath = [rsFullpath, '_', uuid];
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
if ~isempty(inputBbox)
    inSize = inputBbox(4 : 6) - inputBbox(1 : 3) + 1;
end

rs_size = round(inSize ./ [resampleFactor(1), resampleFactor(2), resampleFactor(3)]);

% framework
blockSize = min(rs_size, blockSize);
batchSize = min(rs_size, batchSize);
batchSize = ceil(batchSize ./ blockSize) .* blockSize;
bSubSz = ceil(rs_size ./ batchSize);
numBatch = prod(bSubSz);

% process for each block based on all BlockInfo use distributed computing
fprintf('Process blocks for resampled data...\n')

zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', rsPath, rsfsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end

% initialize zarr file
if exist(rsTmpPath, 'dir')
    bim = blockedImage(rsTmpPath, 'Adapter', CZarrAdapter);
    if any(bim.BlockSize ~= blockSize) || any(bim.Size ~= rs_size)
        rmdir(rsTmpPath, 's');
        rmdir(zarrFlagPath, 's');
        mkdir(zarrFlagPath);
    end
end
if ~exist(rsTmpPath, 'dir')
    dimSeparator = '.';
    if prod(ceil(rs_size ./ blockSize)) > 10000
        dimSeparator = '/';
    end
    createzarr(rsTmpPath, dataSize=rs_size, blockSize=blockSize, dtype=dtype, ...
        dimSeparator=dimSeparator);
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
        '''interpMethod'',''%s'',''inputBbox'',[%s],''borderSize'',[%s],''batchSize'',[%s],', ...
        '''blockSize'',%s)'], strrep(num2str(batchInds, '%d,'), ' ', ''), zarrFullpath, ...
        rsTmpPath, zarrFlagFullpath, strrep(num2str(resampleFactor(:)', '%d,'), ' ', ''), interpMethod, ...
        strrep(num2str(inputBbox(:)', '%d,'), ' ', ''), strrep(num2str(borderSize(:)', '%d,'), ' ', ''), ...
        strrep(num2str(batchSize(:)', '%d,'), ' ', ''), strrep(mat2str(blockSize(:)'), ' ', ','));
end

inputFullpaths = repmat({zarrFullpath}, numTasks, 1);

dtype = getImageDataType(zarrFullpath);
byteNum = dataTypeToByteNumber(dtype);
memAllocate = prod(batchSize) * byteNum / 2^30 * (prod(resampleFactor) + 1) * 3;
is_done_flag= generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
    funcStrs, finalOutFullpath=rsFullpath, cpusPerTask=cpusPerTask, memAllocate=memAllocate, ...
    parseCluster=parseCluster, mccMode=mccMode, configFile=configFile);

if ~all(is_done_flag)
    is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, finalOutFullpath=rsFullpath, cpusPerTask=cpusPerTask, memAllocate=memAllocate, ...
        parseCluster=parseCluster, mccMode=mccMode, configFile=configFile);
end

if ~all(is_done_flag)
    error('%d / %d tasks failed!', sum(~is_done_flag), numel(is_done_flag));
end

if exist(rsFullpath, 'dir')
    rmdir(rsFullpath, 's');
end
movefile(rsTmpPath, rsFullpath);
rmdir(zarrFlagPath, 's');
% fprintf('Done!\n')

end

