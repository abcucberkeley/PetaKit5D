function [] = XR_resampleSingleZarr(zarrFullpath, dsFullpath, dsFactor, varargin)
% 
% Author: Xiongtao Ruan (12/14/2020)
% xruan (11/09/2021): enable arbitray blockSize


if nargin < 1
   zarrFullpath = '/Users/xruan/Images/20201204_p465_p5_LLCPK_488nmSec61_560nm_H2B_RG0p2/tilingTest2_z0p4/matlab_stitch_xcorr_feather_zarr/Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0032542619msecAbs_crop.zarr';
   dsFullpath = [zarrFullpath(1 : end - 5), '_ds.zarr'];
   dsFactor = [2, 2, 2];
end

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath'); 
ip.addRequired('dsFullpath'); 
ip.addRequired('dsFactors'); 
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @isnumeric); % size to process in one batch 
ip.addParameter('BorderSize', [5, 5, 5], @isnumeric); % padded boarder for each batch
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear', 'nearest'})));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('cpuOnlyNodes', ~true, @islogical);
ip.addParameter('cpusPerTask', 1, @islogical);
ip.addParameter('uuid', '', @ischar);

ip.parse(zarrFullpath, dsFullpath, dsFactor, varargin{:});

pr = ip.Results;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
BorderSize = pr.BorderSize;
Interp = pr.Interp;
parseCluster = pr.parseCluster;
cpuOnlyNodes = pr.cpuOnlyNodes;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;

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

try
    bim = blockedImage(zarrFullpath, "Adapter", CZarrAdapter);
catch ME
    disp(ME);
    bim = blockedImage(zarrFullpath, "Adapter", ZarrAdapter);
end    
dtype = bim.ClassUnderlying;
sz = bim.Size;
init_val = zeros(1, dtype);

ds_size = round(sz ./ [dsFactor(1), dsFactor(2), dsFactor(3)]);
try
    ds_bim = blockedImage(dsTmppath, ds_size, blockSize, init_val, "Adapter", CZarrAdapter, 'Mode', 'w');
catch ME
    disp(ME);
    ds_bim = blockedImage(dsTmppath, ds_size, blockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');
end
ds_bim.Adapter.close();

% framework
blockSize = min(ds_size, blockSize);
batchSize = min(ds_size, batchSize);
batchSize = ceil(batchSize ./ blockSize) .* blockSize;
bSubSz = ceil(ds_size ./ batchSize);
numBatch = prod(bSubSz);

% process for each block based on all BlockInfo use distributed computing
fprintf('Process blocks for downsampled data...\n')

zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', dsPath, dsfsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
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
    
    funcStrs{i} = sprintf(['resampleZarrBlock([%s],''%s'',''%s'',''%s'',[%s],''Interp'',''%s'',', ...
        '''BorderSize'',[%s],''batchSize'',[%s],''blockSize'',%s)'], strrep(num2str(batchInds, '%d,'), ' ', ''), ...
        zarrFullpath, dsTmppath, zarrFlagFullpath, strrep(num2str(dsFactor(:)', '%d,'), ' ', ''), Interp, ...
        strrep(num2str(BorderSize(:)', '%d,'), ' ', ''), strrep(num2str(batchSize(:)', '%d,'), ' ', ''), ...
        strrep(mat2str(blockSize(:)'), ' ', ','));
end

inputFullpaths = repmat({zarrFullpath}, numTasks, 1);

is_done_flag= slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
    funcStrs, 'cpusPerTask', cpusPerTask, 'parseCluster', parseCluster, 'cpuOnlyNodes', cpuOnlyNodes);

if ~all(is_done_flag)
    slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', cpusPerTask * 2, 'parseCluster', parseCluster, ...
        'cpuOnlyNodes', cpuOnlyNodes);
end    

if exist(dsFullpath, 'dir')
    rmdir(dsFullpath, 's');
end
movefile(dsTmppath, dsFullpath);
rmdir(zarrFlagPath, 's');

end

