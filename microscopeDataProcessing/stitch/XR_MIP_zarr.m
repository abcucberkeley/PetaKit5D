function XR_MIP_zarr(zarrFullpath, varargin)
% save MIP for large scale zarr file. The idea is to first generate MIPs
% for each batch for all three axis, and then generate final MIPs. 
% 
% 
% Author: Xiongtao Ruan (02/17/2022)



ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath', @(x) ischar(x) || iscell(x));
ip.addParameter('axis', [0, 0, 1], @isnumeric); % y, x, z
ip.addParameter('BatchSize', [2048, 2048, 2048] , @isvector); % in y, x, z
ip.addParameter('BlockSize', [2048, 2048, 2048] , @isvector); % in y, x, z
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpuOnlyNodes', ~true, @islogical);
ip.addParameter('cpusPerTask', 3, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(zarrFullpath, varargin{:});

pr = ip.Results;
axis = pr.axis;
BatchSize = pr.BatchSize;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpuOnlyNodes = pr.cpuOnlyNodes;
cpusPerTask = pr.cpusPerTask;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;

zarrFullpath = strip(zarrFullpath, 'right', filesep);
[dataPath, fsname, ext] = fileparts(zarrFullpath);

MIPPath = [dataPath, '/MIPs/'];
if ~exist(MIPPath, 'dir')
    mkdir(MIPPath);
end

axis_strs = {'y', 'x', 'z'};
MIPFullpaths = cellfun(@(x) sprintf('%s/%s_MIP_%s.tif', MIPPath, fsname, x), axis_strs, 'unif', 0);

done_flag = false(3, 1);
for i = 1 : 3
    done_flag(i) = (axis(i) == 0) | exist(MIPFullpaths{i}, 'file');
end
if all(done_flag)
    disp('The output results exist, skip it!');
    return;
end

fprintf('Start Large-file MIP for %s...\n', fsname);

if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
end

tic
zarrFlagPath = sprintf('%s/zarr_flag/%s_%s/', MIPPath, fsname, uuid);
if ~exist(zarrFlagPath, 'dir')
    mkdir_recursive(zarrFlagPath);
end 

try
    bim = blockedImage(zarrFullpath, 'Adapter', CZarrAdapter);
catch ME
    disp(ME);
    bim = blockedImage(zarrFullpath, 'Adapter', ZarrAdapter);    
end
imSize = bim.Size;
dtype = bim.ClassUnderlying;

% MIPs for each block
BatchSize = min(imSize, BatchSize);
bSubSz = ceil(imSize ./ BatchSize);
numBatch = prod(bSubSz);

[Y, X, Z] = ndgrid(1 : bSubSz(1), 1 : bSubSz(2), 1 : bSubSz(3));
bSubs = [Y(:), X(:), Z(:)];
clear Y X Z

batchBBoxes = zeros(numBatch, 6);
batchBBoxes(:, 1 : 3) = (bSubs - 1) .* BatchSize + 1; 
batchBBoxes(:, 4 : 6) = min(batchBBoxes(:, 1 : 3) + BatchSize - 1, imSize);

MIPZarrpaths = cellfun(@(x) sprintf('%s/%s_%s.zarr', MIPPath, fsname, x), axis_strs, 'unif', 0);
MIPZarrTmppaths = cellfun(@(x) sprintf('%s/%s_%s_%s.zarr', MIPPath, fsname, x, uuid), axis_strs, 'unif', 0);

zarr_done_flag = false(3, 1);
for i = 1 : 3
    zarr_done_flag(i) = (axis(i) == 0) | exist(MIPZarrpaths{i}, 'dir');
end

% if all mip zarr files exist, directly generate final MIPs
if all(zarr_done_flag)
    % collect results and generate MIPs 
    for i = 1 : 3
        if axis(i) == 0
            continue;
        end

        saveMIP_zarr(MIPZarrpaths{i}, MIPFullpaths{i}, dtype, (1 : 3) == i);
    end
    return;
end

% initialize zarr files
init_val = zeros(1, dtype);
for i = 1 : 3
    if exist(MIPZarrTmppaths{i}, 'dir')
        continue;
    end
    axis_flag = false(3, 1);
    axis_flag(i) = true;
    
    outSize = imSize;
    outSize(axis_flag) = bSubSz(axis_flag);
    BlockSize = BatchSize;
    BlockSize(axis_flag) = 1;
    
    try
        mip_bim = blockedImage(MIPZarrTmppaths{i}, outSize, BlockSize, init_val, "Adapter", CZarrAdapter, 'Mode', 'w');
    catch ME
        disp(ME);
        mip_bim = blockedImage(MIPZarrTmppaths{i}, outSize, BlockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');
    end
    mip_bim.Adapter.close();
end

% set up parallel computing 
numBatch = size(batchBBoxes, 1);
taskSize = max(2, min(10, round(numBatch / 5000))); % the number of batches a job should process
numTasks = ceil(numBatch / taskSize);

maxJobNum = inf;
taskBatchNum = 1;
SlurmParam = '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M';

% get the function string for each batch
funcStrs = cell(numTasks, 1);
outputFullpaths = cell(numTasks, 1);
for i = 1 : numTasks
    batchInds = (i - 1) * taskSize + 1 : min(i * taskSize, numBatch);
    batchBBoxes_i = batchBBoxes(batchInds, :);
    bSubs_i = bSubs(batchInds, :);
    
    zarrFlagFullpath = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    outputFullpaths{i} = zarrFlagFullpath;
    MIPZarrTmppaths_str = sprintf('{''%s''}', strjoin(MIPZarrTmppaths, ''','''));
    
    funcStrs{i} = sprintf(['MIP_block([%s],''%s'',%s,''%s'',%s,%s,''uuid'',''%s'',', ...
        '''debug'',%s)'], strrep(num2str(batchInds, '%d,'), ' ', ''), zarrFullpath, ...
        MIPZarrTmppaths_str, zarrFlagFullpath, strrep(mat2str(batchBBoxes_i), ' ', ','), ...
        strrep(mat2str(bSubs_i), ' ', ','), uuid, string(debug));
end

% submit jobs 
inputFullpaths = repmat({zarrFullpath}, numTasks, 1);
if parseCluster || ~parseParfor 
    is_done_flag= slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
        'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);

    if ~all(is_done_flag)
        slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
            'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);
    end
elseif parseParfor
    GPUJob = false;
    nworker = 12;
    is_done_flag= matlab_parfor_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'taskBatchNum', taskBatchNum, 'GPUJob', GPUJob, 'nworker', nworker, 'uuid', uuid);
end

if ~all(is_done_flag)
    error('Some blocks are not finished!')
end

for i = 1 : 3
    if exist(MIPZarrpaths{i}, 'dir') && exist(MIPZarrTmppaths{i}, 'dir')
        rmdir(MIPZarrpaths{i}, 's');
    end
    if exist(MIPZarrTmppaths{i}, 'dir')
        movefile(MIPZarrTmppaths{i}, MIPZarrpaths{i});
    end
end

% collect results and generate MIPs 
for i = 1 : 3
    if axis(i) == 0
        continue;
    end
    
    saveMIP_zarr(MIPZarrpaths{i}, MIPFullpaths{i}, dtype, (1 : 3) == i);
end

end

