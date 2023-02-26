function [] = XR_parallel_rsync_wrapper(source, dest, varargin)
% rsync with multiple workers (12/21/2020)
%
% Author: Xiongtao Ruan (12/21/2020)
% xruan (11/04/2021): add support for file shuffling
% xruan (09/08/2022): add support for user defined nice factor for priority

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('source', @(x) ischar(x));
ip.addRequired('dest', @(x) ischar(x));
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('numStream', inf, @isnumeric); % number of streams for the transfer
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar); % slurm parameters 
ip.addParameter('niceFactor', 0, @isnumeric); % nice factor for priority 
ip.addParameter('includeSubdir', true, @islogical); % also transfer subfolders

ip.parse(source, dest, varargin{:});

pr = ip.Results;
cpusPerTask = pr.cpusPerTask; 
numStream = pr.numStream;
SlurmParam = pr.SlurmParam;
niceFactor = pr.niceFactor;
includeSubdir = pr.includeSubdir;
SlurmParam = sprintf('%s --nice=%d', SlurmParam, niceFactor);

tic
% go to the repo root folder
setupFn = which('setup.m');
[codePath, ~, ~] = fileparts(setupFn);
cd(codePath);
if ~contains(codePath, 'LLSM5DTools')
    cd('LLSM5DTools');
end
codePath = pwd;

uuid = get_uuid();

% first make the corresponding directory structures in the destination
% cmd = sprintf("rsync -av --include='*/' --exclude='*' %s %s", source, dest);
% system(cmd, '-echo');
% parfeval(@system, 2, cmd, '-echo');

source = strrep(source, '//', '/');
if includeSubdir
    inputFullpaths = dir_recursive(source, 'file', true);
else
    dir_info = dir([source, '/*']);
    fsns = {dir_info.name}';
    fsns = fsns(~[dir_info.isdir]);
    inputFullpaths = cellfun(@(x) [source, '/', x], fsns, 'UniformOutput', false);
end
nF = numel(inputFullpaths);
% obtain the proposed output paths
dest = strrep(dest, '//', '/');
if ~strcmp(dest(end), '/')
    dest = [dest, '/'];
end

if ~strcmp(source(end), '/')
    [p, d, ext] = fileparts(source);
    source = [source, '/'];
    dest = [dest, d, ext filesep];
end

outputFullpaths = strrep(inputFullpaths, source, dest);
outputPaths = fileparts(outputFullpaths);

% handle file/folder with space
inputFullpaths = strrep(inputFullpaths, ' ', '\ ');
outputFullpaths = strrep(outputFullpaths, ' ', '\ ');

% split into batches (argument list cannot be too long for a single batch)
batchSize = max(1, round(min(nF / 500, min(max(5, nF / 10000), 200))));
% batchSize = 1;
numBatch = ceil(nF / batchSize);

% set flag files
tmpDir = sprintf('%s/tmp/%s/', dest, uuid);
mkdir(tmpDir);

% random shuffling of the files so that large files are distributed to
% different workers evenly.
rng(1);
rand_inds = randperm(nF);
inputFullpaths = inputFullpaths(rand_inds);
outputFullpaths = outputFullpaths(rand_inds);
outputPaths = outputPaths(rand_inds);

% check if all files are transferred, if the number of files is 100k or
% fewer, direct check if they are exist in destination, otherwise, use
% cluster based check
for n = 1 : 5
    if numel(outputPaths) < 100000
        out_exist_mat = batch_file_exist(outputFullpaths);
    else
        t = datetime('now', 'format', 'yMMddHHmmss');
        exist_info_dir = sprintf('%s/%s/', tmpDir, t);        
        out_exist_mat = batch_file_exist_cluster(outputFullpaths, exist_info_dir);
    end

    if all(out_exist_mat)
        break;
    end
    
    inputFullpaths = inputFullpaths(~out_exist_mat);
    outputFullpaths = outputFullpaths(~out_exist_mat);
    outputPaths = outputPaths(~out_exist_mat);
    flagFullpaths = outputFullpaths;
    
    nF = numel(inputFullpaths);

    % split into batches
    batchSize = max(1, round(min(nF / 500, min(max(5, nF / 10000), 200))));
    % batchSize = 1;
    numBatch = ceil(nF / batchSize);

    func_strs = cell(numBatch, 1);
    for i = 1 : numBatch
        s = (i - 1) * batchSize + 1;
        t = min(i * batchSize, nF);

        inputFullpath_str = sprintf('%s', strjoin(inputFullpaths(s : t), ','));
        outPath_str = sprintf('%s', strjoin(outputPaths(s : t), ','));

        func_strs{i} = sprintf(['bash %s/microscopeDataProcessing/io/rsync_wrapper/rsync_batch_files.sh ', ...
            '''''%s'''' ''''%s'''' '], codePath, inputFullpath_str, outPath_str);
    end

    inputFlagpaths = inputFullpaths(max(1, batchSize - 1) : batchSize : end);
    outputFlagpaths = flagFullpaths(max(1, batchSize - 1) : batchSize : end);
    maxJobNum = numStream;

    is_done_flag = slurm_cluster_generic_computing_wrapper(inputFlagpaths, ...
        outputFlagpaths, func_strs, 'tmpDir', tmpDir, 'language', 'bash', ...
        'SlurmParam', SlurmParam, 'cpusPerTask', cpusPerTask, 'maxJobNum', maxJobNum, ...
        'masterCompute', ~false);
end

% final check with all file rsync 
cmd = sprintf("rsync -avP %s %s", source, dest);
system(cmd, '-echo');

rmdir(tmpDir, 's');
toc


end

