function [] = XR_MIP_wrapper(dataPaths, varargin)   
% generate MIP for dataset

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @(x) ischar(x) || iscell(x));
ip.addParameter('axis', [0, 0, 1], @isnumeric); % y, x, z
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('largeZarr', false, @islogical); % use zarr file as input
ip.addParameter('Save16bit', true, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpuOnlyNodes', ~true, @islogical);
ip.addParameter('cpusPerTask', 3, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
axis = pr.axis;
ChannelPatterns =  pr.ChannelPatterns;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
Save16bit = pr.Save16bit;
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

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

if Save16bit
    dtype = 'uint16';
else
    dtype = 'single';
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    resultPath = [dataPath, '/MIPs/'];
    mkdir(resultPath);
    
    dataPaths{d} = simplifyPath(dataPath);
    resultPaths{d} = simplifyPath(resultPath);    
    fileattrib(resultPath, '+w', 'g');
    save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str, jobLogDir] = checkSlurmCluster(dataPaths{1}, jobLogDir, cpuOnlyNodes);
end

fnames = cell(nd, 1);

for d = 1 : nd
    dataPath = dataPaths{d};
    if zarrFile
        dir_info = dir([dataPath, '/', '*.zarr']);        
    else
        dir_info = dir([dataPath, '/', '*.tif']);
    end
    fnames_d = {dir_info.name}';
    nF = numel(fnames_d);
    sz = getImageSize([dataPath, '/', fnames_d{1}]);    
    nc = numel(ChannelPatterns);
    ch_inds = false(nF, nc);
    for c = 1 : nc
        ch_inds(:, c) = contains(fnames_d, ChannelPatterns{c}, 'IgnoreCase', true) | contains(fnames_d, regexpPattern(ChannelPatterns{c}), 'IgnoreCase', true);
    end
    fnames_d = fnames_d(any(ch_inds, 2));
    fnames{d} = fnames_d;
end

fd_inds = arrayfun(@(x) ones(numel(fnames{x}), 1) * x, 1 : nd, 'unif', 0);
fnames = cat(1, fnames{:});
[~, fsns] = fileparts(fnames);
fd_inds = cat(1, fd_inds{:});
nF = numel(fnames);
if nF == 1
    fsns = {fsns};
end

%% use generic framework for the MIP computing

frameFullpaths = cell(nF, 1);
MIPFullpaths = cell(nF, 1);
func_strs = cell(nF, 1);

% use the last select axis as flag if more than one axis for MIP
axis_strs = {'y', 'x', 'z'};
fidx = find(axis, 1, 'last');
axis_str = axis_strs{fidx};

for f = 1 : nF
    fname = fnames{f};
    d = fd_inds(f);
    dataPath = dataPaths{d};
    resultPath = resultPaths{d};
    
    frameFullpath = [dataPath, '/', fname];
    frameFullpaths{f} = frameFullpath;
    MIPFullpath = sprintf('%s/%s_MIP_%s.tif', resultPath, fsns{f}, axis_str);
    MIPFullpaths{f} = MIPFullpath;
    
    if zarrFile
        if largeZarr || any(axis(1 : 2))
            func_strs{f} = sprintf(['XR_MIP_zarr(''%s'',''axis'',%s,''parseCluster'',%s,', ...
                '''parseParfor'',%s)'], frameFullpath, strrep(mat2str(axis), ' ', ','), ...
                string(parseCluster), string(parseParfor));       
        else
            func_strs{f} = sprintf(['saveMIP_zarr(''%s'',''%s'',''%s'')'], frameFullpath, ...
                MIPFullpath, dtype);
        end
    else
        func_strs{f} = sprintf(['saveMIP_tiff(''%s'',''%s'',''dtype'',''%s'',''axis'',%s)'], ...
            frameFullpath, MIPFullpath, dtype, strrep(mat2str(axis), ' ', ','));
    end
end

if parseParfor && ~largeZarr
    matlab_parfor_generic_computing_wrapper(frameFullpaths, MIPFullpaths, func_strs, 'GPUJob', false, 'nworker', 12)
end
taskBatchNum = max(1, round(nF /1000));
cpusPerTask = max(min(ceil(prod(sz) * 4 / 1024^3 * 4 / 20), 24), cpusPerTask);
slurm_cluster_generic_computing_wrapper(frameFullpaths, MIPFullpaths, func_strs, ...
    'parseCluster', parseCluster, 'masterCompute', masterCompute, 'taskBatchNum', taskBatchNum, ...
    'cpusPerTask', cpusPerTask);

end

