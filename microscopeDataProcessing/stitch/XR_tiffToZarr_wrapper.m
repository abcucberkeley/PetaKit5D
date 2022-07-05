function [] = XR_tiffToZarr_wrapper(tiffFullpaths, varargin)
% The wrapper for convert a list of tiff files (a data folder) to zarr
%
% Author: Xiongtao Ruan (10/02/2020)
% 
% xruan (10/11/2020): add function handle for processing before saving to zarr
% xruan (07/26/2021): add support for flipped tile
% xruan (07/27/2021): add support for resampling; simplify cluster job submission.
% xruan (08/25/2021): add support for channel-specific user functions
% xruan (09/23/2021): add support for including partial files 
% xruan (10/13/2021): add support for cropping data
% xruan (02/16/2022): accelerate the code by first get filenames for every data folder.
% xruan (07/05/2022): add support for single tiff file (char) conversion

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tiffFullpaths', @(x) iscell(x) || ischar(x));
ip.addParameter('zarrPathstr', 'zarr', @ischar);
ip.addParameter('blockSize', [500, 500, 250], @isnumeric);
ip.addParameter('flippedTile', [], @(x) isempty(x) || islogical(x));
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('partialFile', false, @islogical);
ip.addParameter('ChannelPatterns', {'tif'}, @iscell);
ip.addParameter('InputBbox', [], @isnumeric); % crop input tile before processing
ip.addParameter('CropToSize', [], @isnumeric); % size cropped to 
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x) || iscell(x));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);


ip.parse(tiffFullpaths, varargin{:});

pr = ip.Results;
 % Resolution = pr.Resolution;
zarrPathstr = pr.zarrPathstr;
blockSize = pr.blockSize;
flippedTile = pr.flippedTile;
resample = pr.resample;
partialFile = pr.partialFile;
ChannelPatterns = pr.ChannelPatterns;
InputBbox = pr.InputBbox;
CropToSize = pr.CropToSize;
usrFcn = pr.usrFcn;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes = pr.cpuOnlyNodes;


if ischar(tiffFullpaths)
    if exist(tiffFullpaths, 'dir')
        dir_info = dir([tiffFullpaths, filesep, '*.tif']);
        fnames = {dir_info.name}';
        tiffFullpaths = cellfun(@(x) sprintf('%s/%s', tiffFullpaths, x), fnames, 'unif', 0);
    else
        tiffFullpaths = {tiffFullpaths};
    end
end

[dataPath, ~] = fileparts(tiffFullpaths{1});
    
% check if a slurm-based computing cluster exist
if parseCluster 
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
end

nC = numel(ChannelPatterns);
usrFcn_strs = repmat({''}, nC, 1);
if ~isempty(usrFcn)
    if isa(usrFcn,'function_handle')
        usrFcn_strs = repmat({func2str(usrFcn)}, nC, 1);
    elseif ischar(usrFcn)
        usrFcn_strs = repmat({usrFcn}, nC, 1);
    elseif iscell(usrFcn)
        if isa(usrFcn{1},'function_handle')
            usrFcn_strs = cellfun(@(x) func2str(x), usrFcn, 'unif', 0);
        elseif ischar(usrFcn{1})
            usrFcn_strs = usrFcn;
        end
    end
end

nF = numel(tiffFullpaths);
if partialFile
    dataPaths = fileparts(tiffFullpaths);
    udataPaths = unique(dataPaths);
    for d = 1 : numel(udataPaths)
        dir_info = dir([udataPaths{d}, '/*tif']);
        fsnames = {dir_info.name}';
        fnames = cellfun(@(x) [udataPaths{d}, filesep, x], fsnames, 'unif', 0);
    end
end

zarrFullpaths = cell(nF, 1);
func_strs = cell(nF, 1);
for i = 1 : nF
    tiffFullpath_i = tiffFullpaths{i};
    [dataPath, fsname_i] = fileparts(tiffFullpath_i);
    
    if partialFile
        tiffFullpath_group_i = fnames(contains(fnames, dataPath) & contains(fnames, fsname_i));
    else
        tiffFullpath_group_i = {tiffFullpath_i};        
    end

    zarrPath = [dataPath, filesep, zarrPathstr];
    if ~exist(zarrPath, 'dir')
        mkdir(zarrPath);
    end

    zarrFullpaths{i} = sprintf('%s/%s/%s.zarr/', dataPath, zarrPathstr, fsname_i);
    if ~isempty(flippedTile)
        flipZstack = flippedTile(i);
    else
        flipZstack = false;
    end
    
    cind = cellfun(@(x) contains(tiffFullpath_i, x), ChannelPatterns);
    if ~any(cind)
        error('The file %s does not match any channel patterns %s', tiffFullpath_i, string(ChannelPatterns));
    end
    
    func_strs{i} = sprintf(['tiffToZarr(%s,''%s'',[],''BlockSize'',%s,''flipZstack'',%s,', ...
        '''resample'',%s,''InputBbox'',%s,''CropToSize'',%s,''usrFcn'',''%s'')'], ...
        sprintf('{''%s''}', strjoin(tiffFullpath_group_i, ''',''')), zarrFullpaths{i}, ...
        strrep(mat2str(blockSize), ' ', ','), string(flipZstack), strrep(mat2str(resample), ' ', ','), ...
        strrep(mat2str(InputBbox), ' ', ','), strrep(mat2str(CropToSize), ' ', ','), usrFcn_strs{cind});
end

[estMem, estGPUMem, rawImageSize] = XR_estimateComputingMemory(tiffFullpaths{1}, {'deconvolution'}, 'cudaDecon', false);
if cpusPerTask * 20 < rawImageSize * 4 * numel(tiffFullpath_group_i)
    cpusPerTask = min(24, ceil(rawImageSize * 4 * numel(tiffFullpath_group_i) / 20));
end
maxTrialNum = 2;

MatlabLaunchStr = 'module load matlab/r2021a; matlab -nodisplay -nosplash -nodesktop -r';
is_done_flag = slurm_cluster_generic_computing_wrapper(tiffFullpaths, zarrFullpaths, ...
    func_strs, 'parseCluster', parseCluster, 'masterCompute', true, 'maxTrialNum', maxTrialNum, ...
    'MatlabLaunchStr', MatlabLaunchStr, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes);
if ~all(is_done_flag)
    slurm_cluster_generic_computing_wrapper(tiffFullpaths, zarrFullpaths, ...
        func_strs, 'parseCluster', parseCluster, 'masterCompute', true, ...
        'maxTrialNum', maxTrialNum, 'MatlabLaunchStr', MatlabLaunchStr, ...
        'cpusPerTask', min(24, cpusPerTask * 2), 'cpuOnlyNodes', cpuOnlyNodes);
end


end


