function XR_imaris_conversion_data_wrapper(dataPaths, varargin)
% wrapper to convert dataset to ims file
%
% Author: Xiongtao Ruan (06/23/2022)
% 
% xruan (07/28/2022): add support for user defined chunk size



ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths');
ip.addParameter('imsPathstr', 'imaris',  @(x) ischar(x));
ip.addParameter('Overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 2) && islogical(x));
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('pixelSizes', [0.108, 0.108, 0.108], @isnumeric); % y, x, z
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('blockSize', [64, 64, 64], @isnumeric); % y, x, z
ip.addParameter('timepoints', [], @isnumeric); % number of time points included
ip.addParameter('ImsConverter', '/clusterfs/fiona/matthewmueller/imarisWriter/writeImarisParallel', @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 24, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);

ip.parse(dataPaths, varargin{:});

% parameters
pr = ip.Results;
imsPathstr = pr.imsPathstr;
Overwrite = pr.Overwrite;
ChannelPatterns = pr.ChannelPatterns;
pixelSizes = pr.pixelSizes;
zarrFile = pr.zarrFile;
blockSize = pr.blockSize;
timepoints = pr.timepoints;
ImsConverter = pr.ImsConverter;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
SlurmParam = pr.SlurmParam;

% suppress directory exists warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), filesep)
        dataPaths{d} = [dataPath, filesep];
    end
end

if numel(Overwrite) == 1
    Overwrite = repmat(Overwrite, 1, 2);
end

% check if a slurm-based computing cluster exists
if parseCluster
    cpuOnlyNodes = false;
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str, jobLogDir] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
end

% handle output directories.
imsName = 'imaris';
if ~isempty(imsPathstr)
    imsName = imsPathstr;
end
    
imsPaths = cell(nd, 1);
inputFullpaths = cell(nd, 1);
outputFullpaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    dataPath_fst = dataPath;        
    if iscell(dataPath)
        dataPath_fst = dataPath{1};
    end
    
    dir_info = dir([dataPath_fst, '*tif']);
    if zarrFile
        dir_info = dir([dataPath_fst, '*zarr']);
    end
    fsn = {dir_info.name}';
    if ~isempty(fsn)
        inputFullpaths{d} = sprintf('%s/%s', dataPath_fst, fsn{1});
    end

    imsPath = [dataPath_fst, imsName, filesep];
    if Overwrite(1) && exist(imsPath, 'dir')
        rmdir(imsPath, 's');
    end
    if ~exist(imsPath, 'dir')
        mkdir(imsPath);
        fileattrib(imsPath, '+w', 'g');
    end
    imsPaths{d} = imsPath;
    outputFullpaths{d} = [imsPath, 'output.ims'];
    
    % save decon parameters
    save('-v7.3', [imsPath, '/parameters.mat'], 'pr');
    writetable(struct2table(pr, 'AsArray', true), [imsPath, '/parameters.txt'])
end

%% setup jobs and computing

ChannelPatterns_str = strjoin(ChannelPatterns, ',');
% voxelsize_str = strrep(mat2str(pixelSizes), ' ', ',');
voxelsize_str = num2str(pixelSizes, '%.20f,');
type_str = 'tiff';
if zarrFile
    type_str = 'zarr';
end
% optional flags, also include flag str
time_str = '';
if ~isempty(timepoints)
    time_str = sprintf('-t %d', timepoints);
end
blockSize_str = num2str(blockSize, '%d,'); 

func_strs = cell(nd, 1);
for d = 1 : nd
    dataPath_str = dataPaths{d};
    if iscell(dataPaths{d})
        dataPath_str = strjoin(dataPaths{d}, ',');
    end 

    func_strs{d} = sprintf('%s -P %s -F %s -r %s -v %s -o %s %s -b %s', ImsConverter, ...
        ChannelPatterns_str, dataPath_str, type_str, voxelsize_str, imsPaths{d}, time_str, ...
        blockSize_str);
end

% use slurm job wrapper for computing
slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, func_strs, ...
    'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', false, 'parseCluster', parseCluster, ...
    'SlurmParam', SlurmParam, 'language', 'bash');

end


