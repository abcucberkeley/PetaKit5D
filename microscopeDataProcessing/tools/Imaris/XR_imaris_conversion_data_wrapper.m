function XR_imaris_conversion_data_wrapper(dataPaths, varargin)
% wrapper to convert dataset to ims file
%
% Author: Xiongtao Ruan (06/23/2022)
% 
% xruan (07/28/2022): add support for user defined chunk size



ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'imaris',  @(x) ischar(x));
ip.addParameter('overwrite', false,  @(x) islogical(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('pixelSizes', [0.108, 0.108, 0.108], @isnumeric); % y, x, z
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('blockSize', [64, 64, 64], @isnumeric); % y, x, z
ip.addParameter('inputBbox', [], @isnumeric); % ymin, xmin, zmin, ymax, xmax, zmax
ip.addParameter('timepoints', [], @isnumeric); % number of time points included
ip.addParameter('converterPath', '', @ischar);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 24, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

% parameters
pr = ip.Results;
resultDirName = pr.resultDirName;
overwrite = pr.overwrite;
channelPatterns = pr.channelPatterns;
pixelSizes = pr.pixelSizes;
zarrFile = pr.zarrFile;
blockSize = pr.blockSize;
inputBbox = pr.inputBbox;
timepoints = pr.timepoints;
converterPath = pr.converterPath;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

% if converter path is not provided, use the default one in the repo
if isempty(converterPath)
    mfilePath = fileparts(which(mfilename));
    if ispc
        mfilePath = strrep(mfilePath, '\', '/');
        converterPath = sprintf('%s/Parallel_Imaris_Writer/windows/parallelimariswriter', mfilePath);        
    elseif ismac
        converterPath = sprintf('%s/Parallel_Imaris_Writer/mac/parallelimariswriter', mfilePath);
    else
        converterPath = sprintf('%s/Parallel_Imaris_Writer/linux/parallelimariswriter', mfilePath);
    end
end

% suppress directory exists warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ispc
        dataPath = strrep(dataPath, '\', '/');
    end
    if ~strcmp(dataPath(end), '/')
        dataPaths{d} = [dataPath, '/'];
    end
end

if numel(overwrite) == 1
    overwrite = repmat(overwrite, 1, 2);
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str] = checkSlurmCluster(dataPath, jobLogDir);
end

% handle output directories.
imsName = 'imaris';
if ~isempty(resultDirName)
    imsName = resultDirName;
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

    imsPath = [dataPath_fst, imsName, '/'];
    if overwrite(1) && exist(imsPath, 'dir')
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

channelPatterns_str = strjoin(channelPatterns, ',');
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

bbox_str = '';
if ~isempty(inputBbox)
    bbox_str = sprintf('-B %s', strrep(num2str(inputBbox, '%d,'), ' ', ''));
end

func_strs = cell(nd, 1);
for d = 1 : nd
    dataPath_str = dataPaths{d};
    if iscell(dataPaths{d})
        dataPath_str = strjoin(dataPaths{d}, ',');
    end 
    
    func_strs{d} = sprintf('"%s" -P "%s" -F "%s" -r %s -v %s -o "%s" %s -b %s %s', converterPath, ...
        channelPatterns_str, dataPath_str, type_str, voxelsize_str, imsPaths{d}, time_str, ...
        blockSize_str, bbox_str);
end

% use slurm job wrapper for computing
generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, func_strs, ...
    'cpusPerTask', cpusPerTask, 'parseCluster', parseCluster, 'language', 'bash', ...
    mccMode=mccMode, configFile=configFile);

end

