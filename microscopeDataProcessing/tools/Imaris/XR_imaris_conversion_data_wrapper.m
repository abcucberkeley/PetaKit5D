function XR_imaris_conversion_data_wrapper(dataPaths, varargin)
% Dataset level wrapper to convert dataset to Imaris file.
%
%
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%
% Parameters (as 'specifier'-value pairs):
%       resultDirName : char (default: 'imaris'). Result directory under data paths.
%           overwrite : true|false (default: false). Overwrite existing results
%     channelPatterns : a cell array (default: {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}).  Channel identifiers for included channels. 
%          pixelSizes : 1x3 vector (default: [0.108, 0.108, 0.108]). Pixel sizes in um, in yxz order.
%            zarrFile : true|false (default: false). Use Zarr file as input.
%           blockSize : 1x3 vector (default: [64, 64, 64]). Block/chunk size for Imaris output.
%           inputBbox : empty or 1x6 vector (default: []). Input bounding box for crop. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%       converterPath : empty or char (default: ''). Converter binary path. If empty, use the default one within the software.
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%           jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%         cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%
%
% Author: Xiongtao Ruan (06/23/2022)


ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'imaris',  @(x) ischar(x));
ip.addParameter('overwrite', false,  @(x) islogical(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('pixelSizes', [0.108, 0.108, 0.108], @isnumeric);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('blockSize', [64, 64, 64], @isnumeric);
ip.addParameter('inputBbox', [], @isnumeric);
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
converterPath = pr.converterPath;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if isempty(uuid)
    if ispc
        uuid = num2str(randi(1024));
    else
        uuid = get_uuid();
    end
end

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
    else
        dataPaths{d} = dataPath;
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
tmpFsns = cell(nd, 1);
tmpFullpaths = cell(nd, 1);
outputFullpaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    dataPath_fst = dataPath;        
    if iscell(dataPath)
        dataPath_fst = dataPath{1};
    end
    
    if zarrFile
        dir_info = dir([dataPath_fst, '*zarr']);
    else
        dir_info = dir([dataPath_fst, '*tif']);        
    end
    fsn = {dir_info.name}';
    if ~isempty(fsn)
        inputFullpaths{d} = sprintf('%s/%s', dataPath_fst, fsn{1});
    end
    [~, ffsn] = fileparts(fsn{1});

    imsPath = [dataPath_fst, imsName, '/'];
    if overwrite(1) && exist(imsPath, 'dir')
        rmdir(imsPath, 's');
    end
    if ~exist(imsPath, 'dir')
        mkdir(imsPath);
        fileattrib(imsPath, '+w', 'g');
    end
    imsPaths{d} = imsPath;
    if ~ispc
        if strcmp(imsPath(1 : 2), '~/')
            homedir = getenv('HOME');
            imsPath = sprintf('%s%s', homedir, imsPath(2 : end));
        end
    end
    tmpFsns{d} = sprintf('%s_%s', ffsn, uuid);
    tmpFullpaths{d} = sprintf('%s%s_%s.ims', imsPath, ffsn, uuid);
    outputFullpaths{d} = sprintf('%s%s.ims', imsPath, ffsn);
    
    % save decon parameters
    save('-v7.3', [imsPath, 'parameters.mat'], 'pr');
    writetable(struct2table(pr, 'AsArray', true), [imsPath, 'parameters.txt'])
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
    if ispc
        rename_str = sprintf('powershell -command "mv "%s" "%s""', tmpFullpaths{d}, outputFullpaths{d});
    else
        rename_str = sprintf('mv "%s" "%s"', tmpFullpaths{d}, outputFullpaths{d});
    end

    func_strs{d} = sprintf('"%s" -P "%s" -F "%s" -n "%s" -r %s -v %s -o "%s" -b %s %s && %s', converterPath, ...
        channelPatterns_str, dataPath_str, tmpFsns{d}, type_str, voxelsize_str, imsPaths{d}, ...
        blockSize_str, bbox_str, rename_str);
end

% use slurm job wrapper for computing
generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, func_strs, ...
    'cpusPerTask', cpusPerTask, 'parseCluster', parseCluster, 'language', 'bash', ...
    mccMode=mccMode, configFile=configFile);

end

