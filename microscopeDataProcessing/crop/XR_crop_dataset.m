function [] = XR_crop_dataset(dataPaths, resultPaths, bbox, varargin)
% crop a dataset via cluster computing
% 
%
% Author: Xiongtao Ruan (02/18/2020)
% 
% xruan (08/22/2020): add group write permission for the result folder
% xruan (08/23/2020): add option for not using tmp directory if write to a new folder.
% xruan (01/15/2021): refactor code and add option for moving bounding box.
% xruan (06/25/2021): add support for crop from center, the size is from bbox
% xruan (07/13/2021): add option to pad data if it is outside of the bbox
% xruan (07/13/2021): add support for multiple datasets
% xruan (01/25/2022): add support for zarr read and write
% xruan (06/03/2022): add support for large zarr files


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('resultPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('bbox', @isnumeric);
ip.addParameter('cropType', 'fixed', @ischar); % fixed or moving or center
ip.addParameter('pad', false, @islogical); % pad region that is outside the bbox
ip.addParameter('lastStart', [], @isnumeric); % start coordinate of the last time point
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamB_ch0'}, @iscell);
ip.addParameter('zarrFile', false , @islogical); % read zarr
ip.addParameter('largeZarr', false, @islogical); % use zarr file as input
ip.addParameter('saveZarr', false , @islogical); % save as zarr
ip.addParameter('BlockSize', [500, 500, 500] , @isnumeric); % save as zarr
ip.addParameter('Save16bit', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPaths, resultPaths, bbox, varargin{:});

warning('off', 'MATLAB:MKDIR:DirectoryExists');

pr = ip.Results;
% Overwrite = pr.Overwrite;
cropType = pr.cropType;
pad = pr.pad;
lastStart = pr.lastStart;
ChannelPatterns = pr.ChannelPatterns;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
saveZarr = pr.saveZarr;
BlockSize = pr.BlockSize;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes = pr.cpuOnlyNodes;

% temporary directory for intermediate results
if ischar(dataPaths)
    dataPaths = {dataPaths};
    resultPaths = {resultPaths};
end

nd = numel(dataPaths);
for d = 1 : nd
    dataPath = dataPaths{d};
    resultPath = resultPaths{d};
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

% processing file paths and bbox for moving option. 
fnames = cell(nd, 1);
if strcmp(cropType, 'moving')
    max_tp_mat = zeros(nd, 1);
end

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
    fnames_d = sort(fnames_d); % filenames are in order of time points
    nc = numel(ChannelPatterns);
    ch_inds = false(nF, nc);
    for c = 1 : nc
        ch_inds(:, c) = contains(fnames_d, ChannelPatterns{c}, 'IgnoreCase', true) | contains(fnames_d, regexpPattern(ChannelPatterns{c}), 'IgnoreCase', true);
    end
    fnames_d = fnames_d(any(ch_inds, 2));
    fnames{d} = fnames_d;
    
    if strcmp(cropType, 'moving')
        ch_inds = ch_inds(any(ch_inds, 2), :);
        ch_tps = cumsum(ch_inds);
        ch_tps = sum(ch_tps .* ch_inds, 2);
        max_tp = max(ch_tps);
        max_tp_mat(d) = max_tp;
        if any(bbox(4 : 6) - bbox(1 : 3) + lastStart > sz)
            error('The moving bounding box moves out of the image for the last time point!')
        end
    end
end

fd_inds = arrayfun(@(x) ones(numel(fnames{x}), 1) * x, 1 : nd, 'unif', 0);
fnames = cat(1, fnames{:});
fd_inds = cat(1, fd_inds{:});
nF = numel(fnames);


%% use generic framework for the cropping do computing
% first define the cropping function strings
frameFullpaths = cell(nF, 1);
cropFullpaths = cell(nF, 1);
func_strs = cell(nF, 1);
for f = 1 : nF
    fname = fnames{f};
    d = fd_inds(f);
    dataPath = dataPaths{d};
    resultPath = resultPaths{d};
    
    frameFullpath = [dataPath, '/', fname];
    frameFullpaths{f} = frameFullpath;
    cropFullpath = [resultPath, '/', fname];
    if zarrFile && ~saveZarr
        cropFullpath = [resultPath, '/', fname(1 : end - 5), '.tif'];
    end
    cropFullpaths{f} = cropFullpath;

    if strcmp(cropType, 'moving')
        max_tp = max_tp_mat(d);
        bbox_size = bbox(4 : 6) - bbox(1 : 3) + 1;
        bbox_start_f = bbox(1 : 3) + round((lastStart - bbox(1 : 3)) / (max_tp - 1) * (ch_tps(f) - 1));
        bbox_f = [bbox_start_f, bbox_start_f + bbox_size - 1];
    elseif strcmp(cropType, 'center')
        bbox_size = bbox(4 : 6) - bbox(1 : 3) + 1;
        % half size
        hsize = floor(bbox_size / 2);
        imSize = getImageSize(frameFullpath);
        % image center
        center = round((imSize + 1) / 2);
        bbox_start_f = center - hsize;
        bbox_f = [bbox_start_f, bbox_start_f + bbox_size - 1];        
    else
        bbox_f = bbox;
    end
    
    func_strs{f} = sprintf(['XR_crop_frame(''%s'',''%s'',[%s],''pad'',%s,''zarrFile'',%s,', ...
        '''largeZarr'',%s,''saveZarr'',%s,''BlockSize'',%s)'], frameFullpath, cropFullpath, ...
        strrep(num2str(bbox_f, '%.10f,'), ' ', ''), string(pad), string(zarrFile), ...
        string(largeZarr), string(saveZarr), strrep(mat2str(BlockSize), ' ', ','));
end

cpusPerTask = max(min(ceil(prod(sz) * 4 / 1024^3 * 8 / 20), 24), cpusPerTask);
slurm_cluster_generic_computing_wrapper(frameFullpaths, cropFullpaths, func_strs, ...
    'masterCompute', masterCompute, 'cpusPerTask', cpusPerTask);

end
