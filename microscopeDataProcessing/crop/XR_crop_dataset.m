function [] = XR_crop_dataset(dataPaths, inputBbox, varargin)
% Dataset level wrapper for the cropping tool with a given cropping bounding box.
% 
% 
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%           inputBbox : 1x6 vector. Input bounding box for crop. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%
% Parameters (as 'specifier'-value pairs):
%       resultDirName : char (default: 'Cropped'). Crop result directory under data path.
%            cropType : 'fixed'|'moving'|'center' (default: 'fixed'). 'fixed': fixed bounding box; 'moving': crop the first time point with inputBbox, and linearly move the bbox according to the last time point coordinates (lastStartCoords); 'center': crop the center of the images with the size defined by inputBbox.
%                 pad : true|false (default: false). Pad an empty region if the bounding box is outside of the image.
%     lastStartCoords : empty or 1x3 vector. The start coordinates for the last time point for 'moving' crop type.
%     channelPatterns : a cell array (default: {'CamA_ch0', 'CamB_ch0'}).  Channel identifiers for included channels. 
%            zarrFile : true|false (default: false). Use Zarr file as input.
%           largeFile : true|false (default: false). Use large-scale crop strategy with in-place read and write cropped regions. Only for Zarr files.
%            saveZarr : true|false (default: false). Save results as Zarr files.
%           batchSize : 1x3 vector (default: [1024, 1024, 1024]). Batch size per task for large scale cropping.
%           blockSize : 1x3 vector (default: [256, 256, 256]). Block/chunk size for zarr output.
%           save16bit : true|false (default: true). Save 16bit result for deskew/rotate and stitch. 
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%           jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%         cpusPerTask : a number (default: 2). The number of cpu cores per task for slurm job submission.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%
%
% Author: Xiongtao Ruan (02/18/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('inputBbox', @isnumeric);
ip.addParameter('resultDirName', 'Cropped', @ischar);
ip.addParameter('cropType', 'fixed', @ischar);
ip.addParameter('pad', false, @islogical);
ip.addParameter('lastStartCoords', [], @isnumeric);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamB_ch0'}, @iscell);
ip.addParameter('zarrFile', false , @islogical);
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('saveZarr', false , @islogical);
ip.addParameter('batchSize', [1024, 1024, 1024] , @isnumeric);
ip.addParameter('blockSize', [256, 256, 256] , @isnumeric);
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, inputBbox, varargin{:});

warning('off', 'MATLAB:MKDIR:DirectoryExists');

pr = ip.Results;
% Overwrite = pr.Overwrite;
resultDirName = pr.resultDirName;
cropType = pr.cropType;
pad = pr.pad;
lastStartCoords = pr.lastStartCoords;
channelPatterns = pr.channelPatterns;
zarrFile = pr.zarrFile;
largeFile = pr.largeFile;
saveZarr = pr.saveZarr;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if isempty(uuid)
    uuid = get_uuid();
end

% temporary directory for intermediate results
if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    resultPath = sprintf('%s/%s/', dataPath, resultDirName);
    if ~exist(resultPath, 'dir')
        mkdir(resultPath);
    end
    
    dataPaths{d} = simplifyPath(dataPath);
    resultPaths{d} = simplifyPath(resultPath);
    try
        fileattrib(resultPath, '+w', 'g');
    catch ME
        disp(ME);
    end
    save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
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
    nc = numel(channelPatterns);
    ch_inds = false(nF, nc);
    for c = 1 : nc
        ch_inds(:, c) = contains(fnames_d, channelPatterns{c}, 'IgnoreCase', true) | contains(fnames_d, regexpPattern(channelPatterns{c}), 'IgnoreCase', true);
    end
    fnames_d = fnames_d(any(ch_inds, 2));
    fnames{d} = fnames_d;
    
    if strcmp(cropType, 'moving')
        ch_inds = ch_inds(any(ch_inds, 2), :);
        ch_tps = cumsum(ch_inds);
        ch_tps = sum(ch_tps .* ch_inds, 2);
        max_tp = max(ch_tps);
        max_tp_mat(d) = max_tp;
        if any(inputBbox(4 : 6) - inputBbox(1 : 3) + lastStartCoords > sz)
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
        bbox_size = inputBbox(4 : 6) - inputBbox(1 : 3) + 1;
        bbox_start_f = inputBbox(1 : 3) + round((lastStartCoords - inputBbox(1 : 3)) / (max_tp - 1) * (ch_tps(f) - 1));
        bbox_f = [bbox_start_f, bbox_start_f + bbox_size - 1];
    elseif strcmp(cropType, 'center')
        bbox_size = inputBbox(4 : 6) - inputBbox(1 : 3) + 1;
        % half size
        hsize = floor(bbox_size / 2);
        imSize = getImageSize(frameFullpath);
        % image center
        center = round((imSize + 1) / 2);
        bbox_start_f = center - hsize;
        bbox_f = [bbox_start_f, bbox_start_f + bbox_size - 1];        
    else
        bbox_f = inputBbox;
    end
    
    func_strs{f} = sprintf(['XR_crop_frame(''%s'',''%s'',[%s],''pad'',%s,''zarrFile'',%s,', ...
        '''largeFile'',%s,''saveZarr'',%s,''batchSize'',%s,''blockSize'',%s,''uuid'',''%s'',''parseCluster'',%s,', ...
        '''mccMode'',%s,''configFile'',''%s'')'], frameFullpath, cropFullpath, ...
        strrep(num2str(bbox_f, '%.10f,'), ' ', ''), string(pad), string(zarrFile), ...
        string(largeFile), string(saveZarr), mat2str_comma(batchSize), strrep(mat2str(blockSize), ' ', ','), ...
        uuid, string(parseCluster), string(mccMode), configFile);
end

memAllocate = prod(sz) * 4 / 1024^3 * 8;
generic_computing_frameworks_wrapper(frameFullpaths, cropFullpaths, func_strs, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, memAllocate=memAllocate, mccMode=mccMode, configFile=configFile);

end
