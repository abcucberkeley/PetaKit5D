function [] = XR_crop_dataset(dataPaths, inputBbox, varargin)
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
ip.addRequired('inputBbox', @isnumeric);
ip.addParameter('resultDirName', 'Cropped', @ischar);
ip.addParameter('cropType', 'fixed', @ischar); % fixed or moving or center
ip.addParameter('pad', false, @islogical); % pad region that is outside the bbox
ip.addParameter('lastStartCoords', [], @isnumeric); % start coordinate of the last time point
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamB_ch0'}, @iscell);
ip.addParameter('zarrFile', false , @islogical); % read zarr
ip.addParameter('largeZarr', false, @islogical); % use zarr file as input
ip.addParameter('saveZarr', false , @islogical); % save as zarr
ip.addParameter('blockSize', [500, 500, 500] , @isnumeric);
ip.addParameter('save16bit', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
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
largeZarr = pr.largeZarr;
saveZarr = pr.saveZarr;
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
        '''largeZarr'',%s,''saveZarr'',%s,''BlockSize'',%s,''uuid'',''%s'',''parseCluster'',%s,', ...
        '''mccMode'',%s,''configFile'',''%s'')'], frameFullpath, cropFullpath, ...
        strrep(num2str(bbox_f, '%.10f,'), ' ', ''), string(pad), string(zarrFile), ...
        string(largeZarr), string(saveZarr), strrep(mat2str(blockSize), ' ', ','), ...
        uuid, string(parseCluster), string(mccMode), configFile);
end

memAllocate = prod(sz) * 4 / 1024^3 * 8;
generic_computing_frameworks_wrapper(frameFullpaths, cropFullpaths, func_strs, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, memAllocate=memAllocate, mccMode=mccMode, configFile=configFile);

end
