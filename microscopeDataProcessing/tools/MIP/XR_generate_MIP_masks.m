function [] = XR_generate_MIP_masks(dataPaths, varargin)
% generate MIP masks based on intensity thresholding and object size. The
% main purpose is to only keep the major objects in the MIPs (for big data) for
% post-processing, i.e., deskew and rotate, deconvolution. 
%
% Author: Xiongtao Ruan (08/02/2023)



ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('ResultDirStr', 'masks', @ischar); % y, x, z
ip.addParameter('axis', [1, 1, 1], @isnumeric); % y, x, z
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @iscell);
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('saveZarr', false, @islogical); % use zarr file as output
ip.addParameter('blockSize', [256, 256, 1], @isnumeric); % zarr output block size
ip.addParameter('intThresh', [100, 100, 100], @isnumeric); % intensity threshold, y, x, z
ip.addParameter('volThresh', 100, @isnumeric); % volume threshold
ip.addParameter('dilateSize', 100, @isnumeric); % dilate the mask to add some buffer room
ip.addParameter('Save16bit', true, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 3, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
ResultDirStr = pr.ResultDirStr;
axis = pr.axis;
ChannelPatterns =  pr.ChannelPatterns;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
intThresh = pr.intThresh;
volThresh = pr.volThresh;
dilateSize = pr.dilateSize;
Save16bit = pr.Save16bit;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;

parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    resultPath = [dataPath, '/', ResultDirStr, '/'];
    mkdir(resultPath);
    
    dataPaths{d} = simplifyPath(dataPath);
    resultPaths{d} = simplifyPath(resultPath);    
    fileattrib(resultPath, '+w', 'g');
    save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
end

fnames = cell(nd, 1);
mip_axes = cell(nd, 1);
axis_strs = {'y', 'x', 'z'};

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
    
    mip_axes_d = zeros(nF, 1);
    nF = numel(fnames_d);
    axis_inds = false(nF, 3);
    for i = 1 : 3
        axis_pattern = sprintf('MIP_%s', axis_strs{i});
        axis_inds(:, i) = contains(fnames_d, axis_pattern, 'IgnoreCase', true);
        mip_axes_d(axis_inds(:, i)) = i;

        if axis(i) > 0
            if ~any(axis_inds(:, i))
                warning('MIPs for axis %s for data %s do not exist!', axis_strs{i}, dataPath);
            end
        else
            axis_inds(:, i) = false;
        end
    end

    fnames_d = fnames_d(any(axis_inds, 2));
    
    fnames{d} = fnames_d;
    mip_axes{d} = mip_axes_d;
end

fd_inds = arrayfun(@(x) ones(numel(fnames{x}), 1) * x, 1 : nd, 'unif', 0);
fnames = cat(1, fnames{:});
[~, fsns] = fileparts(fnames);
fd_inds = cat(1, fd_inds{:});
mip_axes = cat(1, mip_axes{:});
nF = numel(fnames);
if nF == 1
    fsns = {fsns};
end

%% use generic framework for the MIP mask computing

if numel(intThresh) == 1
    intThresh = repmat(intThresh, 1, 3);
end

if numel(dilateSize) == 1
    dilateSize = repmat(dilateSize, 1, 3);
end

frameFullpaths = cell(nF, 1);
outFullpaths = cell(nF, 1);
func_strs = cell(nF, 1);

for f = 1 : nF
    fname = fnames{f};
    d = fd_inds(f);
    dataPath = dataPaths{d};
    resultPath = resultPaths{d};
    
    fn = [dataPath, '/', fname];
    frameFullpaths{f} = fn;
    if saveZarr
        fnout = sprintf('%s/%s.zarr', resultPath, fsns{f});        
    else
        fnout = sprintf('%s/%s.tif', resultPath, fsns{f});
    end
    outFullpaths{f} = fnout;
    
    dilateSize_f = dilateSize;
    dilateSize_f(mip_axes(f)) = [];    
    
    intThresh_f = intThresh(mip_axes(f));

    func_strs{f} = sprintf(['gererate_single_MIP_mask(''%s'',''%s'',%d,%d,%s,%s,%s,''%s'')'], ...
        fn, fnout, intThresh_f, volThresh, strrep(mat2str(dilateSize_f), ' ', ','), strrep(mat2str(blockSize), ' ', ','), ...
        string(Save16bit), uuid);
end

taskBatchNum = max(1, round(nF /1000));
memAllocate = prod(sz) * 4 / 1024^3 * 20;
generic_computing_frameworks_wrapper(frameFullpaths, outFullpaths, func_strs, ...
    parseCluster=parseCluster, masterCompute=masterCompute, taskBatchNum=taskBatchNum, ...
    cpusPerTask=cpusPerTask, memAllocate=memAllocate, mccMode=mccMode, ConfigFile=ConfigFile);


end

