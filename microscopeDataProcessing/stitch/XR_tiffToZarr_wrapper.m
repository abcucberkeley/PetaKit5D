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
% xruan (08/25/2022): change CropToSize to tileOutBbox (more generic)
% xruan (05/26/2023): add support for loc specific pocessing for multiLoc
% datasets, support InputBbox and tileOutBbox for now

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tiffFullpaths', @(x) iscell(x) || ischar(x));
ip.addParameter('zarrPathstr', 'zarr', @ischar);
ip.addParameter('locIds', [], @isnumeric); % location ids for the tiles
ip.addParameter('blockSize', [500, 500, 250], @isnumeric);
ip.addParameter('flippedTile', [], @(x) isempty(x) || islogical(x));
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('partialFile', false, @islogical);
ip.addParameter('ChannelPatterns', {'tif'}, @iscell);
ip.addParameter('InputBbox', [], @isnumeric); % crop input tile before processing
ip.addParameter('tileOutBbox', [], @isnumeric); % crop output tile after processing
ip.addParameter('processFunPath', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x) || isstring(x) || iscell(x));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('bigData', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);


ip.parse(tiffFullpaths, varargin{:});

pr = ip.Results;
 % Resolution = pr.Resolution;
zarrPathstr = pr.zarrPathstr;
locIds = pr.locIds;
blockSize = pr.blockSize;
flippedTile = pr.flippedTile;
resample = pr.resample;
partialFile = pr.partialFile;
ChannelPatterns = pr.ChannelPatterns;
InputBbox = pr.InputBbox;
tileOutBbox = pr.tileOutBbox;
processFunPath = pr.processFunPath;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
bigData = pr.bigData;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

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
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath, jobLogDir);
end

nC = numel(ChannelPatterns);
usrFcn_strs = repmat({''}, nC, size(processFunPath, 2));
fprintf('Process function paths:\n')
disp(processFunPath');
if ~isempty(processFunPath)
    for i = 1 : size(processFunPath, 1)
        for j = 1 : size(processFunPath, 2)
            if isempty(processFunPath{i, j})
                continue;
            end
            a = load(processFunPath{i, j});
            usrFun = a.usrFun;
            if ~isempty(usrFun)
                usrFcn_strs{i, j} = usrFun;
            end
        end
    end
end

nF = numel(tiffFullpaths);
if partialFile
    dataPaths = fileparts(tiffFullpaths);
    udataPaths = unique(dataPaths);
    fnames_cell = cell(numel(udataPaths), 1);
    for d = 1 : numel(udataPaths)
        dir_info = dir([udataPaths{d}, '/*tif']);
        fsnames = {dir_info.name}';
        fnames_cell{d} = cellfun(@(x) [udataPaths{d}, filesep, x], fsnames, 'unif', 0);
    end
end

if bigData
    compressor = 'zstd';
else
    compressor = 'lz4';
end

uniq_locIds = unique(locIds);
nLoc = numel(uniq_locIds);
locSpecific = false;
if numel(uniq_locIds) > 1 && (size(InputBbox, 1) == nLoc || size(tileOutBbox, 1) == nLoc || ...
        size(usrFcn_strs, 2) == nLoc)
    locSpecific = true;
    if size(InputBbox, 1) <= 1
        InputBbox = repmat({InputBbox}, nLoc, 1);
    elseif size(InputBbox, 1) == nLoc
        InputBbox = mat2cell(InputBbox, ones(1, size(InputBbox, 1)), size(InputBbox, 2));
    else
        error('The number of InputBbox (%d) does not match the number of tile locations (%d)!', size(InputBbox, 1), numel(uniq_locIds));
    end
    if size(tileOutBbox, 1) <= 1
        tileOutBbox = repmat({tileOutBbox}, nLoc, 1);
    elseif size(tileOutBbox, 1) == nLoc
        tileOutBbox = mat2cell(tileOutBbox, ones(1, size(tileOutBbox, 1)), size(tileOutBbox, 2));
    else
        error('The number of tileOutBbox (%d) does not match the number of tile locations (%d)!', size(tileOutBbox, 1), numel(uniq_locIds));
    end
    if size(usrFcn_strs, 2) <= 1
        usrFcn_strs = repmat(usrFcn_strs, 1, nLoc);
    elseif size(usrFcn_strs, 2) == nLoc
    else
        error('The number of usrFcn (%d) does not match the number of tile locations (%d)!', size(usrFcn_strs, 2), numel(uniq_locIds));
    end
end

zarrFullpaths = cell(nF, 1);
func_strs = cell(nF, 1);
for i = 1 : nF
    tiffFullpath_i = tiffFullpaths{i};
    [dataPath, fsname_i] = fileparts(tiffFullpath_i);
    
    if partialFile
        if numel(udataPaths) > 1
            did = ismember(udataPaths, dataPath);
        else
            did = 1;
        end
        fnames = fnames_cell{did};
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
    
    if locSpecific
        locInd = uniq_locIds == locIds(i);
        InputBbox_i = InputBbox{locInd};
        tileOutBbox_i = tileOutBbox{locInd};
        usrFcn_str_i = usrFcn_strs{cind, locInd};        
    else
        InputBbox_i = InputBbox;
        tileOutBbox_i = tileOutBbox;
        usrFcn_str_i = usrFcn_strs{cind, 1};
    end
    
    func_strs{i} = sprintf(['tiffToZarr(%s,''%s'',[],''BlockSize'',%s,''flipZstack'',%s,', ...
        '''resample'',%s,''InputBbox'',%s,''tileOutBbox'',%s,''compressor'',''%s'',''usrFcn'',"%s")'], ...
        sprintf('{''%s''}', strjoin(tiffFullpath_group_i, ''',''')), zarrFullpaths{i}, ...
        strrep(mat2str(blockSize), ' ', ','), string(flipZstack), strrep(mat2str(resample), ' ', ','), ...
        strrep(mat2str(InputBbox_i), ' ', ','), strrep(mat2str(tileOutBbox_i), ' ', ','), ...
        compressor, usrFcn_str_i);
end

imSizes = zeros(numel(tiffFullpath_group_i), 3);
for i = 1 : numel(tiffFullpath_group_i)
    imSizes(i, :) = getImageSize(tiffFullpath_group_i{i});
end
imSize = [imSizes(1, 1 : 2), sum(imSizes(:, 3))];
memAllocate = prod(imSize) * 4 / 1024^3 * 2.5;
if ~bigData
    memAllocate = memAllocate * 2;
end
maxTrialNum = 2;

is_done_flag = generic_computing_frameworks_wrapper(tiffFullpaths, zarrFullpaths, ...
    func_strs, 'parseCluster', parseCluster, 'masterCompute', masterCompute, ...
    'maxTrialNum', maxTrialNum,  'cpusPerTask', cpusPerTask, 'memAllocate', memAllocate, ...
    'mccMode', mccMode, 'ConfigFile', ConfigFile);
if ~all(is_done_flag)
    generic_computing_frameworks_wrapper(tiffFullpaths, zarrFullpaths, func_strs, ...
        'parseCluster', parseCluster, 'masterCompute', masterCompute, 'maxTrialNum', maxTrialNum, ...
        'cpusPerTask', cpusPerTask * 2, 'memAllocate', memAllocate * 2, 'mccMode', mccMode, ...
        'ConfigFile', ConfigFile);
end


end


