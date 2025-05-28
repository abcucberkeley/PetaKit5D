function [] = XR_tiffToZarr_wrapper(dataPaths, varargin)
% Dataset level wrapper for convert a list of Tiff files to Zarr.
%
%
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%
% Parameters (as 'specifier'-value pairs):
%       tiffFullpaths : empty or char or cell array (default: ''). A single Tiff file path in char, or a list of Tiff file paths in a cell array. If empty, get the Tiff files from dataPaths.
%       resultDirName : char (default: 'zarr'). Result directory under data paths.
%              locIds : empty or 1x#Tile vector (default: []). Location ids for the tiles (only used in stitching).
%           blockSize : 1x3 vector (default: [256, 256, 256]). Block/chunk size for zarr output.
%           shardSize : empty or 1x3 vector (default: []). If not empty, set shard size in zarr output. 
%         flippedTile : empty or 1x#Tile vector (default: []). Define if a tile is flipped, if so, flip the flipped tiles back. If empty, tiles are not flipped.
%      resampleFactor : empty or 1x1, 1x2 or 1x3 vector (default: []). Resampling factor after rotation. Empty: no resampling; axis order yxz.
%         partialFile : true|false (default: false). If true, check if containing partial files ending with *_part[0-9][0-9][0-9].tif. Append partial files to the main ones if there are.
%     channelPatterns : a cell array (default: {'tif'}).  Channel identifiers for included channels. 
%           inputBbox : empty or 1x6 vector (default: []). Input bounding box for crop. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%         tileOutBbox : empty or 1x6 vector (default: []). Crop tiles after preprocessing. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%       readWholeTiff : true|false (default: true). Read the entire image for the conversion. Set it to false when the image is too large to fit to the RAM.
%      processFunPath : empty or char (default: ''). Path for user-defined process function handle. Support .mat or .txt formats.
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%             bigData : true|false (default: true). Big dataset, use fewer resources for each task to all more parallel task running.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%           jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%         cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%         maxTrialNum : a number (default: 3). The max number of retries for a task.
%        unitWaitTime : a number (default: 1). The wait time per file in minutes to check whether the computing is done.
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%
%
% Author: Xiongtao Ruan (10/02/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) iscell(x) || ischar(x));
ip.addParameter('tiffFullpaths', '', @(x) iscell(x) || ischar(x));
ip.addParameter('resultDirName', 'zarr', @ischar);
ip.addParameter('locIds', [], @isnumeric);
ip.addParameter('blockSize', [256, 256, 256], @isnumeric);
ip.addParameter('shardSize', [], @isnumeric);
ip.addParameter('flippedTile', [], @(x) isempty(x) || islogical(x));
ip.addParameter('resampleFactor', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('partialFile', false, @islogical);
ip.addParameter('channelPatterns', {'tif'}, @iscell);
ip.addParameter('inputBbox', [], @isnumeric);
ip.addParameter('tileOutBbox', [], @isnumeric);
ip.addParameter('readWholeTiff', true, @islogical);
ip.addParameter('processFunPath', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x) || isstring(x) || iscell(x));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('bigData', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 3, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
 % Resolution = pr.Resolution;
tiffFullpaths = pr.tiffFullpaths;
resultDirName = pr.resultDirName;
locIds = pr.locIds;
blockSize = pr.blockSize;
shardSize= pr.shardSize;
flippedTile = pr.flippedTile;
resampleFactor = pr.resampleFactor;
partialFile = pr.partialFile;
channelPatterns = pr.channelPatterns;
inputBbox = pr.inputBbox;
tileOutBbox = pr.tileOutBbox;
readWholeTiff = pr.readWholeTiff;
processFunPath = pr.processFunPath;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
bigData = pr.bigData;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ~isempty(dataPaths)
    if ischar(dataPaths)
        dataPaths = {dataPaths};
    end
    
    nd = numel(dataPaths);
    resultPaths = cell(nd, 1);
    for d = 1 : nd
        dataPaths{d} = simplifyPath(dataPaths{d});
        resultPath = [dataPaths{d}, '/' resultDirName, '/'];
        resultPaths{d} = resultPath;
        mkdir(resultPath);
        fileattrib(resultPath, '+w', 'g');
    
        save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
        writeJsonFile(pr, [resultPath, '/parameters.json']);
    end
    
    % parse image filenames
    zarrFile = false;
    [~, ~, ~, filepaths] = parseImageFilenames(dataPaths, zarrFile, channelPatterns);
    tiffFullpaths = filepaths;
end

if ischar(tiffFullpaths)
    if exist(tiffFullpaths, 'dir')
        dir_info = dir([tiffFullpaths, '/', '*.tif']);
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

nC = numel(channelPatterns);
usrFcn_strs = repmat({''}, nC, max(1, size(processFunPath, 2)));
keywords = repmat({{''}}, nC, max(1, size(processFunPath, 2)));
if ~isempty(processFunPath)
    fprintf('Process function paths:\n')
    disp(processFunPath');
end
if ~isempty(processFunPath)
    for i = 1 : size(processFunPath, 1)
        for j = 1 : size(processFunPath, 2)
            if isempty(processFunPath{i, j})
                continue;
            end
            [~, ~, ext] = fileparts(processFunPath{i, j});
            switch ext
                case '.mat'
                    a = load(processFunPath{i, j});
                    usrFun = a.usrFun;
                case '.txt'
                    a = readTextFile(processFunPath{i, j});
                    usrFun = a{1};
                otherwise
                    error('Unknown file type for the processed function file!')
            end
            if ~isempty(usrFun)
                usrFcn_strs{i, j} = usrFun;
            end
            if isfield(a, 'keywords')
                keywords{i, j} = a.keywords;
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
        fnames_cell{d} = cellfun(@(x) [udataPaths{d}, '/', x], fsnames, 'unif', 0);
    end
end

% if bigData
%     compressor = 'zstd';
% else
%     compressor = 'lz4';
% end
compressor = 'zstd';

uniq_locIds = unique(locIds);
nLoc = numel(uniq_locIds);
locSpecific = false;
if numel(uniq_locIds) > 1 && (size(inputBbox, 1) == nLoc || size(tileOutBbox, 1) == nLoc)
    locSpecific = true;
    if size(inputBbox, 1) <= 1
        inputBbox = repmat({inputBbox}, nLoc, 1);
    elseif size(inputBbox, 1) == nLoc
        inputBbox = mat2cell(inputBbox, ones(1, size(inputBbox, 1)), size(inputBbox, 2));
    else
        error('The number of InputBbox (%d) does not match the number of tile locations (%d)!', size(inputBbox, 1), numel(uniq_locIds));
    end
    if size(tileOutBbox, 1) <= 1
        tileOutBbox = repmat({tileOutBbox}, nLoc, 1);
    elseif size(tileOutBbox, 1) == nLoc
        tileOutBbox = mat2cell(tileOutBbox, ones(1, size(tileOutBbox, 1)), size(tileOutBbox, 2));
    else
        error('The number of tileOutBbox (%d) does not match the number of tile locations (%d)!', size(tileOutBbox, 1), numel(uniq_locIds));
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

    zarrPath = [dataPath, '/', resultDirName];
    if ~exist(zarrPath, 'dir')
        mkdir(zarrPath);
        fileattrib(zarrPath, '+w', 'g');
    end

    zarrFullpaths{i} = sprintf('%s/%s/%s.zarr/', dataPath, resultDirName, fsname_i);
    if ~isempty(flippedTile)
        flipZstack = flippedTile(i);
    else
        flipZstack = false;
    end
    
    cind = cellfun(@(x) contains(tiffFullpath_i, x), channelPatterns);
    if ~any(cind)
        error('The file %s does not match any channel patterns %s', tiffFullpath_i, string(channelPatterns));
    end
    
    if locSpecific
        locInd = uniq_locIds == locIds(i);
        inputBbox_i = inputBbox{locInd};
        tileOutBbox_i = tileOutBbox{locInd};
        % usrFcn_str_i = usrFcn_strs{cind, i, locInd};        
    else
        inputBbox_i = inputBbox;
        tileOutBbox_i = tileOutBbox;
        % usrFcn_str_i = usrFcn_strs{cind, i};
    end

    % map user defined function to the tile handle based on the keywords
    usrFcn_strs_c = usrFcn_strs(cind, :);
    keywords_c = keywords(cind, :);
    usrFcn_str_i = usrFcn_strs_c{1};
    for j = 1 : numel(usrFcn_strs_c)
        break_loops = false;
        keywords_cj = keywords_c{j};
        for k = 1 : numel(keywords_cj)
            if contains(tiffFullpath_i, keywords_cj{k})
                break_loops = true;
                break;
            end
        end
        if break_loops
            usrFcn_str_i = usrFcn_strs_c{j};
            break;
        end
    end

    func_strs{i} = sprintf(['tiffToZarr(%s,''%s'',[],''BlockSize'',%s,''shardSize'',%s,', ...
        '''flipZstack'',%s,''resampleFactor'',%s,''inputBbox'',%s,''tileOutBbox'',%s,', ...
        '''readWholeTiff'',%s,''compressor'',''%s'',''usrFcn'',"%s")'], sprintf('{''%s''}', ...
        strjoin(tiffFullpath_group_i, ''',''')), zarrFullpaths{i}, strrep(mat2str(blockSize), ' ', ','), ...
        strrep(mat2str(shardSize), ' ', ','), string(flipZstack), strrep(mat2str(resampleFactor), ' ', ','), ...
        strrep(mat2str(inputBbox_i), ' ', ','), strrep(mat2str(tileOutBbox_i), ' ', ','), ...
        string(readWholeTiff), compressor, usrFcn_str_i);
end

imSizes = zeros(numel(tiffFullpath_group_i), 3);
for i = 1 : numel(tiffFullpath_group_i)
    imSizes(i, :) = getImageSize(tiffFullpath_group_i{i});
end
dtype = getImageDataType(tiffFullpath_group_i{1});
byteNum = dataTypeToByteNumber(dtype);

imSize = [imSizes(1, 1 : 2), sum(imSizes(:, 3))];
if all(isempty(usrFcn_strs))
    memFactor = 2.5;
else
    % for our own flat field code, set memFactor as 2.75, and other user defined function, set memFactor as 4
    if contains(usrFcn_strs{1}, 'processFFCorrectionFrame')
        memFactor = 2.75;
    else
        memFactor = 4;
    end
end
memAllocate = prod(imSize) * byteNum / 1024^3 * memFactor;
if ~bigData
    memAllocate = memAllocate * 2;
end
maxTrialNum = 2;

is_done_flag = generic_computing_frameworks_wrapper(tiffFullpaths, zarrFullpaths, ...
    func_strs, 'parseCluster', parseCluster, 'masterCompute', masterCompute, ...
    'maxTrialNum', maxTrialNum,  'cpusPerTask', cpusPerTask, 'memAllocate', memAllocate, ...
    'mccMode', mccMode, 'configFile', configFile);
if ~all(is_done_flag)
    generic_computing_frameworks_wrapper(tiffFullpaths, zarrFullpaths, func_strs, ...
        'parseCluster', parseCluster, 'masterCompute', masterCompute, 'maxTrialNum', maxTrialNum, ...
        'cpusPerTask', cpusPerTask * 2, 'memAllocate', memAllocate * 2, 'mccMode', mccMode, ...
        'configFile', configFile);
end

end

