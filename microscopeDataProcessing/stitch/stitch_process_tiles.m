function [] = stitch_process_tiles(inputFullpaths, varargin)
% The wrapper for convert a list of tiff files (a data folder) to zarr or
% process zarr tiles
% 
% copied from XR_tiffToZarr_wrapper.m
%
% Author: Xiongtao Ruan (04/29/2023)
% 
% xruan (05/26/2023): add support for loc specific pocessing for multiLoc
% datasets, support InputBbox and tileOutBbox for now

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpaths', @(x) iscell(x) || ischar(x));
ip.addParameter('resultDirName', 'zarr', @ischar);
ip.addParameter('zarrFile', false, @islogical); 
ip.addParameter('locIds', [], @isnumeric); % location ids for the tiles
ip.addParameter('blockSize', [500, 500, 250], @isnumeric);
ip.addParameter('shardSize', [], @isnumeric);
ip.addParameter('flippedTile', [], @(x) isempty(x) || islogical(x));
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('partialFile', false, @islogical);
ip.addParameter('channelPatterns', {'.'}, @iscell);
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
ip.addParameter('configFile', '', @ischar);

ip.parse(inputFullpaths, varargin{:});

pr = ip.Results;
 % Resolution = pr.Resolution;
resultDirName = pr.resultDirName;
zarrFile = pr.zarrFile;
locIds = pr.locIds;
blockSize = pr.blockSize;
shardSize = pr.shardSize;
flippedTile = pr.flippedTile;
resample = pr.resample;
partialFile = pr.partialFile;
channelPatterns = pr.channelPatterns;
InputBbox = pr.InputBbox;
tileOutBbox = pr.tileOutBbox;
processFunPath = pr.processFunPath;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
bigData = pr.bigData;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ~zarrFile
    fprintf('\nConvert tiles from Tiff to Zarr...\n');
    XR_tiffToZarr_wrapper('', 'tiffFullpaths', inputFullpaths, 'resultDirName', resultDirName, ...
        'locIds', locIds, 'blockSize', blockSize, 'shardSize', shardSize, 'flippedTile', flippedTile, ...
        'resample', resample, 'partialFile', partialFile, 'channelPatterns', channelPatterns, ...
        'InputBbox', InputBbox, 'tileOutBbox', tileOutBbox, 'processFunPath', processFunPath, ...
        'jobLogDir', jobLogDir, 'parseCluster', parseCluster, 'masterCompute', masterCompute, ...
        'bigData', bigData, 'cpusPerTask', cpusPerTask, 'mccMode', mccMode, 'configFile', configFile);
    return;
end

if isempty(resultDirName)
    return;
end

fprintf('\nPreprocess Zarr tiles...\n');

[dataPath, ~] = fileparts(inputFullpaths{1});
    
% check if a slurm-based computing cluster exist
if parseCluster 
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath, jobLogDir);
end

nC = numel(channelPatterns);
usrFcn_strs = repmat({''}, nC, size(processFunPath, 2));
keywords = repmat({{''}}, nC, max(1, size(processFunPath, 2)));
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
            if isfield(a, 'keywords')
                keywords{i, j} = a.keywords;
            end            
        end
    end
end

if isempty(InputBbox) && isempty(tileOutBbox) && isempty(usrFcn_strs)
    return;
end

nF = numel(inputFullpaths);

% if bigData
%     compressor = 'zstd';
% else
%     compressor = 'lz4';
% end
compressor = 'zstd';

uniq_locIds = unique(locIds);
nLoc = numel(uniq_locIds);
locSpecific = false;
if numel(uniq_locIds) > 1 && (size(InputBbox, 1) == nLoc || size(tileOutBbox, 1) == nLoc)
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
end

zarrFullpaths = cell(nF, 1);
func_strs = cell(nF, 1);
for i = 1 : nF
    inputFullpath_i = inputFullpaths{i};
    [dataPath, fsname_i] = fileparts(inputFullpath_i);
    
    inputFullpath_group_i = {inputFullpath_i};        

    zarrPath = [dataPath, '/', resultDirName];
    if ~exist(zarrPath, 'dir')
        mkdir(zarrPath);
    end

    zarrFullpaths{i} = sprintf('%s/%s/%s.zarr/', dataPath, resultDirName, fsname_i);
    if ~isempty(flippedTile)
        flipZstack = flippedTile(i);
    else
        flipZstack = false;
    end
    
    cind = cellfun(@(x) contains(inputFullpath_i, x), channelPatterns);
    if ~any(cind)
        error('The file %s does not match any channel patterns %s', inputFullpath_i, string(channelPatterns));
    end
    
    if locSpecific
        locInd = uniq_locIds == locIds(i);
        InputBbox_i = InputBbox{locInd};
        tileOutBbox_i = tileOutBbox{locInd};
        % usrFcn_str_i = usrFcn_strs{cind, locInd};
    else
        InputBbox_i = InputBbox;
        tileOutBbox_i = tileOutBbox;
        % usrFcn_str_i = usrFcn_strs{cind, 1};
    end

    % map user defined function to the tile handle based on the keywords
    usrFcn_strs_c = usrFcn_strs(cind, :);
    keywords_c = keywords(cind, :);
    usrFcn_str_i = usrFcn_strs_c{1};
    for j = 1 : numel(usrFcn_strs_c)
        break_loops = false;
        keywords_cj = keywords_c{j};
        for k = 1 : numel(keywords_cj)
            if contains(inputFullpath_i, keywords_cj{k})
                break_loops = true;
                break;
            end
        end
        if break_loops
            usrFcn_str_i = usrFcn_strs_c{j};
            break;
        end
    end

    func_strs{i} = sprintf(['stitch_process_zarr_tile(%s,''%s'',[],''BlockSize'',%s,', ...
        '''shardSize'',%s,''flipZstack'',%s,''resample'',%s,''InputBbox'',%s,', ...
        '''tileOutBbox'',%s,''compressor'',''%s'',''usrFcn'',"%s")'], ...
        sprintf('{''%s''}', strjoin(inputFullpath_group_i, ''',''')), zarrFullpaths{i}, ...
        strrep(mat2str(blockSize), ' ', ','), strrep(mat2str(shardSize), ' ', ','), ...
        string(flipZstack), strrep(mat2str(resample), ' ', ','), strrep(mat2str(InputBbox_i), ' ', ','), ...
        strrep(mat2str(tileOutBbox_i), ' ', ','), compressor, usrFcn_str_i);
end

imSizes = zeros(numel(inputFullpath_group_i), 3);
for i = 1 : numel(inputFullpath_group_i)
    imSizes(i, :) = getImageSize(inputFullpath_group_i{i});
end
imSize = [imSizes(1, 1 : 2), sum(imSizes(:, 3))];
memAllocate = prod(imSize) * 4 / 1024^3 * 2.5;
if ~bigData
    memAllocate = memAllocate * 2;
end
maxTrialNum = 2;

is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, zarrFullpaths, ...
    func_strs, 'parseCluster', parseCluster, 'masterCompute', masterCompute, ...
    'maxTrialNum', maxTrialNum,  'cpusPerTask', cpusPerTask, 'memAllocate', memAllocate, ...
    'mccMode', mccMode, 'configFile', configFile);
if ~all(is_done_flag)
    generic_computing_frameworks_wrapper(inputFullpaths, zarrFullpaths, func_strs, ...
        'parseCluster', parseCluster, 'masterCompute', masterCompute, 'maxTrialNum', maxTrialNum, ...
        'cpusPerTask', cpusPerTask * 2, 'memAllocate', memAllocate * 2, 'mccMode', mccMode, ...
        'configFile', configFile);
end


end


