function [] = XR_zarrToTiff_wrapper(dataPaths, varargin)
% The wrapper for convert all zarr files in dataPaths to tiff
%
% Author: Xiongtao Ruan (01/19/2021)
% 
% xruan (10/11/2020): add function handle for processing before saving to zarr

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) iscell(x) || ischar(x));
ip.addParameter('resultDirName', 'tiffs', @ischar);
ip.addParameter('channelPatterns', {'CamA', 'CamB'}, @iscell);
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
channelPatterns = pr.channelPatterns;
usrFcn = pr.usrFcn;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
allZarrFullpaths = cell(nd, 1);
allTiffFullpaths = cell(nd, 1);

% get all zarr files and the proposed tiff files
for d = 1 : nd
    dataPath = dataPaths{d};
    dir_info = dir([dataPath, '/', '*.zarr']);
    fnames = {dir_info.name}';
    mkdir([dataPath, '/', resultDirName]);
    
    if numel(fnames) > 0
        zarrFullpaths = cellfun(@(x) [dataPath, '/', x], fnames, 'unif', 0);
        tiffFullpaths = cellfun(@(x) sprintf('%s/%s/%s.tif', dataPath, resultDirName, x(1 : end - 5)), fnames, 'unif', 0);
        allZarrFullpaths{d} = zarrFullpaths;
        allTiffFullpaths{d} = tiffFullpaths;
    end
end

allZarrFullpaths = cat(1, allZarrFullpaths{:});
allTiffFullpaths = cat(1, allTiffFullpaths{:});

% filter filenames by channel patterns
include_flag = false(numel(allZarrFullpaths), 1);
for c = 1 : numel(channelPatterns)
    include_flag = include_flag | contains(allZarrFullpaths, channelPatterns{c}) | contains(allZarrFullpaths, regexpPattern(channelPatterns{c}));
end
allZarrFullpaths = allZarrFullpaths(include_flag);
allTiffFullpaths = allTiffFullpaths(include_flag);
nF = numel(allZarrFullpaths);

func_strs = arrayfun(@(x) sprintf(['zarrToTiff(''%s'',''%s'',''usrFcn'',''%s'')'], ...
    allZarrFullpaths{x}, allTiffFullpaths{x}, usrFcn), 1 : nF, 'unif', 0);

imSizes = zeros(numel(allZarrFullpaths), 3);
for i = 1 : numel(allZarrFullpaths)
    imSizes(i, :) = getImageSize(allZarrFullpaths{i});
end
imSize = [imSizes(1, 1 : 2), sum(imSizes(:, 3))];
memAllocate = prod(imSize) * 4 / 1024^3 * 2.5;

is_done_flag = generic_computing_frameworks_wrapper(allZarrFullpaths, allTiffFullpaths, ...
    func_strs, 'maxTrialNum', 2, 'memAllocate', memAllocate, 'mccMode', mccMode, 'configFile', configFile);
if ~all(is_done_flag)
    is_done_flag = generic_computing_frameworks_wrapper(allZarrFullpaths, ...
        allTiffFullpaths, func_strs, 'maxTrialNum', 1, 'memAllocate', memAllocate, ...
        'mccMode', mccMode, 'configFile', configFile);
end


end
