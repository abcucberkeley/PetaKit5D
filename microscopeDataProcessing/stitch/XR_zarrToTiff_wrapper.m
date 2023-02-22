function [] = XR_zarrToTiff_wrapper(dataPaths, varargin)
% The wrapper for convert all zarr files in dataPaths to tiff
%
% Author: Xiongtao Ruan (01/19/2021)
% 
% xruan (10/11/2020): add function handle for processing before saving to zarr

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) iscell(x) || ischar(x));
ip.addParameter('ChannelPatterns', {'CamA', 'CamB'}, @iscell);
ip.addParameter('resultDirStr', 'tiffs/', @ischar);
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2022a; matlab -nodisplay -nosplash -nodesktop -r', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);


ip.parse(dataPaths, varargin{:});

pr = ip.Results;
% Resolution = pr.Resolution;
ChannelPatterns = pr.ChannelPatterns;
resultDirStr = pr.resultDirStr;
usrFcn = pr.usrFcn;

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
    mkdir([dataPath, '/', resultDirStr]);
    
    if numel(fnames) > 0
        zarrFullpaths = cellfun(@(x) [dataPath, '/', x], fnames, 'unif', 0);
        tiffFullpaths = cellfun(@(x) sprintf('%s/%s/%s.tif', dataPath, resultDirStr, x(1 : end - 5)), fnames, 'unif', 0);    
        allZarrFullpaths{d} = zarrFullpaths;
        allTiffFullpaths{d} = tiffFullpaths;
    end
end

allZarrFullpaths = cat(1, allZarrFullpaths{:});
allTiffFullpaths = cat(1, allTiffFullpaths{:});

func_strs = arrayfun(@(x) sprintf(['zarrToTiff(''%s'',''%s'',''usrFcn'',''%s'')'], ...
    allZarrFullpaths{x}, allTiffFullpaths{x}, usrFcn), 1 : numel(allZarrFullpaths), 'unif', 0);

imSize = getImageSize(allZarrFullpaths{1});
cpusPerTask = min(24, ceil(prod(imSize) * 2 / 2^30 / 20 * 2));
is_done_flag = slurm_cluster_generic_computing_wrapper(allZarrFullpaths, allTiffFullpaths, ...
    func_strs, 'maxTrialNum', 2, 'cpusPerTask', cpusPerTask);
if ~all(is_done_flag)
    is_done_flag = slurm_cluster_generic_computing_wrapper(allZarrFullpaths, ...
        allTiffFullpaths, func_strs, 'maxTrialNum', 1, 'cpusPerTask', min(24, cpusPerTask * 2));
end


end
