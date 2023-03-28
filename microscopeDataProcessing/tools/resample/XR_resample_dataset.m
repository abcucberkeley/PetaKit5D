function [] = XR_resample_dataset(dataPath, resultPath, rsfactor, varargin)
% resample a dataset via cluster computing
% 
%
% Author: Xiongtao Ruan (01/15/2021)
% 
% xruan (03/09/2022): add support for zarr


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @ischar);
ip.addRequired('resultPath', @ischar);
ip.addRequired('rsfactor', @isnumeric);
ip.addParameter('Interp', 'linear', @ischar);
ip.addParameter('Save16bit', true, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('saveZarr', false, @islogical); % use zarr file as output
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPath, resultPath, rsfactor, varargin{:});

warning('off', 'MATLAB:MKDIR:DirectoryExists');

pr = ip.Results;
Interp = pr.Interp;
Save16bit = pr.Save16bit;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

mkdir(resultPath);
fileattrib(resultPath, '+w', 'g');
save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
% temporary directory for intermediate results
dataPath = simplifyPath(dataPath);
resultPath = simplifyPath(resultPath);

% processing file paths and bbox for moving option. 
ext = '.tif';
if zarrFile
    ext = '.zarr';
end
dir_info = dir([dataPath, '/', '*', ext]);
fnames = {dir_info.name}';
nF = numel(fnames);

%% use generic framework for the cropping do computing

frameFullpaths = cellfun(@(x) [dataPath, '/', x], fnames, 'unif', 0);

ext = '.tif';
if saveZarr
    ext = '.zarr';
end
[~, fsns] = fileparts(fnames);
if nF == 1
    fsns = {fsns};
end

resultFullpaths = cellfun(@(x) [resultPath, '/', x, ext], fsns, 'unif', 0);
    
func_strs = arrayfun(@(x) sprintf(['XR_resampleFrame(''%s'',''%s'',[%s],''zarrFile'',%s,', ...
    '''saveZarr'',%s,''Interp'',''%s'',''Save16bit'',%s)'], frameFullpaths{x}, ...
    resultFullpaths{x}, strrep(num2str(rsfactor, '%.10f,'), ' ', ''), ...
    string(zarrFile), string(saveZarr), Interp, string(Save16bit)), 1 : nF, 'unif', false);

sz = getImageSize(frameFullpaths{1});
memAllocate = prod(sz) * 4 / 1024^3 * (2 + 2 / prod(rsfactor));
generic_computing_frameworks_wrapper(frameFullpaths, resultFullpaths, func_strs, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, memAllocate=memAllocate, mccMode=mccMode, ConfigFile=ConfigFile);

end
