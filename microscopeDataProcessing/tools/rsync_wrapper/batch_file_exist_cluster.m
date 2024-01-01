function [file_exist_mat] = batch_file_exist_cluster(fileFullpaths, tmpDir, varargin)
% check large number of files exist with multiple number of jobs
% 
% Author: Xiongtao Ruan (03/19/2022)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fileFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('tmpDir', @(x) ischar(x));
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(fileFullpaths, tmpDir, varargin{:});

pr = ip.Results;

batchSize = 50000;
nF = numel(fileFullpaths);
numBatch = ceil(nF / batchSize);

% uuid = get_uuid();
flag_dir = sprintf('%s/file_exist_flags/', tmpDir);
mkdir_recursive(flag_dir);

paramFullpath = [flag_dir, '/parameters.mat'];
save('-v7.3', paramFullpath, 'pr');

flagFullpaths = cell(numBatch, 1);
func_strs = cell(numBatch, 1);
for i = 1 : numBatch
    s = (i - 1) * batchSize + 1;
    t = min(i * batchSize, nF);
    
    inputFullpath_str = sprintf('{''%s''}', strjoin(fileFullpaths(s : t), ''','''));
    
    flagFullpath = sprintf('%s/file_exist_%d_%d.mat', flag_dir, s, t);
    flagFullpaths{i} = flagFullpath;

    func_strs{i} = sprintf(['batch_file_exist(%s,''%s'',true)'], ...
        inputFullpath_str, flagFullpath);
end

inputFlagpaths = repmat({paramFullpath}, numBatch, 1);
cpusPerTask = 1;
memAllocate = 20;
maxJobNum = inf;
is_done_flag = generic_computing_frameworks_wrapper(inputFlagpaths, flagFullpaths, ...
    func_strs, 'tmpDir', flag_dir, 'language', 'matlab', 'MatlabLaunchStr', MatlabLaunchStr, ...
    'cpusPerTask', cpusPerTask, 'memAllocate', memAllocate, 'maxJobNum', maxJobNum, ...
    'masterCompute', ~false, 'mccMode', mccMode, 'ConfigFile', ConfigFile);

% collect results
file_exist_mat = cell(numBatch, 1);
if all(is_done_flag)
    for i = 1 : numBatch
        a = load(flagFullpaths{i}, 'file_exist_mat');
        file_exist_mat{i} = a.file_exist_mat;
    end
    file_exist_mat = cat(1, file_exist_mat{:});
end

end
