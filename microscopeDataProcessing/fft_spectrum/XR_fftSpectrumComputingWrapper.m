function XR_fftSpectrumComputingWrapper(dataPaths, varargin)
% data wrapper for fft spectrum computing
% 
% Author: Xiongtao Ruan (11/25/2020)
% xruan (11/11/2021): add support for rescale to isotropic and crop/pad to
% given size (1001 in each dimension)
% xruan (12/13/2021): add background subtraction
% xruan (01/03/2022): add support for output voxel size
% xruan (07/18/2022): add support for different interpolation methods for resampling


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x)); % data structure from loadConditionData
ip.addParameter('resultDirName', 'FFT',  @ischar);
ip.addParameter('overwrite', false,  @(x) (numel(x) == 1) && islogical(x));
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric); % actual pixel size in z
ip.addParameter('zarrFile', false, @islogical); % input as zarr
ip.addParameter('outPixelSize', [], @(x) isnumeric(x) || isempty(x)); % output pixel size
ip.addParameter('outSize', [1001, 1001, 1001], @isnumeric);
ip.addParameter('channelPatterns', {}, @(x) ischar(x) || iscell(x));
ip.addParameter('save3DStack', false, @islogical);
ip.addParameter('background', 0, @isnumeric);
ip.addParameter('interpMethod', 'linear', @ischar);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 3, @isscalar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root of XR_Repository. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../../'];
cd(repo_rt);

pr = ip.Results;
resultDirName = pr.resultDirName;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
zarrFile = pr.zarrFile;
outPixelSize = pr.outPixelSize;
outSize = pr.outSize;
channelPatterns = pr.channelPatterns;
save3DStack = pr.save3DStack;
background = pr.background;
interpMethod = pr.interpMethod;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(channelPatterns)
    channelPatterns = {channelPatterns};
end

% suppress directory exists warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');

if isempty(outPixelSize)
    outPixelSize = min(xyPixelSize, dz);
    fprintf('Out pixel size: %d um\n', outPixelSize);
end

if ischar(dataPaths)
    dataPaths = {dataPaths};
end
nd = numel(dataPaths);
if ispc
    for d = 1 : nd
        dataPaths{d} = strrep(dataPaths{d}, '\', '/');
    end
end
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), '/')
        dataPaths{d} = [dataPath, '/'];
    end
end

% retrieve files and make directory for each data folder. 
allFullpaths = cell(numel(dataPaths), 1);
allSpctmFullpaths = cell(numel(dataPaths), 1);
for i = 1 : numel(dataPaths)
    dataPath = dataPaths{i};
    
    dPath = dataPath;
    
    FFTPath = [dPath, resultDirName, '/'];
    mkdir(FFTPath);
    try
        fileattrib(FFTPath, '+w', 'g');
    catch ME
        disp(ME);
    end
    save('-v7.3', [FFTPath, '/parameters.mat'], 'pr');    
    
    if zarrFile
        dir_info = dir([dPath, '*zarr']);
    else
        dir_info = dir([dPath, '*tif']);
    end

    fnames = {dir_info.name}';
    dFullpaths = cellfun(@(x) [dPath, x], fnames, 'unif', 0);
    spctmFullnames = cellfun(@(x) [FFTPath, x], fnames, 'unif', 0);

    if ~isempty(channelPatterns)
        include_flag = false(numel(fnames), 1);
        for c = 1 : numel(channelPatterns)
            include_flag = include_flag | contains(dFullpaths, channelPatterns{c}) | ...
                contains(dFullpaths, regexpPattern(channelPatterns{c}));
        end
        dFullpaths = dFullpaths(include_flag);
        spctmFullnames = spctmFullnames(include_flag);
    end

    allFullpaths{i, 1} = dFullpaths;
    allSpctmFullpaths{i, 1} = spctmFullnames;   
end

allFullpaths = cat(1, allFullpaths{:});
allSpctmFullpaths = cat(1, allSpctmFullpaths{:});

func_strs = arrayfun(@(x) sprintf(['XR_fftSpectrumComputingFrame(''%s'',''%s'',', ...
    '''xyPixelSize'',%0.20d,''dz'',%0.20d,''zarrFile'',%s,''outPixelSize'',%0.20d,', ...
    '''outSize'',%s,''save3DStack'',%s,''background'',%d,''interpMethod'',''%s'')'], ...
    allFullpaths{x}, allSpctmFullpaths{x}, xyPixelSize, dz, string(zarrFile), ...
    outPixelSize, strrep(mat2str(outSize), ' ', ','), string(save3DStack), ...
    background, interpMethod), 1 : numel(allFullpaths), 'unif', 0);

sz = getImageSize(allFullpaths{1});
memAllocate = prod(max(sz, outSize)) * 4 / 1024^3 * 15;
[is_done_flag] = generic_computing_frameworks_wrapper(allFullpaths, allSpctmFullpaths, ...
    func_strs, memAllocate=memAllocate, parseCluster=parseCluster, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, mccMode=mccMode, configFile=configFile);
if ~all(is_done_flag)
    [is_done_flag] = generic_computing_frameworks_wrapper(allFullpaths, allSpctmFullpaths, ...
        func_strs, memAllocate=memAllocate * 2, parseCluster=parseCluster, masterCompute=masterCompute, ...
        cpusPerTask=cpusPerTask, mccMode=mccMode, configFile=configFile);
end

end

