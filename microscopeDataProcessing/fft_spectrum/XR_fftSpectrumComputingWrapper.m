function XR_fftSpectrumComputingWrapper(dataPaths, varargin)
% Dataset level wrapper for fft spectrum computing.
% For each image, the default outputs are central slices and MIPs of power 
% spectrum with gamma 0.5. If save 3D stack, the log10 power spectrum is saved. 
% 
%
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%
% Parameters (as 'specifier'-value pairs):
%       resultDirName : char (default: 'FFT'). Result directory under data path.
%           overwrite : true|false (default: false). Overwrite existing results
%         xyPixelSize : a number (default: 0.108). Pixel size in um.
%                  dz : a number (default: 0.5). Scan interval in um.
%            zarrFile : true|false (default: false). Use Zarr file as input.
%        outPixelSize : empty or a number (default: []). Voxel size for output. This number applies to all three axes. If empty, the minimum of xyPixelSize and dz is used.
%             outSize : 1x3 vector (default: [1001, 1001, 1001]). Output size for the FFT spectrum computing. The image is cropped/padded to this size. Must be the same across all axes.
%     channelPatterns : a cell array (default: {}).  Channel identifiers for included channels. If not provided, process all images in dataPaths.
%         save3DStack : true|false (default: false). Save 3D stack for FFT power spectrum (in log10 scale). 
%          background : a number (default: 0). Background subtraction before FFT computing.
%        interpMethod : 'linear'|'cubic'|'lanczos2'|'lanczos3' (default: 'linear'). Interpolation method for resampling to out pixel size.
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%         cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission.
%               debug : true|false (default: false). Debug mode for simplified and OMW method. If true, save the intermediate steps every n iterations (n is defined by saveStep).
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%
%
% Author: Xiongtao Ruan (11/25/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'FFT',  @ischar);
ip.addParameter('overwrite', false,  @(x) (numel(x) == 1) && islogical(x));
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('outPixelSize', [], @(x) isnumeric(x) || isempty(x));
ip.addParameter('outSize', [1001, 1001, 1001], @isnumeric);
ip.addParameter('channelPatterns', {}, @(x) ischar(x) || iscell(x));
ip.addParameter('save3DStack', false, @islogical);
ip.addParameter('background', 0, @isnumeric);
ip.addParameter('interpMethod', 'linear', @ischar);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
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

