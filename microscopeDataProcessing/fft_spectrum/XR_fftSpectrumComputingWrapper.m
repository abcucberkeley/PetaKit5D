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
ip.addRequired('dataPaths'); % data structure from loadConditionData
ip.addParameter('Overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 5) && islogical(x));
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric); % actual pixel size in z
ip.addParameter('outPixelSize', [], @(x) isnumeric(x) || isempty(x)); % output pixel size
ip.addParameter('N', [1001, 1001, 1001], @isnumeric);
ip.addParameter('ChannelPatterns', {}, @(x) ischar(x) || iscell(x));
ip.addParameter('save3DStack', false, @islogical);
ip.addParameter('background', 0, @isnumeric);
ip.addParameter('Interp', 'linear', @ischar);

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root of XR_Repository. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../../'];
cd(repo_rt);

pr = ip.Results;
Overwrite = pr.Overwrite;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
outPixelSize = pr.outPixelSize;
N = pr.N;
ChannelPatterns = pr.ChannelPatterns;
save3DStack = pr.save3DStack;
background = pr.background;
Interp = pr.Interp;

if ischar(ChannelPatterns)
    ChannelPatterns = {ChannelPatterns};
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
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), filesep)
        dataPaths{d} = [dataPath, filesep];
    end
end

% retrieve files and make directory for each data folder. 
allFullpaths = cell(numel(dataPaths), 1);
allSpctmFullpaths = cell(numel(dataPaths), 1);
for i = 1 : numel(dataPaths)
    dataPath = dataPaths{i};
    
    dPath = dataPath;
    
    FFTPath = [dPath, 'FFT/'];
    mkdir(FFTPath);
    fileattrib(FFTPath, '+w', 'g');
    save('-v7.3', [FFTPath, '/parameters.mat'], 'pr');    

    dir_info = dir([dPath, '*tif']);
    fnames = {dir_info.name}';
    dFullpaths = cellfun(@(x) [dPath, x], fnames, 'unif', 0);
    spctmFullnames = cellfun(@(x) [FFTPath, x], fnames, 'unif', 0);

    if ~isempty(ChannelPatterns)
        include_flag = false(numel(fnames), 1);
        for c = 1 : numel(ChannelPatterns)
            include_flag = include_flag | contains(dFullpaths, ChannelPatterns{c}) | ...
                contains(dFullpaths, regexpPattern(ChannelPatterns{c}));
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
    '''xyPixelSize'',%0.20d,''dz'',%0.20d,''outPixelSize'',%0.20d,''N'',%s,', ...
    '''save3DStack'',%s,''background'',%d,''Interp'',''%s'')'], ...
    allFullpaths{x}, allSpctmFullpaths{x}, xyPixelSize, dz, outPixelSize, ...
    strrep(mat2str(N), ' ', ','), string(save3DStack), background, Interp), ...
    1 : numel(allFullpaths), 'unif', 0);

slurm_cluster_generic_computing_wrapper(allFullpaths, allSpctmFullpaths, func_strs, 'cpusPerTask', 6);


end

