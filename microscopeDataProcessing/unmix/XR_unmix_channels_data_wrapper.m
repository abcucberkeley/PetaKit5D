function [] = XR_unmix_channels_data_wrapper(dataPaths, varargin)
% Unmixing a channel by subtracting a weighted fraction of other channels
%
% Required inputs:
%           dataPaths : char or cell array. Path(s) to the dataset(s) to be processed. Can be a single path as a char string or multiple paths in a cell array.
%
% Parameters (as 'specifier'-value pairs):
%        unmixFactors : numeric matrix (default: []). The unmixing factors required for spectral unmixing.
%                mode : char (default: 'linear'). The unmixing mode to be used. Valid options are 'linear' or 'gaussian'.
%         unmixSigmas : numeric array (default: []). Sigmas for Gaussian unmixing when 'mode' is set to 'gaussian'.
%       resultDirName : char (default: 'Unmixed'). Name of the subdirectory within each dataPath where results will be saved.
%     channelPatterns : cell array (default: {'CamA', 'CamB'}). File extensions or patterns used to identify image files.
%          channelInd : numeric scalar (default: 1). The index of the channel to which the unmixing will be applied.
%        FFCorrection : true|false (default: false). Flat-field correction. Only for small files (largeFile = false).
%          lowerLimit : a number between 0 and 1 (default: 0.4). Lower limit to cap the flat field image.
%        FFImagePaths : empty or a cell array of paths for corresponding channels (default: {'', ''}). Flat field image paths.
%     backgroundPaths : empty or a cell array of paths for corresponding channels (default: {'', ''}). Background image paths.
%     constBackground : empty or a scalar, or a numerical array for corresponding channels (default: []). Background subtraction with constant background values.
%         constOffset : empty or a scalar (default: []). If not empy, add the const offset after unmixing.
%            zarrFile : true|false (default: false). Specifies if the input data is in Zarr format.
%           largeFile : true|false (default: false). Flag to indicate if the input files are large and require special handling.
%            saveZarr : true|false (default: false). If true, the output is saved in Zarr format; otherwise, saved in Tiff format.
%           save16bit : true|false (default: true). If true, the output is saved as 16-bit integers; otherwise, save as single.
%           batchSize : 1x3 numeric vector (default: [1024, 1024, 1024]). The size of data batches [y, x, z] for processing.
%           blockSize : 1x3 numeric vector (default: [256, 256, 256]). The chunk/block size [y, x, z] for writing data, especially for Zarr.
%          borderSize : 1x3 numeric vector (default: [0, 0, 0]). The size of the border [y, x, z] to consider between adjacent blocks/batches.
%        parseCluster : true|false (default: true). If true, sets up the parallel computing cluster environment.
%       masterCompute : true|false (default: true). If true, the master node will also participate in computation tasks.
%           jobLogDir : char (default: '../job_logs'). Directory for saving logs from cluster jobs.
%         cpusPerTask : numeric scalar (default: 3). The number of CPU cores to allocate for each parallel task.
%          configFile : char (default: ''). Path to an optional configuration file for slurm job submission settings.
%             mccMode : true|false (default: false). Set to true if running the code as a compiled application (using MATLAB Compiler).
%                uuid : char (default: ''). A unique identifier string for this processing run, often used for temporary files.
%               debug : true|false (default: false). If true, enables debug mode which may provide more verbose output.
% 
% Author: Xiongtao Ruan (09/17/2025)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('unmixFactors', [], @(x) isnumeric(x));
ip.addParameter('mode', 'linear', @ischar); % linear vs gaussian
ip.addParameter('unmixSigmas', [], @isnumeric); 
ip.addParameter('resultDirName', 'Unmixed', @ischar); 
ip.addParameter('channelPatterns', {'CamA', 'CamB'}, @iscell);
ip.addParameter('channelInd', 1, @isnumeric); % unmix for which channel
% flat field parameters
ip.addParameter('FFCorrection', false, @islogical);
ip.addParameter('lowerLimit', 0.4, @isnumeric);
ip.addParameter('FFImagePaths', {'',''}, @iscell);
ip.addParameter('backgroundPaths', {'',''}, @iscell);
% background subtraction with constant background instead of background image
ip.addParameter('constBackground', [], @isnumeric);
% constant offset after unmixing
ip.addParameter('constOffset', [], @(x) isnumeric(x));
% input and output parameters
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('saveZarr', false, @islogical);
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('batchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256] , @isvector); % in y, x, z
ip.addParameter('borderSize', [0, 0, 0] , @isvector); % in y, x, z
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 3, @isnumeric);
ip.addParameter('configFile', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
unmixFactors = pr.unmixFactors;
mode = pr.mode;
unmixSigmas = pr.unmixSigmas;
resultDirName = pr.resultDirName;
channelPatterns = pr.channelPatterns;
channelInd = pr.channelInd;
% flat field parameters
FFCorrection = pr.FFCorrection;
lowerLimit = pr.lowerLimit;
FFImagePaths = pr.FFImagePaths;
backgroundPaths = pr.backgroundPaths;
% background subtraction with constant background instead of background image
constBackground = pr.constBackground;
% constant offset after unmixing
constOffset = pr.constOffset;
% input and output parameters
zarrFile = pr.zarrFile;
largeFile = pr.largeFile;
saveZarr = pr.saveZarr;
save16bit = pr.save16bit;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
borderSize = pr.borderSize;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
configFile = pr.configFile;
mccMode = pr.mccMode;
uuid = pr.uuid;

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), filesep)
        dataPaths{d} = [dataPath, filesep];
    end

    dataPaths{d} = [simplifyPath(dataPaths{d}), '/'];
    resultPath = [dataPaths{d}, '/' resultDirName, '/'];
    resultPaths{d} = resultPath;
    mkdir(resultPath);
    fileattrib(resultPath, '+w', 'g');
end

if isempty(uuid)
    uuid = get_uuid();
end

% parse image filenames
[fnames, fsns, fd_inds, filepaths, ch_inds] = parseImageFilenames(dataPaths, zarrFile, channelPatterns);
nc = numel(channelPatterns);

uniq_ch = unique(ch_inds);
if numel(uniq_ch) ~= nc
    error('Number of channels of the images does not match the number of channel patterns!');
end
if nc ~= numel(unmixFactors)
    error('The number of images does not match that of the unmixing factors!');
end
if FFCorrection && (nc ~= numel(FFImagePaths) || nc ~= numel(backgroundPaths))
    error('The number of flat field or background images does not match number of channels!');
end
if ~isempty(constBackground)
    if isscalar(constBackground)
        constBackground = ones(1, nc) * constBackground;
    elseif(nc ~= numel(constBackground))
        error('The number of constBackground values does not match number of channels!');
    end
end

% group filenames by channels. If each folder contain a channel, group all
% files, otherwise, group by each data folder.
single_channel_per_folder = nd > 1 && (nc == nd) && all(arrayfun(@(x) isscalar(unique(ch_inds(fd_inds == x))), 1 : nd));
if single_channel_per_folder
    fd_ch_groups = groupImageFilenamesByChannels(fnames, channelPatterns, ch_inds);
else
    fd_ch_groups_cell = cell(nd, 1);
    for d = 1 : nd
        fd_ind_d = find(fd_inds == d);
        fnames_d = fnames(fd_ind_d);
        fd_ch_group_d = groupImageFilenamesByChannels(fnames_d, channelPatterns, ch_inds);
        if size(fd_ch_group_d, 1) == 1
            fd_ch_groups_cell{d} = fd_ind_d(fd_ch_group_d)';
        else
            fd_ch_groups_cell{d} = fd_ind_d(fd_ch_group_d);
        end
    end
    fd_ch_groups = cat(1, fd_ch_groups_cell{:});
end

% exclude groups that don't have all channels
incomplete_groups = find(any(fd_ch_groups == 0, 2), 1);
if ~isempty(incomplete_groups)
    warning('Some images do not have all channels and are excluded for unmixing!');
end

primary_fd_inds = fd_ch_groups(:, channelInd);

ext = '.tif';
if saveZarr
    ext = '.zarr';
end

imSize = getImageSize([dataPaths{fd_inds(1)}, '/', fnames{1}]);

ng = numel(primary_fd_inds);
frameFullpaths = cell(ng, 1);
resultFullpaths = cell(ng, 1);
func_strs = cell(ng, 1);
for f = 1 : ng
    frameFullpaths_f = arrayfun(@(x) sprintf('%s/%s', dataPaths{fd_inds(x)}, fnames{x}), ...
        fd_ch_groups(f, :), 'UniformOutput', false);
    frameFullpaths_str = sprintf('{''%s''}', strjoin(frameFullpaths_f, ''','''));
    frameFullpaths{f} = frameFullpaths_f{channelInd};
    p_ind = fd_ch_groups(f, channelInd);
    resultFullpaths{f} = sprintf('%s/%s%s', resultPaths{fd_inds(p_ind)}, fsns{p_ind}, ext);

    if largeFile    
        func_strs{f} = sprintf(['XR_unmix_channels_zarr(%s,%s,''mode'',''%s'',', ...
            '''unmixSigmas'',%s,''resultDirName'',''%s'',''channelInd'',%d,''save16bit'',%s,', ...
            '''batchSize'',%s,''blockSize'',%s,''borderSize'',%s,''parseCluster'',%s,', ...
            '''masterCompute'',%s,''jobLogDir'',''%s'',''cpusPerTask'',%d,''configFile'',''%s'',', ...
            '''mccMode'',%s,''uuid'',''%s'')'], frameFullpaths_str, mat2str_comma(unmixFactors, 10), ...
            mode, mat2str_comma(unmixSigmas, 10), resultDirName, channelInd, ...
            string(save16bit), mat2str_comma(batchSize), mat2str_comma(blockSize), ...
            mat2str_comma(borderSize), string(parseCluster), string(masterCompute), ...
            jobLogDir, cpusPerTask, configFile, string(mccMode), uuid);
    else
        FFImagePaths_str = sprintf('{''%s''}', strjoin(FFImagePaths, ''','''));
        backgroundPaths_str = sprintf('{''%s''}', strjoin(backgroundPaths, ''','''));

        func_strs{f} = sprintf(['XR_unmix_channels_frame(%s,%s,''mode'',''%s'',', ...
            '''unmixSigmas'',%s,''resultDirName'',''%s'',''channelInd'',%d,', ...
            '''FFCorrection'',%s,''lowerLimit'',%.20d,''FFImagePaths'',%s,', ...
            '''backgroundPaths'',%s,''constBackground'',%s,''constOffset'',%s,', ...
            '''save16bit'',%s,''zarrFile'',%s,''saveZarr'',%s,''blockSize'',%s,', ...
            '''uuid'',''%s'')'], frameFullpaths_str, mat2str_comma(unmixFactors, 10), ...
            mode, mat2str_comma(unmixSigmas, 10), resultDirName, channelInd, ...
            string(FFCorrection), lowerLimit, FFImagePaths_str, backgroundPaths_str, ...
            mat2str_comma(constBackground), mat2str_comma(constOffset), string(save16bit), ...
            string(zarrFile), string(saveZarr), mat2str_comma(blockSize), uuid);
    end
end

if largeFile
    memAllocate = prod(batchSize) / 1024^3 * 4 * 2.5;
else
    memAllocate = prod(imSize) / 1024^3 * 4 * 2.5;
end

generic_computing_frameworks_wrapper(frameFullpaths, resultFullpaths, func_strs, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, memAllocate=memAllocate, mccMode=mccMode, configFile=configFile);

end
