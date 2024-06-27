function [] = XR_psf_detection_and_analysis_wrapper(dataPaths, varargin)
% Check psf and crop them from the raw image that contain many PSF patterns. 
% The idea is to detect local maximum in radon image for the skewed angle in xz to
% get the peaks of psfs. Then, remove peaks if they are close to each other
% or to the boarder. Finally, crop the psfs for kept peaks for given crop
% size. 
%
%
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%
% Parameters (as 'specifier'-value pairs):
%         xyPixelSize : a number (default: 0.108). Pixel size in um.
%                  dz : a number (default: 0.5). Scan interval in um.
%           skewAngle : a number (default: 32.45). Skew angle (in degree) of the stage.
%            cropSize : 1x3 vector (default: [256, 128, 201]). Size to crop isolated PSFs.
%          flipZstack : true|false (default: false). Flip z stacks.
%          distThresh : 1x3 vector (default: [256, 128, 201]). Distance threshold between detected PSF centers in yxz. Only keep isolated PSFs beyond the distance thresholds.
%     channelPatterns : a cell array (default: {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}).  Channel identifiers for included channels. 
%            channels : 1x#Channels (default: [488, 560]). Wavelength for the channels.
%                RWFn : a cell array with 1x#channels (default: {'', ''}). Richards Wolf theoretical Widefield PSF paths mapped to the channels.
%           sourceStr : char (default: 'test'). Source of the PSF shown in the title of the PSF analysis plots.
%             visible : true|false (default: true). Make figure visible; otherwise, plot the figure in the background.
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%         cpusPerTask : a number (default: 8). The number of cpu cores per task for slurm job submission.
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%
% 
% Author: Xiongtao Ruan (07/2021)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric);
ip.addParameter('skewAngle', 32.45, @isnumeric);
ip.addParameter('cropSize', [256, 128, 201], @isnumeric);
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('distThresh', [256, 128, 201], @isnumeric);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamB_ch0'}, @iscell);
ip.addParameter('channels', [488, 560], @isnumeric);
ip.addParameter('RWFn', {'', ''}, @iscell);
ip.addParameter('sourceStr', 'test', @ischar);
ip.addParameter('visible', true, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', false, @islogical);
ip.addParameter('cpusPerTask', 8, @isscalar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);
ip.parse(dataPaths, varargin{:});

pr = ip.Results;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
skewAngle = pr.skewAngle;
cropSize = pr.cropSize;
flipZstack = pr.flipZstack;
distThresh = pr.distThresh;
channelPatterns = pr.channelPatterns;
channels = pr.channels;
RWFn = pr.RWFn;
sourceStr = pr.sourceStr;
visible = pr.visible;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
configFile = pr.configFile;

tic
% rt = '/Users/xruan/Images/20210607_PSFs_L15_37C/';
if ischar(dataPaths)
    dataPaths = {dataPaths};
end
if ispc
    dataPaths = cellfun(@(x) strrep(x, '\', '/'), dataPaths, 'unif', 0);
end

dataPath_exps = dataPaths;
disp(dataPath_exps);

fns_cell = cell(numel(dataPath_exps), 1);
prefix_cell = cell(numel(dataPath_exps), 1);
for i = 1 : numel(dataPath_exps)
    dir_info = dir([dataPath_exps{i}, '/*.tif']);
    fns = {dir_info.name}';
    % only keep main filenames without _partXXX
    prefixes = fns(~cellfun(@isempty, regexp(fns, '^(?!.*_part[0-9]*.tif).*$')));
    prefixes = cellfun(@(x) x(1 : end - 4), prefixes, 'unif', 0);
    prefixes = unique(prefixes);
    prefix_cell{i} = prefixes;
    fns_cell{i} = cellfun(@(x) [dataPath_exps{i}, x], fns, 'unif', 0);
end
prefixes = cat(1, prefix_cell{:});
fns_cell = cat(1, fns_cell{:});

include_flag = false(numel(fns_cell), 1);
for c = 1 : numel(channelPatterns)
    include_flag = include_flag | contains(fns_cell, channelPatterns{c});
end
fns_cell = fns_cell(include_flag);

frameFullpaths = cell(numel(prefixes), 1);
outFullpaths = cell(numel(prefixes), 1);
func_strs = cell(numel(prefixes), 1);
for i = 1 : numel(prefixes)
    prefix = prefixes{i};
    fns_i = fns_cell(contains(fns_cell, prefix));
    if isempty(fns_i)
        continue;
    end
    
    dn = fileparts(fns_i{1});
    % cropSize(3) = round((cropSize_80(3) - 1) * 0.08 / dz) + 1; 
    % distThresh(3) = round((distThresh_80(3) - 1) * 0.08 / dz); 
    
    result_dir = [dn, '/Cropped/'];
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end

    frameFullpaths{i} = fns_i{1};
    outFullpaths{i} = [result_dir, prefix, '.mat'];
    % XR_psf_detection_and_cropping(fns_i, result_dir, 'xyPixelSize', xyPixelSize, 'dz', dz, ...
    %    'angle', angle, 'cropSize', cropSize, 'distThresh', distThresh, 'prefix', prefix);
   
    fns_i_str = sprintf('{''%s''}', strjoin(fns_i, ''','''));
    func_strs{i} = sprintf(['XR_psf_detection_and_cropping(%s,''%s'',''xyPixelSize'',%.20f,''dz'',%.20f,' ...
        '''skewAngle'',%.20f,''cropSize'',[%s],''distThresh'',[%s],''prefix'',''%s'')'], fns_i_str, ...
        result_dir, xyPixelSize, dz, skewAngle, strrep(num2str(cropSize, '%.10f,'), ' ', ''), ...
        strrep(num2str(distThresh, '%.10f,'), ' ', ''), prefix);
end

empty_inds = cellfun(@isempty, frameFullpaths);
frameFullpaths(empty_inds) = [];
outFullpaths(empty_inds) = [];
func_strs(empty_inds) = [];


% use cluster computing for the psf detection and cropping
memAllocate = prod(getImageSize(frameFullpaths{1})) * 4 / 1024^3 * 100;
is_done_flag = generic_computing_frameworks_wrapper(frameFullpaths, outFullpaths, ...
    func_strs, 'masterCompute', true, 'cpusPerTask', cpusPerTask, memAllocate=memAllocate, ...
    mccMode=mccMode, configFile=configFile);
if ~all(is_done_flag)
    generic_computing_frameworks_wrapper(frameFullpaths, outFullpaths, ...
        func_strs, 'masterCompute', true, 'cpusPerTask', cpusPerTask, memAllocate=memAllocate, ...
        mccMode=mccMode, configFile=configFile);
end   


% deskew psfs and psf analysis
% parameters for psf analysis
% deskew psf if it is in skewed space
Deskew = true;
% z stage scan 
zStageScan = false;
% objective scan or not
objectiveScan = false;

dataPath_exps = cellfun(@(x) [x, '/Cropped/'], dataPaths, 'unif', 0);
disp(dataPath_exps);

XR_psf_analysis_wrapper(dataPath_exps, 'dz', dz, 'skewAngle', skewAngle, 'channelPatterns', channelPatterns, ...
    'Channels', channels, 'Deskew', Deskew, 'flipZstack', flipZstack, 'objectiveScan', objectiveScan, ...
    'zStageScan', zStageScan, 'sourceStr', sourceStr, 'RWFn', RWFn, visible=visible, parseCluster=parseCluster, ...
    masterCompute=masterCompute, cpusPerTask=cpusPerTask, mccMode=mccMode, configFile=configFile);

end

