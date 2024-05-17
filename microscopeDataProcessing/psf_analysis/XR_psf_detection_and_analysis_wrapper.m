function [] = XR_psf_detection_and_analysis_wrapper(dataPaths, varargin)
% check psf and crop them from the raw calibration image. 
%
% The idea is to detect local maximum in radon image for the skewed angle in xz to
% get the peaks of psfs. Then, remove peaks if they are close to each other
% or to the boarder. Finally, crop the psfs for kept peaks for given crop
% size. 
% 
% xruan (07/29/2021): change deskew and analysis step with
% XR_psf_analysis_wrapper.m

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
ip.addParameter('RWFn', {'/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif', '/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_605em_128_128_101_100nmSteps.tif'}, @iscell);
ip.addParameter('sourceStr', 'test', @ischar);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', false, @islogical);
% ip.addParameter('prefix', 'test_', @ischar);
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
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
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
cpusPerTask = 24;
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
    'zStageScan', zStageScan, 'sourceStr', sourceStr, 'RWFn', RWFn, parseCluster=parseCluster, ...
    masterCompute=masterCompute, mccMode=mccMode, configFile=configFile);

end

