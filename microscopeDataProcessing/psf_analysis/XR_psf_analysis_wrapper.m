function [] = XR_psf_analysis_wrapper(dataPaths, varargin)
% psf analysis wrapper
%
% xruan (07/27/2021): add support for z-stage scan
% xruan (07/28/2021): save RW line cut info to avoid the computing in each iteration, 
% and add parallel computing for plotting
% xruan (08/16/2021): add support for flipped psfs
% xruan (12/21/2021): add support for background subtraction factor
% xruan (04/12/2022): update z size in PSFsubpix to match it for RW image


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric);
ip.addParameter('skewAngle', 32.45, @isnumeric);
ip.addParameter('deskew', true, @islogical);
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('zStageScan', false, @islogical);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamB_ch0'}, @iscell);
ip.addParameter('channels', [488, 560], @isnumeric);
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('bgFactor', 1.5, @isnumeric);
ip.addParameter('RWFn', {'/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif', '/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_605em_128_128_101_100nmSteps.tif'}, @iscell);
ip.addParameter('sourceStr', 'test', @ischar);
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
deskew = pr.deskew;
flipZstack = pr.flipZstack;
objectiveScan = pr.objectiveScan;
zStageScan = pr.zStageScan;
channelPatterns = pr.channelPatterns;
channels = pr.channels;
save16bit = pr.save16bit;
bgFactor = pr.bgFactor;
RWFn = pr.RWFn;
sourceStr = pr.sourceStr;
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

dataPath_exps = cellfun(@(x) [x, '/'], dataPaths, 'unif', 0);
disp(dataPath_exps);


%% deskew psfs

if deskew
    % dataPath_exps = cellfun(@(x) [x, '/'], dataPaths, 'unif', 0);
    disp(dataPath_exps);

    % save16bit = true;
    Reverse = true;

    general_options = {'xyPixelSize', xyPixelSize, ...
                       'dz' dz, ...
                       'SkewAngle', skewAngle, ...
                       'Reverse', Reverse, ...
                       'zStageScan', zStageScan, ...                       
                       'channelPatterns', channelPatterns, ...
                       'save16bit', save16bit...
                       'Overwrite', false, ...
                       'Streaming', false, ...
                       'cpusPerTask', cpusPerTask, ...
                       'parseCluster', parseCluster, ...
                       'masterCompute', masterCompute, ...
                       'mccMode', mccMode, ...
                       'configFile', configFile, ...
                       };

    % dsr
    % Rotate is for DSR, set Rotate as true if DSR is needed.
    dsr_options = {'Deskew', true, ...
                   'Rotate', false, ...
                   'DSRCombined', false, ...
                   'parseSettingFile', false, ...  
                   'flipZstack', flipZstack, ...
                   'FFCorrection', false,...
                  };

    % stitch
    stitch_options = {};

    XR_microscopeAutomaticProcessing(dataPath_exps, general_options{:}, ...
        dsr_options{:}, stitch_options{:});

end


%% psf analysis

if deskew
    dataPath_exps = cellfun(@(x) [x, '/DS/'], dataPaths, 'unif', 0);
end
    
disp(dataPath_exps);

NAdet = 1.0;
index = 1.33;
gamma = 0.5;
source_descrip = sourceStr;

xypixsize= xyPixelSize * 1000;
if objectiveScan
    zpixsize = dz * 1000;    
    % PSFsubpix = [128, 128, round((501 - 1) * 0.04 / dz * sind(angle)) + 1];   
elseif zStageScan
    zpixsize = dz * cosd(skewAngle) * 1000;
    % PSFsubpix = [128, 128, round((501 - 1) * 0.04 / dz) + 1];                
else
    zpixsize = dz * sind(skewAngle) * 1000;
    % PSFsubpix = [128, 128, round((501 - 1) * 0.04 / dz) + 1];        
end
PSFsubpix = [128, 128, round(100 * 100 / zpixsize) + 1];

zpixsize_RW = 0.1 * 1000;
PSFsubpix_RW = [128, 128, 101];
bgFactor_RW = 0;

RW_info = cell(numel(channelPatterns), 1);

% run psf analysis for RW images
for c = 1 : numel(channelPatterns)
    Channel_k = channels(c);
    RWFn_k = RWFn{c};
    switch Channel_k
        case 488
            exc_lambda = 488;
            det_lambda = 515;
        case 560
            exc_lambda = 560;
            det_lambda = 605;
        case 642
            exc_lambda = 642;
            det_lambda = 680;                
    end

    [xy_exp_PSF_RW, xz_exp_PSF_RW, yz_exp_PSF_RW, xy_exp_OTF_RW, xz_exp_OTF_RW, ...
        yz_exp_OTF_RW, xOTF_linecut_RW, yOTF_linecut_RW, zOTF_linecut_RW, zOTF_bowtie_linecut_RW, ...
        zOTF_bowtie_linecut_yz] = Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update(RWFn_k, source_descrip, ...
        xypixsize, zpixsize_RW, NAdet, index, exc_lambda, det_lambda, PSFsubpix_RW, gamma, bgFactor_RW);
    
    RW_info{c} = {xz_exp_PSF_RW, xz_exp_OTF_RW, xOTF_linecut_RW, yOTF_linecut_RW, zOTF_linecut_RW, zOTF_bowtie_linecut_RW};  
end

close all;

nd = numel(dataPath_exps);

% save RW info to each psf analysis folder
RW_info_fullnames = cell(nd, 1);
for d = 1 : nd
    rtd = dataPath_exps{d};    
    result_dir = [rtd, 'PSFAnalysis', filesep];
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end
    
    RW_info_fullnames{d} = sprintf('%s/RW_info.mat', result_dir);
    save('-v7.3', RW_info_fullnames{d}, 'RW_info');
end

frameFullpaths = cell(nd, 1);
figureFullpaths = cell(nd, 1);
func_strs = cell(nd, 1);

for d = 1 : numel(dataPath_exps)
    rtd = dataPath_exps{d};
    result_dir = [rtd, 'PSFAnalysis', filesep];
    
    dir_info = dir([rtd '*.tif']);
    fsn = {dir_info.name}';
    fn = cellfun(@(x) [rtd, x], fsn, 'unif', 0);
    
    include_flag = false(numel(fn), 1);
    for c = 1 : numel(channelPatterns)
        include_flag = include_flag | contains(fn, channelPatterns{c}) | contains(fn, regexpPattern(channelPatterns{c}));
    end
    fn = fn(include_flag);
    fsn = fsn(include_flag);
    
    frameFullpaths{d} = fn;
    figureFullpaths{d} = cellfun(@(x) [result_dir, 'wT_', x(1 : end - 4), '.png'], fsn, 'unif', 0);
    
    func_strs{d} = cell(numel(fn), 1);
    for f = 1 : numel(fn)
        ch_ind = cellfun(@(x) contains(fn{f}, x), channelPatterns);
        if ~any(ch_ind)
            continue;
        end
        
        Channel_k = channels(ch_ind);
        switch Channel_k
            case 488
                exc_lambda = 488;
                det_lambda = 515;
            case 560
                exc_lambda = 560;
                det_lambda = 605;
            case 642
                exc_lambda = 642;
                det_lambda = 680;                
        end
        
        func_strs{d}{f} = sprintf(['XR_psf_analysis_plot(''%s'',''%s'',''%s'',', ...
            '%d,''%s'',%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%s,%.20f,%.20f)'], fn{f}, ...
            figureFullpaths{d}{f}, RW_info_fullnames{d}, find(ch_ind), source_descrip, xypixsize, ...
            zpixsize, NAdet, index, exc_lambda, det_lambda, strrep(mat2str(PSFsubpix), ' ', ','), ...
            gamma, bgFactor);
    end
end

frameFullpaths = cat(1, frameFullpaths{:});
figureFullpaths = cat(1, figureFullpaths{:});
func_strs = cat(1, func_strs{:});

% use cluster computing for the psf analysis
cpusPerTask = 4;
maxTrialNum = 2;
is_done_flag = false;
memAllocate = prod(getImageSize(frameFullpaths{1})) * 4 / 1024^3 * 100;

for i = 1 : 3
    if all(is_done_flag)
        break;
    end

    is_done_flag = generic_computing_frameworks_wrapper(frameFullpaths, figureFullpaths, ...
        func_strs, maxTrialNum=maxTrialNum, parseCluster=parseCluster, masterCompute=masterCompute, ...
        cpusPerTask=cpusPerTask, memAllocate=memAllocate * i, mccMode=mccMode, ...
        configFile=configFile);
end


end

