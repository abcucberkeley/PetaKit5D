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
ip.addRequired('dataPaths');
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric);
ip.addParameter('angle', 32.45, @isnumeric);
ip.addParameter('Deskew', true, @islogical);
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('ZstageScan', false, @islogical);
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamB_ch0'}, @iscell);
ip.addParameter('Channels', [488, 560], @isnumeric);
ip.addParameter('Save16bit', true, @islogical);
ip.addParameter('bgFactor', 1.5, @isnumeric);
ip.addParameter('RWFn', {'/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif', '/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_605em_128_128_101_100nmSteps.tif'}, @iscell);
ip.addParameter('sourceStr', 'test', @ischar);
ip.addParameter('masterCompute', false, @islogical);
% ip.addParameter('prefix', 'test_', @ischar);
ip.parse(dataPaths, varargin{:});

pr = ip.Results;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
angle = pr.angle;
Deskew = pr.Deskew;
flipZstack = pr.flipZstack;
ObjectiveScan = pr.ObjectiveScan;
ZstageScan = pr.ZstageScan;
ChannelPatterns = pr.ChannelPatterns;
Channels = pr.Channels;
Save16bit = pr.Save16bit;
bgFactor = pr.bgFactor;
RWFn = pr.RWFn;
sourceStr = pr.sourceStr;
masterCompute = pr.masterCompute;

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

if Deskew
    % dataPath_exps = cellfun(@(x) [x, '/'], dataPaths, 'unif', 0);
    disp(dataPath_exps);

    % Save16bit = true;
    Reverse = true;

    general_options = {'xyPixelSize', xyPixelSize, ...
                       'dz' dz, ...
                       'SkewAngle', angle, ...
                       'Reverse', Reverse, ...
                       'ZstageScan', ZstageScan, ...                       
                       'ChannelPatterns', ChannelPatterns, ...
                       'Save16bit', Save16bit...
                       'Overwrite', false, ...
                       'Streaming', false, ...
                       'cpusPerTask', 8, ...
                       'cpuOnlyNodes', false, ...
                       'parseCluster', ~false
                       };

    % dsr
    % Rotate is for DSR, set Rotate as true if DSR is needed.
    dsr_options = {'Deskew', true, ...
                   'Rotate', ~true, ...
                   'DSRCombined', false, ...
                   'parseSettingFile', ~true, ...  
                   'flipZstack', flipZstack, ...
                   'LLFFCorrection', ~true,...
                  };

    % stitch
    stitch_options = {};

    % decon          
    decon_options = {'Decon', ~true};

    XR_microscopeAutomaticProcessing(dataPath_exps, general_options{:}, ...
        dsr_options{:}, stitch_options{:}, decon_options{:});

end


%% psf analysis

if Deskew
    dataPath_exps = cellfun(@(x) [x, '/DS/'], dataPaths, 'unif', 0);
end
    
disp(dataPath_exps);

NAdet = 1.0;
index = 1.33;
gamma = 0.5;
source_descrip = sourceStr;

xypixsize= xyPixelSize * 1000;
if ObjectiveScan
    zpixsize = dz * 1000;    
    % PSFsubpix = [128, 128, round((501 - 1) * 0.04 / dz * sind(angle)) + 1];   
elseif ZstageScan
    zpixsize = dz * cosd(angle) * 1000;
    % PSFsubpix = [128, 128, round((501 - 1) * 0.04 / dz) + 1];                
else
    zpixsize = dz * sind(angle) * 1000;
    % PSFsubpix = [128, 128, round((501 - 1) * 0.04 / dz) + 1];        
end
PSFsubpix = [128, 128, round(100 * 100 / zpixsize) + 1];

zpixsize_RW = 0.1 * 1000;
PSFsubpix_RW = [128, 128, 101];
bgFactor_RW = 0;

RW_info = cell(numel(ChannelPatterns), 1);

% run psf analysis for RW images
for c = 1 : numel(ChannelPatterns)
    Channel_k = Channels(c);
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
    for c = 1 : numel(ChannelPatterns)
        include_flag = include_flag | contains(fn, ChannelPatterns{c}) | contains(fn, regexpPattern(ChannelPatterns{c}));
    end
    fn = fn(include_flag);
    fsn = fsn(include_flag);
    
    frameFullpaths{d} = fn;
    figureFullpaths{d} = cellfun(@(x) [result_dir, 'wT_', x(1 : end - 4), '.png'], fsn, 'unif', 0);
    
    func_strs{d} = cell(numel(fn), 1);
    for f = 1 : numel(fn)
        ch_ind = cellfun(@(x) contains(fn{f}, x), ChannelPatterns);
        if ~any(ch_ind)
            continue;
        end
        
        Channel_k = Channels(ch_ind);
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
MatlabLaunchStr = 'module load matlab/r2022a; matlab -nodisplay -nosplash -nodesktop -r'; 
is_done_flag = slurm_cluster_generic_computing_wrapper(frameFullpaths, figureFullpaths, ...
    func_strs, 'MatlabLaunchStr', MatlabLaunchStr, 'maxTrialNum', maxTrialNum, ...
    'masterCompute', masterCompute, 'cpusPerTask', cpusPerTask);
if ~all(is_done_flag)
    slurm_cluster_generic_computing_wrapper(frameFullpaths, figureFullpaths, ...
        func_strs, 'MatlabLaunchStr', MatlabLaunchStr, 'maxTrialNum', maxTrialNum, ...
        'masterCompute', masterCompute, 'cpusPerTask', cpusPerTask * 2);
end
if ~all(is_done_flag)
    slurm_cluster_generic_computing_wrapper(frameFullpaths, figureFullpaths, ...
        func_strs, 'MatlabLaunchStr', MatlabLaunchStr, 'masterCompute', masterCompute, 'cpusPerTask', cpusPerTask * 6);
end



end

