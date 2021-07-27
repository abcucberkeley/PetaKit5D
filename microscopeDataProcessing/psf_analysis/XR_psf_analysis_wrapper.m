function [] = XR_psf_analysis_wrapper(dataPaths, varargin)
% psf analysis wrapper
%


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths');
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric);
ip.addParameter('angle', 32.45, @isnumeric);
ip.addParameter('Deskew', true, @islogical);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamB_ch0'}, @iscell);
ip.addParameter('Channels', [488, 560], @isnumeric);
ip.addParameter('RWFn', {'/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif', '/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_605em_128_128_101_100nmSteps.tif'}, @iscell);
ip.addParameter('sourceStr', 'test', @ischar);
% ip.addParameter('prefix', 'test_', @ischar);
ip.parse(dataPaths, varargin{:});

pr = ip.Results;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
angle = pr.angle;
Deskew = pr.Deskew;
ObjectiveScan = pr.ObjectiveScan;
ChannelPatterns = pr.ChannelPatterns;
Channels = pr.Channels;
RWFn = pr.RWFn;
sourceStr = pr.sourceStr;

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


%% deskew psfs

if Deskew
    dataPath_exps = cellfun(@(x) [x, '/'], dataPaths, 'unif', 0);
    disp(dataPath_exps);

    Save16bit = true;
    Reverse = true;

    general_options = {'xyPixelSize', xyPixelSize, ...
                       'dz' dz, ...
                       'Reverse', Reverse, ...
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
                   'flipZstack', ~true, ...
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

for d = 1 : numel(dataPath_exps)
    rtd = dataPath_exps{d};
    fn = dir([rtd '*.tif']);
    fn = {fn.name}';
    
    result_dir = [rtd, 'PSFAnalysis/'];
    mkdir(result_dir);
    
    % rt_RW = '/clusterfs/fiona/Gokul/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif';
    % rt_RW = '/Users/xruan/Images/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif';
    xypixsize= xyPixelSize * 1000;
    if ObjectiveScan
        zpixsize = dz * 1000;    
        PSFsubpix = [128, 128, round((501 - 1) * 0.04 / dz * sind(angle)) + 1];        
    else
        zpixsize = dz * sind(angle) * 1000;
        PSFsubpix = [128, 128, round((501 - 1) * 0.04 / dz) + 1];        
    end

    NAdet = 1.0;
    index = 1.33;
    gamma = 0.5;
    source_descrip = sourceStr;
    
    zpixsize_RW = 0.1 * 1000;
    PSFsubpix_RW = [128, 128, 101];

    for k = 1:numel(fn)
        ch_ind = cellfun(@(x) contains(fn{k}, x), ChannelPatterns);
        Channel_k = Channels(ch_ind);
        RWFn_k = RWFn{ch_ind};
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
        
        if exist([result_dir 'wT_' fn{k}(1:end-4) '.png'], 'file')
            continue;
        end
        
        [xz_exp_PSF_RW, xz_exp_OTF_RW, xOTF_linecut_RW, zOTF_linecut_RW, zOTF_bowtie_linecut_RW] = Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update(RWFn_k, source_descrip, xypixsize, zpixsize_RW, NAdet, index, exc_lambda, det_lambda, PSFsubpix_RW, gamma);

        [xz_exp_PSF, xz_exp_OTF, xOTF_linecut, zOTF_linecut, zOTF_bowtie_linecut] = Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update([rtd fn{k}], source_descrip, xypixsize, zpixsize, NAdet, index, exc_lambda, det_lambda, PSFsubpix, gamma);
        f0 = gcf();
        print(f0, '-painters','-dpng', '-loose',[result_dir 'comp_' fn{k}(1:end-4) '.png']);
        close all

        figure('Renderer', 'painters', 'Position', [10 10 600 600]);
        % figure;
        A = size(xz_exp_OTF);
        D = size(zOTF_linecut);
        plot(log10(zOTF_linecut), 'r', 'LineWidth', 2);hold on
        plot(log10(xOTF_linecut), 'b', 'LineWidth', 2);
        plot(log10(zOTF_bowtie_linecut), 'g', 'LineWidth', 2);
        axis([1 D(2) -3 0]);
        axis square;
        grid on;
        set(gca, 'XTick', [1:(D(2)-1)./10:D(2)]);
        set(gca, 'XTickLabel', [-1:0.2:1]);
        xlabel(['k / (4\pi/\lambda)'], 'FontSize', 14);
        set(gca, 'YTick', [-3:1:0]);
        set(gca, 'YTickLabel', 10.^[-3:1:0]);
        ylabel(['OTF Strength'], 'FontSize', 14);
        % text(-0.1 .*A(2), 0.15, ['Overall OTF linecuts From ', source_descrip], 'FontSize', 14);
        text(0.6.*A(2), -0.15, 'OTF along ky', 'Color', [0 0 1], 'FontSize', 14);
        text(0.6.*A(2), -0.3, 'OTF along kz', 'Color', [1 0 0], 'FontSize', 14);
        text(0.6.*A(2), -0.45, 'Bowtie OTF along kz', 'Color', [0 0.75 0], 'FontSize', 14);

        hold on

        D = size(zOTF_linecut_RW);
        plot(log10(zOTF_linecut_RW), '--r', 'LineWidth', 2);hold on
        plot(log10(xOTF_linecut_RW), '--b', 'LineWidth', 2);
        plot(log10(zOTF_bowtie_linecut_RW), '--g', 'LineWidth', 2);
        axis([1 D(2) -3 0]);
        axis square;
        grid on;
        set(gca, 'XTick', [1:(D(2)-1)./10:D(2)]);
        set(gca, 'XTickLabel', [-1:0.2:1]);
        xlabel(['k / (4\pi/\lambda)'], 'FontSize', 14);
        set(gca, 'YTick', [-3:1:0]);
        set(gca, 'YTickLabel', 10.^[-3:1:0]);
        ylabel(['OTF Strength'], 'FontSize', 14);
        title( ['Overall OTF linecuts From ', source_descrip], 'FontSize', 14);
        legend([{'OTF along kz','OTF along ky','Bowtie OTF along kz', 'RW OTF along kz',  'RW OTF along ky', 'RW Bowtie OTF along kz'}]);

        f0 = gcf();
        print(f0, '-painters','-dpng', '-loose', [result_dir 'wT_' fn{k}(1:end-4) '.png']);

    end
end


end

