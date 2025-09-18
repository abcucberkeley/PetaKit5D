function XR_FSC_analysis_wrapper(dataPaths, varargin)
% Wrapper for datasets level FSC analysis, visualization of mean resolution (0 for x (or y), pi/2 for z).
%
%
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%
% Parameters (as 'specifier'-value pairs):
%       resultDirName : char (default: 'FSCs'). Result directory under data paths.
%         xyPixelSize : a number (default: 0.108). Pixel size in um.
%                  dz : a number (default: 0.5). Scan interval in um.
%                  dr : a number (default: 10). Radial interval in pixel.
%              dtheta : a number between 0 and 2 pi (default: pi/12). Angular interval. 
%     resThreshMethod : 'fixed'|'half-bit'|'one-bit' (default: 'fixed'). Resolution thresholding method. 'fixed': fixed threshold; 'half-bit' and 'one-bit': half or one bit thresholding.
%           resThresh : a number (default: 0.2). Resolution threshold for 'fixed' thresholding method.
%            halfSize : 1x3 vector (default: [251, 251, 251]). Half image size for FSC computing. Must be the same number for all three axes, and no greater than half of the image size.
%           inputBbox : empty or 1x6 vector (default: []). Input bounding box for crop. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%             resAxis : 'xz'|'yz' (default: 'xz'). FSC computing major axes.
%      skipConeRegion : true|false (default: true). Skip the cone region along z axis (containing noisy spectrum).
%     channelPatterns : a cell array (default: {'tif'}).  Channel identifiers for included channels. 
%            channels : 1x#Channels (default: [488, 560]). Wavelength for the channels.
%        multiRegions : true|false (default: false). Select multiple regions for FSC resolution computing.
%      regionInterval : 1x3 vector (default: [50, 50, -1]). Region interval for multi-region FSC analysis in yxz order. -1 means only the center.
%          regionGrid : empty or nx3 vector (default: []). User provided grid for region centers for multi-region FSC analysis.
%             clipPer : empty or a number. Clip the high intensity voxels based on the given percentile. If empty, not clip.
%              suffix : char (default: 'decon'). Suffix for the figure titles.
%        iterInterval : a number (default: 5). Iteration interval for FSC resolution plot.
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%         cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission.
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%
%
% Author: Xiongtao Ruan (12/10/2021)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'FSCs', @ischar);
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric);
ip.addParameter('dr', 10, @isnumeric);
ip.addParameter('dtheta', pi / 12 , @isnumeric);
ip.addParameter('resThreshMethod', 'fixed', @ischar);
ip.addParameter('resThresh', 0.2, @isnumeric);
ip.addParameter('halfSize', [251, 251, 251], @isnumeric);
ip.addParameter('inputBbox', [], @isnumeric);
ip.addParameter('resAxis', 'xz', @ischar);
ip.addParameter('skipConeRegion', true, @islogical);
ip.addParameter('channelPatterns', {'tif'}, @iscell);
ip.addParameter('channels', [488, 560], @isnumeric);
ip.addParameter('multiRegions', false, @islogical);
ip.addParameter('regionInterval', [50, 50, -1], @isnumeric);
ip.addParameter('regionGrid', [], @isnumeric);
ip.addParameter('clipPer', [], @isnumeric);
ip.addParameter('suffix', 'decon', @ischar);
ip.addParameter('iterInterval', 5, @isnumeric);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('cpusPerTask', 4, @isscalar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
dr = pr.dr;
dtheta = pr.dtheta;
halfSize = pr.halfSize;
inputBbox = pr.inputBbox;
resThreshMethod = pr.resThreshMethod;
resThresh = pr.resThresh;
resAxis = pr.resAxis;
skipConeRegion = pr.skipConeRegion;
channelPatterns = pr.channelPatterns;
channels = pr.channels;
multiRegions = pr.multiRegions;
regionInterval = pr.regionInterval;
regionGrid = pr.regionGrid;
clipPer = pr.clipPer;
suffix = pr.suffix;
iterInterval = pr.iterInterval;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
configFile = pr.configFile;

tic
if ischar(dataPaths)
    dataPaths = {dataPaths};
end
if ispc
    dataPaths = cellfun(@(x) strrep(x, '\', '/'), dataPaths, 'unif', 0);
end

dataPath_exps = cellfun(@(x) [x, '/'], dataPaths, 'unif', 0);
disp(dataPath_exps);

inputFullnames = cell(numel(dataPaths), 1);
outputFullnames = cell(numel(dataPaths), 1);
if multiRegions
    imSizes = cell(numel(dataPaths), 1);
end

for i = 1 : numel(dataPath_exps)
    dataPath_i = dataPath_exps{i};
    if ~exist(dataPath_i, 'dir')
        continue;
    end
    
    dir_info = dir([dataPath_i, '*tif']);
    fsn = {dir_info.name}';
    inputFullnames{i} = cellfun(@(x) [dataPath_i, x], fsn, 'unif', 0);

    if multiRegions
        imSizes{i} =  cellfun(@(x) getImageSize([dataPath_i, x]), fsn, 'unif', 0);
    end
    
    fscPath = [dataPath_i, resultDirName, '/'];
    mkdir(fscPath);
    outputFullnames{i} = cellfun(@(x) [fscPath, x(1 : end - 4), '.mat'], fsn, 'unif', 0);
end

inputFullnames = cat(1, inputFullnames{:});
outputFullnames = cat(1, outputFullnames{:});
if multiRegions
    imSizes = cat(1, imSizes{:});
end

% filter filenames by channel patterns
include_flag = false(numel(inputFullnames), 1);
for c = 1 : numel(channelPatterns)
    include_flag = include_flag | contains(inputFullnames, channelPatterns{c}) | contains(inputFullnames, regexpPattern(channelPatterns{c}));
end
inputFullnames = inputFullnames(include_flag);
outputFullnames = outputFullnames(include_flag);
nF = numel(inputFullnames);

% process bbox for multiple regions
bboxes = repmat({inputBbox}, nF, 1);
if multiRegions
    imSizes = imSizes(include_flag);
    imSizes = cat(1, imSizes{:});
    bboxes = cell(nF, 1);
    if ~isempty(regionInterval) && isempty(regionGrid)
        cs = round((imSizes + 1) / 2);
        d = regionInterval;

        for f = 1 : nF
            ycs = cs(1) : -d(1) : halfSize + 1;
            ycs = [ycs(end : -1 : 2), cs(1) : d(1) : imSizes(f, 1) - halfSize + 1];
            if d(1) == -1
                ycs = cs(f, 1);
            end

            xcs = cs(2) : -d(2) : halfSize + 1;
            xcs = [xcs(end : -1 : 2), cs(2) : d(2) : imSizes(f, 2) - halfSize + 1];
            if d(2) == -1
                xcs = cs(f, 2);
            end

            zcs = cs(3) : -d(3) : halfSize + 1;
            zcs = [zcs(end : -1 : 2), cs(3) : d(3) : imSizes(f, 3) - halfSize + 1];
            if d(3) == -1
                zcs = cs(f, 3);
            end
            
            [Y, X, Z] = ndgrid(ycs, xcs, zcs);

            bbox_i = [[Y(:), X(:), Z(:)] - halfSize, [Y(:), X(:), Z(:)] + halfSize - 1];
            if any(bbox_i(:, 1 : 3) < 1 & bbox_i(:, 4 : 6) > imSizes(f, :), 'all')
                bbox_i(:, 1 : 3) = bbox_i(:, 1 : 3) .* (bbox_i(:, 1 : 3) >= 1) + (bbox_i(:, 1 : 3) < 1);
                bbox_i(:, 4 : 6) = bbox_i(:, 4 : 6) .* (bbox_i(:, 4 : 6) <= imSizes(f, :)) ...
                    + (bbox_i(:, 4 : 6) > imSizes(f, :)) .* imSizes(f, :);
            end
            bboxes{f} = bbox_i;
        end
    end

    if ~isempty(regionGrid)
        for f = 1 : nF
            bbox_i = [regionGrid - halfSize, regionGrid + halfSize - 1];
            if any(bbox_i(:, 1 : 3) < 1 & bbox_i(:, 4 : 6) > imSizes(f, :), 'all')
                bbox_i(:, 1 : 3) = bbox_i(:, 1 : 3) .* (bbox_i(:, 1 : 3) >= 1) + (bbox_i(:, 1 : 3) < 1);
                bbox_i(:, 4 : 6) = bbox_i(:, 4 : 6) .* (bbox_i(:, 4 : 6) <= imSizes(f, :)) ...
                    + (bbox_i(:, 4 : 6) > imSizes(f, :)) .* imSizes(f, :);
            end
            
            bboxes{f} = bbox_i;
        end
    end
end

func_strs = arrayfun(@(x) sprintf(['XR_one_image_FSC_analysis_frame(''%s'',''%s'',', ...
                '''xyPixelSize'',%0.20f,''dz'',%0.20f,''dr'',%0.20f,''dtheta'',%0.20f,', ...
                '''halfSize'',%s,''resThreshMethod'',''%s'',''resThresh'',%0.20f,''resAxis'',''%s'',', ...
                '''skipConeRegion'',%s,''bbox'',%s,''clipPer'',[%0.20f])'], inputFullnames{x}, ...
                outputFullnames{x}, xyPixelSize, dz, dr, dtheta, strrep(mat2str(halfSize), ' ', ','), ...
                resThreshMethod, resThresh, resAxis, string(skipConeRegion), ...
                strrep(mat2str(bboxes{x}), ' ', ','), clipPer), ...
                1 : numel(inputFullnames), 'unif', 0);

memAllocate = prod(halfSize) * 4 * 200 / 1024^3;
generic_computing_frameworks_wrapper(inputFullnames, outputFullnames, func_strs, ...
    cpusPerTask=cpusPerTask, memAllocate=memAllocate, parseCluster=parseCluster, masterCompute=masterCompute, ...
    mccMode=mccMode, configFile=configFile);


%% visualize FSC results

iter = iterInterval;
nC = numel(channelPatterns);
ntheta = round(2 * pi / dtheta) + 1;

if ~multiRegions
    for d = 1 : numel(dataPath_exps)
        dataPath_exp_d = dataPath_exps{d};
    
        dir_info = dir([dataPath_exp_d, '*tif']);
        fsns = {dir_info.name}';
        inputFullnames = cellfun(@(x) [dataPath_exp_d, x], {dir_info.name}', 'unif', 0);
        
        figurePath = [dataPath_exp_d, 'figures/'];
        if ~exist(figurePath, 'dir')
            mkdir(figurePath);
        end
    
        fscPath = [dataPath_exp_d, resultDirName, '/'];
    
        for c = 1 : nC
            include_flag = contains(inputFullnames, channelPatterns{c}) | contains(inputFullnames, regexpPattern(channelPatterns{c}));
            fsns_c = fsns(include_flag);
            if isempty(fsns_c)
                continue;
            end
            
            % res_decon = cell(numel(fsns_c), 1);
            res_mat = zeros(numel(fsns_c), ntheta);    
            for i = 1 : numel(fsns_c)
                % fn = [dataPath, fsns{i}];
                fscFn = [fscPath, fsns_c{i}(1:end-4), '.mat'];
                a = load(fscFn);
                
                % res_decon{i} = a.res_mu;
                res_mat(i, :) = a.res_mu.res;
            end
    
            res_mu_mat = mean(res_mat(:, 1 : end - 1, :), 2);
            thetas = a.res_mu.thetas;
    
            % save the mean resolutions to disk
            resFn = sprintf('%s/fsc_res_info_c%d.mat', fscPath, c);
            save('-v7.3', resFn, 'iter', 'thetas', 'res_mu_mat', 'res_mat', 'fsns_c');
            
            % xruan (04/14/2022): fill missing for inf values
            res_mat(isinf(res_mat)) = nan;
            res_mat = fillmissing(res_mat, 'linear', 2);
            res_mu_mat = mean(res_mat(:, 1 : end - 1), 2);
    
            figure('Renderer', 'painters', 'Position', [10 10 1200 800]);
            set(gcf, 'color', 'w')
            set(gcf, 'visible', 'off')
            plot(iter : iter : (numel(res_mu_mat)) * iter, res_mu_mat)
    
            [min_res, mind] = min(res_mu_mat);
            hold on, plot((mind) * iter, res_mu_mat(mind), 'o');
            
            % plot the line of 1.01 * min_res (to better determine the optimal iterations)
            hold on, plot(0 : iter : (numel(res_mu_mat)) * iter, ones(1, numel(res_mu_mat) + 1) * min_res * 1.01, '--');
    
            % annotation('textarrow', [0.3 0.5], [0.6, 0.5], 'String', sprintf('Iter = %d', mind -1));
            text(mind * iter, min(res_mu_mat(mind) * 1.3, res_mu_mat(mind) * 0.9 + max(res_mu_mat) * 0.1), sprintf('Iter = %d', (mind) * iter))
    
            grid on
            xlabel('Iteration')
            ylabel('Resolution (um)');
    
            f0 = gcf();
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s.png', figurePath, channels(c), suffix);
            print(f0, '-painters','-dpng', '-loose', figureFullname);
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s.fig', figurePath, channels(c), suffix);
            saveas(f0, figureFullname);
    
            % rescale y axis 
            ylim([0.2, 0.5])
            yticks(0.2 : 0.025 : 0.5)
    
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s_res_0p2_0p5.png', figurePath, channels(c), suffix);
            print(f0, '-painters','-dpng', '-loose', figureFullname);
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s_res_0p2_0p5.fig', figurePath, channels(c), suffix);
            saveas(f0, figureFullname);
            
            close(f0);
    
            % visualize resolution for different angles for each iteration
            if ~false
                for i = 1 : size(res_mat, 1)
    
                    figure('Renderer', 'painters', 'Position', [10 10 1200 800]);
                    set(gcf, 'color', 'w')
                    set(gcf, 'visible', 'off')
    
                    polarplot(a.res_mu.thetas, res_mat(i, :));
                    grid on
                    rlim([0, 2]);
                    title(sprintf('Iter %d', i * iter))
    
                    f0 = gcf();
                    figureFullname = sprintf('%s/fsc_res_vs_thetas_%s_%s.png', figurePath, fsns_c{i}(1 : end - 4), suffix);
                    print(f0, '-painters','-dpng', '-loose', figureFullname);
                    close(f0);
                end
            end
        end
    end
else
    for d = 1 : numel(dataPath_exps)
        dataPath_exp_d = dataPath_exps{d};
    
        dir_info = dir([dataPath_exp_d, '*tif']);
        fsns = {dir_info.name}';
        inputFullnames = cellfun(@(x) [dataPath_exp_d, x], {dir_info.name}', 'unif', 0);
        
        figurePath = [dataPath_exp_d, 'figures/'];
        if ~exist(figurePath, 'dir')
            mkdir(figurePath);
        end
    
        fscPath = [dataPath_exp_d, resultDirName, '/'];
    
        for c = 1 : nC
            include_flag = contains(inputFullnames, channelPatterns{c}) | contains(inputFullnames, regexpPattern(channelPatterns{c}));
            fsns_c = fsns(include_flag);
            if isempty(fsns_c)
                continue;
            end
            
            % res_decon = cell(numel(fsns_c), 1);
            fscFn = [fscPath, fsns_c{1}(1:end-4), '.mat'];
            a = load(fscFn);
            % number of regions
            nR = size(a.bbox, 1);
            res_mat = zeros(numel(fsns_c), ntheta, nR);
    
            for i = 1 : numel(fsns_c)
                % fn = [dataPath, fsns{i}];
                fscFn = [fscPath, fsns_c{i}(1:end-4), '.mat'];
                a = load(fscFn);
                
                % res_decon{i} = a.res_mu;
                for j = 1 : nR
                    res_mat(i, :, j) = a.res_mu{j}.res;
                end
            end
    
            res_mu_mat = mean(res_mat(:, 1 : end - 1, :), 2);
            thetas = a.res_mu{1}.thetas;
    
            % save the mean resolutions to disk
            resFn = sprintf('%s/fsc_res_info_c%d.mat', fscPath, c);
            save('-v7.3', resFn, 'iter', 'thetas', 'res_mu_mat', 'res_mat', 'fsns_c');
            
            % xruan (04/14/2022): fill missing for inf values
            res_mat(isinf(res_mat)) = nan;
            res_mat = fillmissing(res_mat, 'linear', 2);
            res_mu_mat = squeeze(mean(res_mat(:, 1 : end - 1, :), 2));
            res_overall_mu_mat = mean(res_mu_mat, 2);
    
            figure('Renderer', 'painters', 'Position', [10 10 1200 800]);
            set(gcf, 'color', 'w')
            set(gcf, 'visible', 'off')
            % hold on, plot(0 : 200, mean(res_mat_double(:, 1 : end - 1), 2))
            hold on, plot(iter : iter : (numel(fsns_c)) * iter, res_mu_mat, '--')
            hold on, plot(iter : iter : (numel(fsns_c)) * iter, res_overall_mu_mat)

            [~, mind] = min(res_overall_mu_mat);
            hold on, plot((mind) * iter, res_overall_mu_mat(mind), 'o');

            % point for mean + std (minimum point)
            res_overall_std_mat = std(res_mu_mat, [], 2);
            [~, ms_mind] = min(abs(res_overall_mu_mat(mind : end) - (res_overall_mu_mat(mind) + res_overall_std_mat(mind))));
            ms_mind = ms_mind + mind - 1;
    
            % annotation('textarrow', [0.3 0.5], [0.6, 0.5], 'String', sprintf('Iter = %d', mind -1));
            text(mind * iter, min(res_overall_mu_mat(mind) * 1.3, res_overall_mu_mat(mind) * 0.9 + max(res_overall_mu_mat) * 0.1), sprintf('Iter = %d', (mind) * iter))
            text(ms_mind * iter, min(res_overall_mu_mat(ms_mind) * 1.3, res_overall_mu_mat(ms_mind) * 0.9 + max(res_overall_mu_mat) * 0.1), sprintf('Iter = %d (mu+std)', (ms_mind) * iter))
    
            grid on
            xlabel('Iteration')
            ylabel('Resolution (um)');
    
            f0 = gcf();
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s.png', figurePath, channels(c), suffix);
            print(f0, '-vector','-dpng', '-loose', figureFullname);
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s.fig', figurePath, channels(c), suffix);
            saveas(f0, figureFullname);
    
            % rescale y axis 
            ylim([0.2, 0.5])
            yticks(0.2 : 0.025 : 0.5)
    
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s_res_0p2_0p5.png', figurePath, channels(c), suffix);
            print(f0, '-vector','-dpng', '-loose', figureFullname);
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s_res_0p2_0p5.fig', figurePath, channels(c), suffix);
            saveas(f0, figureFullname);
            
            close(f0);
    
            % visualize resolution for different angles for each iteration
            if ~false
                for i = 1 : size(res_mat, 1)
    
                    figure('Renderer', 'painters', 'Position', [10 10 1200 800]);
                    set(gcf, 'color', 'w')
                    set(gcf, 'visible', 'off')
    
                    polarplot(thetas, mean(res_mat(i, :, :), 3));
                    grid on
                    rlim([0, 2]);
                    title(sprintf('Iter %d', i * iter))
    
                    f0 = gcf();
                    figureFullname = sprintf('%s/fsc_res_vs_thetas_%s_%s.png', figurePath, fsns_c{i}(1 : end - 4), suffix);
                    print(f0, '-vector','-dpng', '-loose', figureFullname);
                    close(f0);
                end
            end
        end
    end
end


end



