function XR_FSC_analysis_wrapper(dataPaths, varargin)
% wrapper for datasets level FSC analysis, visualization of mean resolution
% resolution output: 0 for x, pi/2 for z.
%
% Author: Xiongtao Ruan (12/10/2021)
% 
% xruan (05/12/2022): add support for multi-region FSCs
% xruan (06/09/2022): add support for clipping very bright spots

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths');
ip.addParameter('outDirstr', 'FSCs', @ischar);
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric);
% ip.addParameter('angle', 32.45, @isnumeric);
ip.addParameter('dr', 1 , @isnumeric);
ip.addParameter('dtheta', pi / 12 , @isnumeric);
ip.addParameter('resThreshMethod', 'fixed', @ischar);
ip.addParameter('resThresh', 0.2, @isnumeric);
ip.addParameter('N', [251, 251, 251], @isnumeric);
ip.addParameter('bbox', [], @isnumeric);
ip.addParameter('resAxis', 'xz', @ischar);
ip.addParameter('skipConeRegion', true, @islogical);
% ip.addParameter('Deskew', true, @islogical);
% ip.addParameter('flipZstack', false, @islogical);
% ip.addParameter('ObjectiveScan', false, @islogical);
% ip.addParameter('ZstageScan', false, @islogical);
ip.addParameter('ChannelPatterns', {'tif'}, @iscell);
ip.addParameter('Channels', [488, 560], @isnumeric);
ip.addParameter('multiRegions', false, @islogical);
ip.addParameter('regionInterval', [50, 50, -1], @isnumeric); % yxz, -1 means only center
ip.addParameter('regionGrid', [], @isnumeric); % user provided grid for region centers, N x 3
ip.addParameter('clipPer', [], @isnumeric); % clip intensity higher than the given percentile
ip.addParameter('suffix', 'decon', @ischar);
ip.addParameter('iterInterval', 5, @isnumeric); % iteration interval for FSC resolution plot
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.parse(dataPaths, varargin{:});

pr = ip.Results;
outDirstr = pr.outDirstr;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
% angle = pr.angle;
dr = pr.dr;
dtheta = pr.dtheta;
N = pr.N;
bbox = pr.bbox;
resThreshMethod = pr.resThreshMethod;
resThresh = pr.resThresh;
resAxis = pr.resAxis;
skipConeRegion = pr.skipConeRegion;
ChannelPatterns = pr.ChannelPatterns;
Channels = pr.Channels;
multiRegions = pr.multiRegions;
regionInterval = pr.regionInterval;
regionGrid = pr.regionGrid;
clipPer = pr.clipPer;
suffix = pr.suffix;
iterInterval = pr.iterInterval;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;


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
    
    fscPath = [dataPath_i, outDirstr, '/'];
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
for c = 1 : numel(ChannelPatterns)
    include_flag = include_flag | contains(inputFullnames, ChannelPatterns{c}) | contains(inputFullnames, regexpPattern(ChannelPatterns{c}));
end
inputFullnames = inputFullnames(include_flag);
outputFullnames = outputFullnames(include_flag);
nF = numel(inputFullnames);

% process bbox for multiple regions
bboxes = repmat({bbox}, nF, 1);
if multiRegions
    imSizes = imSizes(include_flag);
    imSizes = cat(1, imSizes{:});
    bboxes = cell(nF, 1);
    if ~isempty(regionInterval) && isempty(regionGrid)
        cs = round((imSizes + 1) / 2);
        d = regionInterval;

        for f = 1 : nF
            ycs = cs(1) : -d(1) : N + 1;
            ycs = [ycs(end : -1 : 2), cs(1) : d(1) : imSizes(f, 1) - N + 1];
            if d(1) == -1
                ycs = cs(f, 1);
            end

            xcs = cs(2) : -d(2) : N + 1;
            xcs = [xcs(end : -1 : 2), cs(2) : d(2) : imSizes(f, 2) - N + 1];
            if d(2) == -1
                xcs = cs(f, 2);
            end

            zcs = cs(3) : -d(3) : N + 1;
            zcs = [zcs(end : -1 : 2), cs(3) : d(3) : imSizes(f, 3) - N + 1];
            if d(3) == -1
                zcs = cs(f, 3);
            end
            
            [Y, X, Z] = ndgrid(ycs, xcs, zcs);

            bbox_i = [[Y(:), X(:), Z(:)] - N, [Y(:), X(:), Z(:)] + N - 1];
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
            bbox_i = [regionGrid - N, regionGrid + N - 1];
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
                '''N'',%s,''resThreshMethod'',''%s'',''resThresh'',%0.20f,''resAxis'',''%s'',', ...
                '''skipConeRegion'',%s,''bbox'',%s,''clipPer'',[%0.20f]);'], inputFullnames{x}, ...
                outputFullnames{x}, xyPixelSize, dz, dr, dtheta, strrep(mat2str(N), ' ', ','), ...
                resThreshMethod, resThresh, resAxis, string(skipConeRegion), ...
                strrep(mat2str(bboxes{x}), ' ', ','), clipPer), ...
                1 : numel(inputFullnames), 'unif', 0);

slurm_cluster_generic_computing_wrapper(inputFullnames, outputFullnames, func_strs, ...
    'cpusPerTask', 4, 'cpuOnlyNodes', false, 'parseCluster', parseCluster, 'masterCompute', masterCompute);


%% visualize FSC results

iter = iterInterval;
nC = numel(ChannelPatterns);
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
    
        fscPath = [dataPath_exp_d, outDirstr, '/'];
    
        for c = 1 : nC
            include_flag = contains(inputFullnames, ChannelPatterns{c}) | contains(inputFullnames, regexpPattern(ChannelPatterns{c}));
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
            % hold on, plot(0 : 200, mean(res_mat_double(:, 1 : end - 1), 2))
    
            [~, mind] = min(res_mu_mat);
            hold on, plot((mind) * iter, res_mu_mat(mind), 'o');
    
            % annotation('textarrow', [0.3 0.5], [0.6, 0.5], 'String', sprintf('Iter = %d', mind -1));
            text(mind * iter, min(res_mu_mat(mind) * 1.3, res_mu_mat(mind) * 0.9 + max(res_mu_mat) * 0.1), sprintf('Iter = %d', (mind) * iter))
    
            grid on
            xlabel('Iteration')
            ylabel('Resolution (um)');
    
            f0 = gcf();
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s.png', figurePath, Channels(c), suffix);
            print(f0, '-painters','-dpng', '-loose', figureFullname);
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s.fig', figurePath, Channels(c), suffix);
            saveas(f0, figureFullname);
    
            % rescale y axis 
            ylim([0.2, 0.5])
            yticks(0.2 : 0.025 : 0.5)
    
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s_res_0p2_0p5.png', figurePath, Channels(c), suffix);
            print(f0, '-painters','-dpng', '-loose', figureFullname);
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s_res_0p2_0p5.fig', figurePath, Channels(c), suffix);
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
    
        fscPath = [dataPath_exp_d, outDirstr, '/'];
    
        for c = 1 : nC
            include_flag = contains(inputFullnames, ChannelPatterns{c}) | contains(inputFullnames, regexpPattern(ChannelPatterns{c}));
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
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s.png', figurePath, Channels(c), suffix);
            print(f0, '-painters','-dpng', '-loose', figureFullname);
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s.fig', figurePath, Channels(c), suffix);
            saveas(f0, figureFullname);
    
            % rescale y axis 
            ylim([0.2, 0.5])
            yticks(0.2 : 0.025 : 0.5)
    
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s_res_0p2_0p5.png', figurePath, Channels(c), suffix);
            print(f0, '-painters','-dpng', '-loose', figureFullname);
            figureFullname = sprintf('%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s_res_0p2_0p5.fig', figurePath, Channels(c), suffix);
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
                    print(f0, '-painters','-dpng', '-loose', figureFullname);
                    close(f0);
                end
            end
        end
    end
end


end



