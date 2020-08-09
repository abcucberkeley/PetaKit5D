%plotLifetimes(lftRes, varargin) displays CCP lifetime distributions
%
% Input: 
%         lftRes : structure returned by runLifetimeAnalysis()
%
% Options ('specifier', value):
%      'PlotAll' : true|{false} detailed display, includes CCP, CS, and visitor distributions

% Francois Aguet (last modified 05/13/2013)

function h = plotLifetimes(lftRes, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('DisplayMode', '');
ip.addParamValue('ShowExpFits', false, @islogical);
ip.addParamValue('ShowStatistics', false, @islogical);
ip.addParamValue('ShowCargoDependent', true, @islogical);
ip.addParamValue('SlaveNames', '');
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('Hues', []);
ip.addParamValue('XTick', 0:20:120);
ip.addParamValue('YTick', 0:0.01:0.04);
ip.addParamValue('SingleChannel', false);
ip.parse(varargin{:});
ya = ip.Results.YTick;
lw = 0.75;
h = [];
chNames = ip.Results.SlaveNames;
if isempty(chNames)
    chNames = '';
end
if ~iscell(chNames)
    chNames = {chNames};
end

fmt = ['%.' num2str(ceil(abs(log10(ya(2))))) 'f'];
yal = ['0' arrayfun(@(x) num2str(x, fmt), ya(2:end), 'UniformOutput', false)];

fset = loadFigureSettings(ip.Results.DisplayMode);
ce = [0 0 0.6;
    1/3 0.3 0.9;
    1/3 1 0.9];
ce = hsv2rgb(ce);

if isstruct(lftRes)
    %============================================================
    % Single-channel data
    %============================================================
    if ~isfield(lftRes, 'lftHistSlaveCCP') || ip.Results.SingleChannel
        h(1) = setupFigure('DisplayMode', ip.Results.DisplayMode, 'Name', 'Lifetime distr.');
        
        if ip.Results.PlotAll
            if isfield(lftRes, 'meanLftHistVisit')
                hp = zeros(4,1);
                hp(4) = plot(lftRes.t, mean(lftRes.pctVisit)*lftRes.meanLftHistVisit, '-', 'Color', fset.cfB, 'LineWidth', lw);
            else
                hp = zeros(3,1);
            end
            hp(3) = plot(lftRes.t, mean(vertcat(lftRes.lftHist_Ia), 1), 'Color', ce(1,:), 'LineWidth', lw);
            hp(2) = plot(lftRes.t, mean(lftRes.pctCS)*lftRes.meanLftHistCS, '-', 'Color', ce(2,:), 'LineWidth', lw);
            if ip.Results.ShowExpFits
                ff = mean(lftRes.pctCS)*lftRes.meanLftHistCS;
                [mu,~,Aexp,~] = fitExpToHist(lftRes.t(5:end), ff(5:end));
                tx = 0:0.1:lftRes.t(end);
                plot(tx, Aexp/mu*exp(-1/mu*tx), 'r--', 'LineWidth', 1)
            end
            hp(1) = plot(lftRes.t, mean(lftRes.pctCCP)*lftRes.meanLftHistCCP, '-', 'Color', ce(3,:), 'LineWidth', lw+0.5);
            
            % Legend:
            ltext = {[' CCPs: ' num2str(mean(lftRes.pctCCP)*100, '%.1f') ' ± ' num2str(std(lftRes.pctCCP)*100, '%.1f') ' %'],...
                [' CSs: ' num2str(mean(lftRes.pctCS)*100, '%.1f') ' ± ' num2str(std(lftRes.pctCS)*100, '%.1f') ' %'],...
                ' All structures'};
            lheight = 1.25;
            if isfield(lftRes, 'meanLftHistVisit')
                ltext = [ltext(1:2) [' Visitors: ' num2str(mean(lftRes.pctVisit)*100, '%.1f') ' ± ' num2str(std(lftRes.pctVisit)*100, '%.1f') ' %'] ltext(3)];
                hp = hp([1 2 4 3]);
                lheight = 1.5;
            end
            hl = legend(hp, ltext{:});
            set(hl, 'Box', 'off');
            if strcmpi(ip.Results.DisplayMode, 'print')
                set(hl, 'Units', 'centimeters', 'Position', [4 5-lheight 1.75 lheight]);
            end
            ylabel('Relative frequency', fset.lfont{:});
        else
            plot(lftRes.t, lftRes.meanLftHistCCP, '-', 'Color', ce(3,:), 'LineWidth', lw+0.5);
            hl = legend(' All CCPs');
            set(hl, 'Box', 'off');
            if strcmpi(ip.Results.DisplayMode, 'print')
                set(hl, 'Position', [5 4.25 1.75 0.75]);
            end
            
            ylabel('Frequency', fset.lfont{:});
        end
        axis([0 min(ip.Results.XTick(end), lftRes.t(end)) 0 ya(end)]);
        set(gca, 'XTick', ip.Results.XTick, 'YTick', ya, 'YTickLabel', yal);
        xlabel('Lifetime (s)', fset.lfont{:});
        
        %if ip.Results.ShowStatistics
        %    fs = loadFigureSettings('print');
        %    h(2) = figure(fs.fOpts{:}, 'Position', [5 5 5 6.5]);
        %    axes(fs.axOpts{:}, 'Position', [1.5 2 3 4]);
        %    boxplot2(lftRes.stats, 'AdjustFigure', false, 'XLabels', {'Raw', 'Below thresh', 'Above thresh'},...
        %        'FaceColor', ce, 'BarWidth', 0.6, 'LineWidth', 1);
        %    ylabel('Lifetime (s)', fset.sfont{:});
        %end
    end
    
    %============================================================
    % Multi-channel data
    %============================================================
    if ip.Results.ShowCargoDependent && isfield(lftRes, 'lftHistSlaveCCP') && ~ip.Results.SingleChannel
        
        % 1) plot all combinations
        %framerate = lftRes.t(2)-lftRes.t(1);
        framerate = 1; % TO DO: add normalization option
        
        ncomb = size(lftRes.slaveCombs,1);
        
        % pct CCPs/CSs (above/below)
        pctCCP = mean(lftRes.pctSlaveCCP,1);
        pctCCPStd = std(lftRes.pctSlaveCCP, [], 1);
        pctCS = mean(lftRes.pctSlaveCS,1);
        
        tmp = double(lftRes.slaveCombs);
        tmp(tmp==1) = '+';
        tmp(tmp==0) = '-';
        labelsA = cell(1,ncomb);
        labelsB = cell(1,ncomb);
        labelsC = cell(1,ncomb);
        for s = 1:ncomb
            labelsA{s} = [' ' tmp(s,1) ' ' chNames{1}];
            labelsB{s} = labelsA{s};
            for c = 2:numel(chNames)
                labelsA{s} = [labelsA{s} '/' tmp(s,c) ' ' chNames{c}];
                labelsB{s} = [labelsB{s} '/' tmp(s,c) ' ' chNames{c}];
            end
            labelsC{s} = [labelsA{s} ' (' num2str(pctCCP(s)/sum(pctCCP)*100, '%.1f') '±' num2str(pctCCPStd(s)/sum(pctCCP)*100, '%.1f') '%)'];
            labelsA{s} = [labelsA{s} ' CCPs (' num2str(pctCCP(s)*100, '%.1f') '%)'];
            labelsB{s} = [labelsB{s} ' CSs (' num2str(pctCS(s)*100, '%.1f') '%)'];
        end
        
        switch ncomb
            case 2
                cmap = [0.33 1 0.9;
                    0.6  1 0.9;
                    0.23 0.7 0.9;
                    0.5  0.7 0.9];
            case 4
                cmap = [0    1 0.9;
                    0.33 1 0.9;
                    0.55 1 0.9;
                    0    0 0.6;
                    0    0.7 0.9;
                    0.28 0.7 0.9;
                    0.5 0.7 0.9;
                    0    0.7 0];
        end
        cmap = hsv2rgb(cmap);
        
        %------------------------------------------------------------
        % 1) Lifetimes distributions for all detected objects
        %------------------------------------------------------------
        if ip.Results.PlotAll
            setupFigure('Name', 'Lifetime dist.', 'DisplayMode', ip.Results.DisplayMode);
            axis([0 min(ip.Results.XTick(end), lftRes.t(end)) 0 ya(end)]);
            
            hp = zeros(1+2*ncomb,1);
            %bAll = 1.96*getSE(lftRes, 'lftHist_Ia');
            %fill([lftRes.t lftRes.t(end:-1:1)], [mu+bAll mu(end:-1:1)-bAll(end:-1:1)], 'r', 'EdgeColor', 'none');
            hp(1) = plot(lftRes.t, mean(vertcat(lftRes.lftHist_Ia), 1)*framerate, 'Color', 0.6*[1 1 1], 'LineWidth', lw);
            for s = ncomb:-1:1 % plot combinations in increasing order of association
                hp(2*(s-1)+3) = plot(lftRes.t, pctCS(s)*mean(lftRes.lftHistSlaveCS{s},1)*framerate, 'Color', cmap(s+ncomb,:), 'LineWidth', lw);
            end
            for s = ncomb:-1:1 % plot combinations in increasing order of association
                hp(2*(s-1)+2) = plot(lftRes.t, pctCCP(s)*mean(lftRes.lftHistSlaveCCP{s},1)*framerate, 'Color', cmap(s,:), 'LineWidth', lw+0.5);
                % missing: histograms of slave classification for all structures (lftRes.lftHistAll{s})
            end
            
            set(gca, 'XTick', ip.Results.XTick, 'YTick', ya, 'YTickLabel', yal);
            xlabel('Lifetime (s)', fset.lfont{:});
            ylabel('Relative frequency', fset.lfont{:});
            
            labels = [labelsA labelsB];
            hl = legend(hp, [' All structures' labels(1:2:end) labels(2:2:end)]);
            set(hl, fset.tfont{:}, 'Box', 'off');
            if strcmpi(ip.Results.DisplayMode, 'print');
                set(hl, 'Units', 'centimeters', 'Position', [3.75 3.25 2 2]);
            end
        end
        
        %------------------------------------------------------------
        % 2) Plot CCP lifetimes only; all combinations
        %------------------------------------------------------------
        setupFigure('Name', 'Lifetime dist.', 'DisplayMode', ip.Results.DisplayMode);
        
        hp = zeros(1,ncomb);
        %tmp = arrayfun(@(s) pctCCP(s)/sum(pctCCP)*mean(lftRes.lftHistSlaveCCP{s},1)*framerate, 1:ncomb, 'unif', 0);
        %tmp = cat(1,tmp{:});
        %hp(1) = plot(lftRes.t, sum(tmp,1), 'k', 'LineWidth', lw);
        hp(1) = plot(lftRes.t, lftRes.meanLftHistCCP, 'k', 'LineWidth', lw);
        for s = ncomb:-1:1 % plot combinations in increasing order of association
            hp(s+1) = plot(lftRes.t, pctCCP(s)/sum(pctCCP)*mean(lftRes.lftHistSlaveCCP{s},1), 'Color', cmap(s,:), 'LineWidth', lw+0.5);
        end
        
        axis([0 min(ip.Results.XTick(end), lftRes.t(end)) 0 ya(end)]);
        set(gca, 'XTick', ip.Results.XTick, 'YTick', ya, 'YTickLabel', yal);
        
        xlabel('Lifetime (s)', fset.lfont{:});
        ylabel('Relative frequency', fset.lfont{:});
        
        hl = legend(hp, [' All CCPs' labelsC], 'Location', 'NorthEast');
        lheight = ncomb+1;%1.5+3.5 = 5 -> 4+1.2
        set(hl, 'Box', 'off');
        if strcmpi(ip.Results.DisplayMode, 'print');
            set(hl, fset.sfont{:}, 'Units', 'centimeters', 'Position', [3.2 5.2-lheight*0.35 2.5 lheight*0.35]);
        end
    end
    
    % for comparisons between multiple conditions
elseif iscell(lftRes)
    if nargin<2
        colorA = hsv2rgb([0.6 1 1;
            1/3 1 1;
            0 1 1;
            1/7 1 1;
            5/9 1 1]);
    end
    colorB = rgb2hsv(colorA);
    colorB(:,2) = 0.5;
    colorB = hsv2rgb(colorB);
    
    nd = numel(lftRes);
    h = figure;
    hold on;
    for i = 1:nd
        hp(2*(i-1)+2) = plot(lftRes{i}.t, lftRes{i}.meanLftHist_B, '.-', 'Color', colorB(i,:), 'LineWidth', 2, 'MarkerSize', 16);
        hp(2*(i-1)+1) = plot(lftRes{i}.t, lftRes{i}.meanLftHist_A, '.-', 'Color', colorA(i,:), 'LineWidth', 2, 'MarkerSize', 16);
        expName = getDirFromPath(getExpDir(lftRes{i}.data));
        legendText{2*(i-1)+1} = [expName ', above threshold (' num2str(mean(lftRes{i}.pctAbove)*100,'%.1f') ' ± ' num2str(std(lftRes{i}.pctAbove)*100,'%.1f') ' %)'];
        legendText{2*(i-1)+2} = [expName ', below threshold (' num2str(mean(1-lftRes{i}.pctAbove)*100,'%.1f') ' ± ' num2str(std(lftRes{i}.pctAbove)*100,'%.1f') ' %)'];
    end
    axis([0 min(ip.Results.XTick(end), lftRes{1}.t(end)) 0 ya(end)]);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    hl = legend(hp, legendText{:}, 'Location', 'NorthEast');
    set(hl, 'Box', 'off', fset.ifont{:}, 'Interpreter', 'none');
    
else
    error('Incompatible input');
end

function SE = getSE(lftData, fieldName)
M = lftData.(fieldName);
nd = size(M,1);
meanM = zeros(size(M));
for i = 1:nd
    meanM(i,:) = mean(M(setdiff(1:nd,i),:),1);
end
SE = std(meanM,[],1) / sqrt(nd);
