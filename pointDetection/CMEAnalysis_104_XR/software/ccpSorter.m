%ccpSorter(data) displays a GUI for selection of CCP cohorts based on maximum intensity thresholds
%
% Inputs:
%             data : data structure from loadConditionData()
%
% Parameters:
%      
%  'DetectionMode' : 'm' selects tracks with independently significant fluorescence signal
%                        in the second channel.
%                    's' selects tracks with signal that is significant relative to background
%                        but requires a master channel for detection. Default. 
%
% Example:
%   ccpSorter(data, 'DetectionMode', 'm');
%
% Copyright (C) 2017, Danuser Lab - UTSouthwestern 
%
% This file is part of CMEAnalysis_Package.
% 
% CMEAnalysis_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CMEAnalysis_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with CMEAnalysis_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Francois Aguet, 10/25/2013

function ccpSorter(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('LftDataName', 'lifetimeData.mat');
ip.addParamValue('Cutoff_f', 5);
ip.addParamValue('Rescale', true, @islogical);
ip.addParamValue('ExcludeVisitors', true, @islogical);
ip.addParamValue('DetectionMode', 's', @(x) any(strcmpi(x, {'m', 's'})));
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 120]);
ip.parse(data, varargin{:});
cohortBounds = ip.Results.CohortBounds_s;

% load data
lftData = getLifetimeData(data,...
    'LifetimeData', ip.Results.LftDataName, 'Scale', ip.Results.Rescale,...
    'Cutoff_f', ip.Results.Cutoff_f, 'ReturnValidOnly', true,...
    'ExcludeVisitors', ip.Results.ExcludeVisitors, 'DisplayScaling', false);

% For scatter plots (pooled data):
%  max. amplitude
maxA = arrayfun(@(i) squeeze(max(i.A,[],2)), lftData, 'unif', 0);
%  lifetime
lftV = arrayfun(@(i) i.lifetime_s, lftData, 'unif', 0);

framerate = data(1).framerate;

% intensity scaling for 2nd channel
chScale = [1 1];%[1 scaleEDFs(maxAall(:,2), maxAall(:,1))];

% calculate cohorts
% cohortBounds = [cohortBounds max(vertcat(lftData.lifetime_s))];
[res, cTime] = getIntensityCohorts(data, lftData, 'CohortBounds_s', cohortBounds);

nCh = numel(data(1).channels);

% apply selection
% if nCh>1
%     if strcmpi(ip.Results.DetectionMode, 'm')
%         for i = 1:numel(lftData)
%             idx = lftData(i).significantMaster(:,2)==1;
%             res(i).aInterp = res(i).aInterp(idx,:,:);
%             res(i).rInterp = res(i).rInterp(idx,:,:);
%             res(i).cidx = res(i).cidx(idx);
%             maxA{i} = maxA{i}(idx,:);
%             lftV{i} = lftV{i}(idx);
%         end
%     elseif strcmpi(ip.Results.DetectionMode, 's')
%         for i = 1:numel(lftData)
%             idx = lftData(i).significantSlave(:,2)==1;
%             res(i).aInterp = res(i).aInterp(idx,:,:);
%             res(i).rInterp = res(i).rInterp(idx,:,:);
%             res(i).cidx = res(i).cidx(idx);
%             maxA{i} = maxA{i}(idx,:);
%             lftV{i} = lftV{i}(idx);
%         end
%     end
% end



% intensity matrix for cohort plots
AMat = cat(1, res.aInterp); % order: long->short tracks
if nCh>1
    AMat(:,:,2) = chScale(2)*AMat(:,:,2);
end
% sort from shortest to longest cohort
cidx = cat(1, res.cidx);
% lftAMat = arrayfun(@(i) lftV{i}(res(i).cidxAll), 1:numel(res), 'unif', 0);
% lftAMat = vertcat(lftAMat{:});
lftAMat = cat(1,res.lifetime);
[cidx, idx] = sort(cidx);
AMat = AMat(idx,:,:);
maxAMat = squeeze(max(AMat,[],2));
lftAMat = lftAMat(idx);

% vectors for scatter plot
lftV = vertcat(lftV{:});
maxA = vertcat(maxA{:});
if nCh>1
    maxA(:,2) = chScale(2)*maxA(:,2);
else
    maxA = maxA(:,1);
end

nCh = numel(data(1).channels);
hues = getFluorophoreHues(data(1).markers);

% initial thresholds
lt0 = 10; % seconds
at0 = prctile(maxA, 10, 1);

% selection indexes for scatter plots
selCh1 = maxA(:,1)>=at0(1) & lftV>=lt0;
selColor1 = hsv2rgb([hues(1) 1 0.8]);
if nCh>1
    selCh2 = maxA(:,2)>=at0(2) & lftV>=lt0;
    selColor2 = hsv2rgb([hues(2) 1 0.8]);
end


% fpos = get(0, 'DefaultFigurePosition');
width = 1000;
height = 450;
figure('Position', [200 200 width height], 'Color', 'w', 'Name', 'ccpSorter',...
    'Toolbar', 'figure');

if nCh>1
    aRange = [0 0; prctile(maxA, 99.9, 1)]';
else
    aRange = [0 prctile(maxA, 99.9, 1)];
end
lRange = [0 prctile(lftV, 99.9, 1)];

% 1) Scatter plot(s)
ha(1) = axes('Position', [0.12/width*height 0.35 0.6/width*height 0.6],'FontSize', 8);
hp1 = zeros(2,1);
hp1(1) = plot(ha(1), lftV(~selCh1), maxA(~selCh1,1), 'o', 'Color', 0.6*[1 1 1], 'HitTest', 'off');
hold(ha(1), 'on');
hp1(2) = plot(ha(1), lftV(selCh1), maxA(selCh1,1), 'o', 'Color', selColor1, 'HitTest', 'off');
axis(ha(1), [lRange aRange(1,:)]);
% thresholds
ht = zeros(nCh,2);
ht(1,1) = plot(ha(1), lt0*[1 1], aRange(1,:), 'r', 'HitTest', 'off');
ht(1,2) = plot(ha(1), lRange, at0(1)*[1 1], 'g', 'HitTest', 'off');
xlabel(ha(1), 'Lifetime (s)', 'FontSize', 12);
ylabel(ha(1), 'Max. fluorescence intensity (A.U.)', 'FontSize', 12);
htxt(1) = text(lRange(1)+0.95*diff(lRange), aRange(1,1)+0.95*diff(aRange(1,:)), [num2str(sum(selCh1)) '/' num2str(numel(selCh1))], 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(mean(lRange), aRange(1,2), ['Channel 1: ' data(1).markers{1}], 'Parent', ha(1),...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
scattercontour(lftV, maxA(:,1), 'Parent', ha(1));

% for lifetime plot
bf = 1; % binning
buffer = size(lftData(1).sbA,2);
cutoff_f = ip.Results.Cutoff_f;
N = data(1).movieLength-2*buffer;
t = ((cutoff_f+(bf-1)/2):bf:N)*framerate;
w = N./(N-cutoff_f+1:-1:1);

if nCh>1
    ha(2) = axes('Position', [0.8/width*height 0.35 0.6/width*height 0.6],'FontSize', 8);
    hp2 = zeros(3,1);
    hp2(1) = plot(ha(2), lftV(~selCh2 & selCh1), maxA(~selCh2 & selCh1,2), 'o', 'Color', selColor1, 'HitTest', 'off');
    hold(ha(2), 'on');
    hp2(2) = plot(ha(2), lftV(selCh2 & selCh1), maxA(selCh2 & selCh1,2), 'o', 'Color', selColor2, 'HitTest', 'off');
    hp2(3) = plot(ha(2), lftV(~selCh1), maxA(~selCh1,2), 'o', 'Color', 0.6*[1 1 1], 'HitTest', 'off');
    axis(ha(2), [lRange aRange(2,:)]);
    ht(2,1) = plot(ha(2), lt0*[1 1], aRange(2,:), 'r', 'HitTest', 'off');
    ht(2,2) = plot(ha(2), lRange, at0(2)*[1 1], 'b', 'HitTest', 'off');
    xlabel(ha(2), 'Lifetime (s)', 'FontSize', 12);
    % ylabel(ha(2), 'Max. fluorescence intensity (A.U.)', 'FontSize', 12);
    htxt(2) = text(lRange(1)+0.95*diff(lRange), aRange(2,1)+0.95*diff(aRange(2,:)), [num2str(sum(selCh1&selCh2)) '/' num2str(sum(selCh1))], 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    text(mean(lRange), aRange(2,2), ['Channel 2: ' data(1).markers{2}], 'Parent', ha(2),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    scattercontour(lftV, maxA(:,2), 'Parent', ha(2));
else
    % Plot lifetime distribution
    ha(2) = axes('Position', [0.8/width*height 0.6 0.6/width*height 0.35],'FontSize', 8);
    
    histCCPs = hist(lftV(maxA(:,1)>=at0(1)), t).*w;
    histCCPs = histCCPs / sum(histCCPs) / framerate / bf;
    hp2(1) = plot(ha(2), t, histCCPs, 'g', 'HitTest', 'off');
    hold(ha(2), 'on');
    
    histCSs = hist(lftV(maxA(:,1)<at0(1)), t).*w;
    histCSs = histCSs / sum(histCSs) / framerate / bf;
    hp2(2) = plot(ha(2), t, histCSs, 'Color', hsv2rgb([0.55 1 0.9]), 'HitTest', 'off');
    
    axis(ha(2), [0 160 0 0.05]);
    xlabel(ha(2), 'Lifetime (s)', 'FontSize', 12);
    
    nCCP = sum(maxA(:,1)>=at0(1));
    nTotal = numel(maxA(:,1));
    hl = legend([' CCPs: ' num2str(nCCP/nTotal*100, '%.1f') '%'], [' CSs: ' num2str((nTotal-nCCP)/nTotal*100, '%.1f') '%']);
    set(hl, 'EdgeColor', 'w');
end
set(ha(1:2), 'ButtonDownFcn', @click_Callback); % after 'hold on'

ha(3) = axes('Position', [1.55/width*height 0.6 0.62/width*height 0.35]);
ha(4) = axes('Position', [1.55/width*height 0.2 0.62/width*height 0.35]);
% linkaxes(ha(3:4));
hold(ha(3), 'on');
hold(ha(4), 'on');


lftThresholdLabel = uicontrol('Style', 'text', 'String', ['Lifetime threshold: ' num2str(lt0, '%.0f') ' s'],...
    'Units', 'pixels', 'Position', [60 90 250 20], 'BackgroundColor', 'w',...
    'HorizontalAlignment', 'left', 'FontSize', 12,'FontUnits','pixels');

c1Label = uicontrol('Style', 'text', 'String', ['Channel 1 threshold: ' num2str(at0(1), '%.0f')],...
    'Units', 'pixels', 'Position', [60 70 250 20], 'BackgroundColor', 'w',...
    'HorizontalAlignment', 'left', 'FontSize', 12,'FontUnits','pixels');

if nCh>1
    c2Label = uicontrol('Style', 'text', 'String', ['Channel 2 threshold: ' num2str(at0(2), '%.0f')],...
        'Units', 'pixels', 'Position', [60 50 250 20], 'BackgroundColor', 'w',...
        'HorizontalAlignment', 'left', 'FontSize', 12,'FontUnits','pixels');
end



% cohort plot parameters
nc = numel(cTime);
chVec = nCh:-1:1;
hb = 0.1;
cmap = cell(1,nCh);
cv = cell(1,nCh);
for ch = 1:nCh
    v = mod(hues(ch)+linspace(-hb, hb, nc)', 1);
    cmap{ch} = hsv2rgb([v ones(nc,1) 0.9*ones(nc,1)]);
    cv{ch} = hsv2rgb([v 0.4*ones(nc,1) ones(nc,1)]);
end

% determine y-scale
ymax = zeros(1,nc);
for ci = nc:-1:1
    ymax(ci) = max(prctile(AMat(cidx==ci,:,1),90,1));
end
ymax = max(ymax);
fset = loadFigureSettings();

% include additional cohort at end
cohortLabels = arrayfun(@(i) [num2str(cohortBounds(i)) '-' num2str(cohortBounds(i+1)-data(1).framerate) 's'], 1:nc, 'Unif', 0);
% XTick = (cohortBounds(1:end-1)+[cohortBounds(2:end-1) cohortBounds(end)-data(1).framerate])/2;
XTick = (cohortBounds(1:end-1) + cohortBounds(2:end))/2;
XTick = [XTick 2*XTick(end)-XTick(end-1)];
% XTick = (cohortBounds(1:end-1)+[cohortBounds(2:end-1) cohortBounds(end)-framerate])/2;

updateCohorts();

set(ha([1 2]), 'XTick', 0:40:250);
        

    function updateCohorts()
        hold(ha(3), 'off');
        hold(ha(4), 'off');

%         if nCh>1
%             selIndex1 = mat2cell(maxA(:,1)>=at0(1) & maxA(:,2)>=at0(2) & lftV>=lt0, tracksPerCohort, 1);
%             selIndex2 = mat2cell(maxA(:,1)>=at0(1) & maxA(:,2)<at0(2)  & lftV>=lt0, tracksPerCohort, 1);
%         else
%             selIndex1 = mat2cell(maxA(:,1)>=at0(1) & lftV>=lt0, tracksPerCohort, 1);
%             selIndex2 = mat2cell(maxA(:,1)<at0(1) & lftV>=lt0, tracksPerCohort, 1);
%         end
        
        % Plot mean/median in front
        nTracksSel = zeros(nc,2);
        for ch = chVec
            for c = nc:-1:1
                if nCh>1
                    Apos = AMat(maxAMat(:,1)>=at0(1) & maxAMat(:,2)>=at0(2) & lftAMat>=lt0 & cidx==c,:,ch);
                    Aneg = AMat(maxAMat(:,1)>=at0(1) & maxAMat(:,2)< at0(2) & lftAMat>=lt0 & cidx==c,:,ch);
                else % above vs. below intensity threshold
                    Apos = AMat(maxAMat(:,1)>=at0(1) & lftAMat>=lt0 & cidx==c,:,1);
                    Aneg = AMat(maxAMat(:,1)<at0(1) & lftAMat>=lt0 & cidx==c,:,1);
                end
                plot(ha(3), cTime{c}, median(Apos(:,1:numel(cTime{c})),1), '-', 'Color', cmap{ch}(c,:), 'LineWidth', 1);
                plot(ha(4), cTime{c}, median(Aneg(:,1:numel(cTime{c})),1), '-', 'Color', cmap{ch}(c,:), 'LineWidth', 1);
                hold(ha(3), 'on');
                hold(ha(4), 'on');
                nTracksSel(c,:) = [size(Apos,1) size(Aneg,1)];
            end
        end

        %YTick = get(ha(1), 'YTick');
        %di = YTick(2)-YTick(1);
        %ya = -di:di:ceil(prctile(maxA(:,1),90)/di)*di;

        set(ha(3:4), 'XLim', [-15 125], 'YLim', [-0.1*ymax ymax]);

        % print # tracks/cohort
        text(-12, ymax, '# obj:', 'Parent',  ha(3), fset.ifont{:},...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        text(-12, ymax, '# obj:', 'Parent',  ha(4), fset.ifont{:},...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        for c = 1:nc
            text(XTick(c), ymax, num2str(nTracksSel(c,1)), 'Parent',  ha(3), fset.ifont{:},...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            text(XTick(c), ymax, num2str(nTracksSel(c,2)), 'Parent',  ha(4), fset.ifont{:},...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        end
        % indicate # objects > last cohort
        %text(XTick(end), ymax, num2str(sum(maxA(:,1)>=at0(1) & lftV>=lt0 & lftV>cohortBounds(end))), 'Parent',  ha(3), fset.ifont{:},...
        %       'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            
        set(ha([3 4]), 'XTick', XTick(1:end-1), 'XTickLabel', [], 'Box', 'off', 'TickDir', 'out');
        set(ha(4), 'XTickLabel', cohortLabels);
        rotateXTickLabels(ha(4), 'AdjustFigure', false);
        xlabel(ha(4), 'Lifetime cohort', fset.sfont{:});
        ylabel(ha(3), 'Median intensity (A.U.)', 'FontSize', 12);
        ylabel(ha(4), 'Median intensity (A.U.)', 'FontSize', 12);

        %set(ha(3), 'XColor', selColor2, 'YColor', selColor2);
        %set(ha(4), 'XColor', selColor1, 'YColor', selColor1);
        
        if nCh==1
            set(ha(2), 'XTick', 0:20:200, 'Box', 'off', 'TickDir', 'out');
        end
        
    end

    function click_Callback(varargin)
        switch gca
            case ha(1)
                a = get(gca, 'CurrentPoint');
                lt0 = round(a(1,1));
                at0(1) = round(a(1,2));
            case ha(2)
                a = get(gca,'CurrentPoint');
                lt0 = round(a(1,1));
                at0(2) = round(a(1,2));
        end
        
        refreshPlots();
        updateCohorts();
        set(gcf, 'WindowButtonMotionFcn', @drag, 'WindowButtonUpFcn', @stopDragging);
    end

    function drag(varargin)
        a = get(gca, 'CurrentPoint');
        switch gca
            case ha(1)
                lt0 = round(a(1,1));
                at0(1) = round(a(1,2));
            case ha(2)
                lt0 = round(a(1,1));
                at0(2) = round(a(1,2));
        end
        refreshPlots();
        updateCohorts();
    end

    function stopDragging(varargin)
        set(gcf, 'WindowButtonMotionFcn', '');
    end

    function refreshPlots()
        set(ht(:,1), 'XData', lt0*[1 1]);
        
        set(ht(1,2), 'YData', at0(1)*[1 1]);
        selCh1 = maxA(:,1)>=at0(1) & lftV>=lt0;

        set(hp1(1), 'XData', lftV(~selCh1), 'YData', maxA(~selCh1,1));
        set(hp1(2), 'XData', lftV(selCh1), 'YData', maxA(selCh1,1));
        set(htxt(1), 'String', [num2str(sum(selCh1)) '/' num2str(numel(selCh1))]);

        set(lftThresholdLabel, 'String', ['Lifetime threshold: ' num2str(lt0, '%.0f') ' s']);
        set(c1Label, 'String', ['Channel 1 threshold: ' num2str(at0(1), '%.0f')]);
        
        if nCh>1
            set(ht(2,2), 'YData', at0(2)*[1 1]);
            selCh2 = maxA(:,2)>=at0(2) & lftV>=lt0;
            
            set(hp2(1), 'XData', lftV(~selCh2 & selCh1), 'YData', maxA(~selCh2 & selCh1,2));
            set(hp2(2), 'XData', lftV(selCh2 & selCh1), 'YData', maxA(selCh2 & selCh1,2));
            set(hp2(3), 'XData', lftV(~selCh1), 'YData', maxA(~selCh1,2));
            set(htxt(2), 'String', [num2str(sum(selCh1&selCh2)) '/' num2str(sum(selCh1))]);
        
            set(c2Label, 'String', ['Channel 2 threshold: ' num2str(at0(2), '%.0f')]);
        else           
            histCCPs = hist(lftV(maxA(:,1)>=at0(1)), t).*w;
            histCCPs = histCCPs / sum(histCCPs) / framerate / bf;
            set(hp2(1), 'YData', histCCPs);
            
            histCSs = hist(lftV(maxA(:,1)<at0(1)), t).*w;
            histCSs = histCSs / sum(histCSs) / framerate / bf;
            set(hp2(2), 'YData', histCSs);
            
            nCCP = sum(maxA(:,1)>=at0(1));
            set(hl, 'String', {[' CCPs: ' num2str(nCCP/nTotal*100, '%.1f') '%'], [' CSs: ' num2str((nTotal-nCCP)/nTotal*100, '%.1f') '%']});
            %nTotal = numel(maxA(:,1));
            %hl = legend(ha(2), ['CCPs: ' num2str(nCCP/nTotal*100, '%.1f') '%'], ['CSs: ' num2str((nTotal-nCCP)/nTotal*100, '%.1f') '%']);
            %set(hl, 'EdgeColor', 'w');
        end
    end

end
