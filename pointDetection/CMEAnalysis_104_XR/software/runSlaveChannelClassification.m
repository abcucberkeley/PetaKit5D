%[tracks] = runSlaveChannelClassification(data, tracks, varargin) identifies trajectories with significant slave channel fluorescence
%
% Input:
%          data : structure returned by loadConditionData()
%
% Options ('specifier', value):
%    'Cutoff_f' : minimum track length to consider for classification (in frames)
%
% Notes: This function modifies the output of runTrackProcessing(), 
%        by default saved in Tracking/ProcessedTracks.mat
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

% Francois Aguet, October 2010 (last modified: 10/09/2012)

function tracks = runSlaveChannelClassification(data, tracks, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('tracks', @isstruct);
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.addParamValue('AlphaSensitivity', 0.05, @isscalar);
ip.addParamValue('AlphaSigMaster', 0.05, @isscalar);
ip.addParamValue('amplitudeRatio', 0, @isscalar);
ip.parse(data, tracks, varargin{:});


% load cell mask
cellmask = logical(getCellMask(data));

% Determine master/slave channels
nc = length(data.channels); % number of channels
mCh = strcmp(data.source, data.channels);
sCh = setdiff(1:nc, mCh);

load([data.source 'Detection' filesep 'detection_v2.mat']);
sigma = frameInfo(1).s;
w = max(ceil(4*sigma(sCh)));


% calculate background distributions over 16 frames sampled across the movie
frameIdx = round(linspace(1, data.movieLength, 16));
nf = numel(frameIdx);

%===============================================================================
% Calculate background probabilities
%===============================================================================
bgA = cell(1,nf);
pSlaveSignal = NaN(nc,nf);
parfor f = 1:nf;
    k = frameIdx(f);
    
    %---------------------------------------------------------------------------
    % Masks: cell mask, CCP mask
    %---------------------------------------------------------------------------
    if iscell(data.framePaths{1}) %#ok<PFBNS>
        ccpMask = double(imread(data.maskPaths{k}));
    else
        ccpMask = double(readtiff(data.maskPaths, k));
    end
    ccpMask(ccpMask~=0) = 1;
    ccpMask = imdilate(ccpMask, strel('disk', w));
    cellBackgroundMask = (cellmask - ccpMask)==1;

    %---------------------------------------------------------------------------
    % Calculate probability of significant slave signal, background distribution
    %---------------------------------------------------------------------------
    bgA{f} = NaN(nc,sum(cellBackgroundMask(:)));
    pDetection = NaN(nc,1);
    for c = sCh

        if iscell(data.framePaths{1})
            frame = double(imread(data.framePaths{c}{k}));
        else
            frame = double(readtiff(data.framePaths{c}, k));
        end
        
        [A_est, ~, pval_Ar] = filterGaussianFit2D(frame, sigma(c)); %#ok<PFBNS>
        
        % background distribution (excluding CCS pixels)
        bgA{f}(c,:) = A_est(cellBackgroundMask(:));
        
        % significant slave signal probability (whole cell)
        pDetection(c) = sum(pval_Ar < ip.Results.Alpha) / sum(cellmask); %#ok<PFBNS>
    end
    pSlaveSignal(:,f) = pDetection;
end
pDetection = mean(pSlaveSignal,2);

for c = sCh
    fprintf('P(random detection in ch. %d) = %.3f\n', c, pDetection(c));
end

%===============================================================================
% Classify tracks in slave channels
%===============================================================================

bg95 = prctile([bgA{:}], 95, 2);

% add new fields
[tracks.significantMaster] = deal([]);
[tracks.significantVsBackground] = deal([]);
[tracks.significantSlave] = deal([]);

% Loops through all the tracks
nt = numel(tracks);
for k = 1:nt
    
    L = numel(tracks(k).t); % track length
    tracks(k).significantMaster = NaN(nc,1);
    tracks(k).significantVsBackground = NaN(nc,L);
    tracks(k).significantSlave = NaN(nc,1);

    for c = sCh % loop through all slave channels
        
        % calculate # of pixels used for detection
        npx = round((tracks(k).sigma_r(c,:) ./ tracks(k).SE_sigma_r(c,:)).^2/2+1);
        A = tracks(k).A(c,:);
        sigma_A = tracks(k).A_pstd(c,:);
        %ts.tracks(k).isDetected(c,:) = pval < 0.05; % == hval_Ar(c,:)
        

        % 1) Test whether the number of significant detections in the slave channel
        % of this track (independent of the master channel) is above chance. The expected
        % number of chance detections for a track of length L is given by the 
        % inverse binomial CDF:
        tracks(k).significantMaster(c) = nansum(tracks(k).hval_Ar(c,:)) > binoinv(0.95, L, pDetection(c));

        % 1bis) Add an optional step to strenghten the condition on the
        % signal
        %tracks(k).significantMaster(c) = nansum(tracks(k).pval_Ar(c,:)<ip.Results.AlphaSigMaster) > binoinv(0.95, L, pDetection(c));

        % 1ter) Add a condition on the amplitude ratio
        ampRatio=max(tracks(k).A(c,:))/max(tracks(k).A(1,:));
        tracks(k).significantMaster(c)=tracks(k).significantMaster(c) & (ampRatio > ip.Results.amplitudeRatio);        
        
        % 2) Test whether the number time points in the slave channel with a signal
        % above the 95th percentile of the background distribution is significant
        
        % 95th percentile and uncertainty
        sigma_r = bg95(c);
        SE_sigma_r = sigma_r ./ sqrt(2*npx-1);
        
        % t-test amplitude against 95th pct of background
        df2 = (npx-1) .* (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
        scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)./npx);
        T = (A - sigma_r) ./ scomb;
        pval = tcdf(-T, df2);
        hval = pval < ip.Results.AlphaSensitivity;
        tracks(k).significantVsBackground(c,:) = hval;
        
        % The background cutoff was set at the 95th percentile, thus there is a
        % 5% chance of a false positive, and the number of points must be greater
        % than the expected number of false positives for a track of length L:
        tracks(k).significantSlave(c) = nansum(hval) > binoinv(0.95, L, 0.05);
    end
end
