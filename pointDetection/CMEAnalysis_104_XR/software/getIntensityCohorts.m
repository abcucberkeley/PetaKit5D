
% Francois Aguet, 2013
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

function [res, ctime] = getIntensityCohorts(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('lftData', [], @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 120]);
ip.addParamValue('Rescale', true, @islogical);
ip.addParamValue('ExcludeVisitors', true, @islogical);
ip.addParamValue('Cutoff_f', 5);
ip.addParamValue('Alpha', 0.05);
ip.addParamValue('DetectionMode', [], @(x) isempty(x) || any(strcmpi(x, {'m', 's'})));
ip.addParamValue('SelIndex', [], @iscell); % selection index for each lftData
ip.addParamValue('LftDataName', 'lifetimeData.mat');
ip.parse(data, varargin{:});
cohortBounds = ip.Results.CohortBounds_s;

lftData = ip.Results.lftData;
if isempty(lftData)
    lftData = getLifetimeData(data,...
        'LifetimeData', ip.Results.LftDataName, 'Scale', ip.Results.Rescale,...
        'Cutoff_f', ip.Results.Cutoff_f, 'ReturnValidOnly', true,...
        'ExcludeVisitors', ip.Results.ExcludeVisitors);
end

framerate = data(1).framerate;

nCh = numel(data(1).channels);
nd = numel(data);
nc = numel(cohortBounds)-1;

kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1);

% selIndex = ip.Results.SelIndex;
% if isempty(selIndex)
%     selIndex = arrayfun(@(i) size(i.A,1), lftData, 'unif', 0);
% end

% # data points in cohort (including buffers)
b = size(lftData(1).sbA,2);
iLength = arrayfun(@(c) floor(mean(cohortBounds([c c+1]))/framerate) + 2*b, 1:nc);

% time vectors for cohorts
ctime = arrayfun(@(i) (-b:i-b-1)*framerate, iLength, 'unif', 0);

res(1:nd) = struct('aInterp', [], 'rInterp', []);
cohortBounds(end) = cohortBounds(end)+framerate;
for i = 1:nd
    
    aInterp = NaN(size(lftData(i).A));
    rInterp = NaN(size(lftData(i).A));
    
    nt = numel(lftData(i).lifetime_s);
    
    % cohort index
    % index of tracks < first cohort is 0; > last cohort is nc+1
    cidx = zeros(nt,1);
    for c = 1:nc
        cidx(cohortBounds(c)<=lftData(i).lifetime_s & lftData(i).lifetime_s<cohortBounds(c+1)) = c;
    end
    cidx(cohortBounds(c+1)<lftData(i).lifetime_s) = nc+1;
    
    % valid cohort index
    vidx = cidx>0 & cidx<=nc;
    
    for ch = 1:nCh % channels
        % loop through tracks and interpolate to cohort time vector
        for t = 1:nt
            if vidx(t)
                tidx = 1:lftData(i).trackLengths(t);
                A = [lftData(i).sbA(t,:,ch) lftData(i).A(t,tidx,ch) lftData(i).ebA(t,:,ch)];
                sigma_r = [lftData(i).sbSigma_r(t,:,ch) lftData(i).sigma_r(t,tidx,ch) lftData(i).ebSigma_r(t,:,ch)];
            
                % interpolate to mean length
                np = iLength(cidx(t));
                xi = linspace(1,lftData(i).trackLengths(t) + 2*b, np);
                aInterp(t,1:np,ch) = binterp(A, xi);
                rInterp(t,1:np,ch) = binterp(sigma_r, xi);
            end
        end
    end
    
    if nCh>1 && ~isempty(ip.Results.DetectionMode) % apply further selection based on slave channel
        if strcmpi(ip.Results.DetectionMode, 'm')
            idx = lftData(i).significantMaster(:,2)==1;
        elseif strcmpi(ip.Results.DetectionMode, 's')
            idx = lftData(i).significantSlave(:,2)==1;
        end
        vidx = vidx & idx;
    end
    res(i).aInterp = aInterp(vidx,:,:);
    res(i).rInterp = kLevel*rInterp(vidx,:,:);
    res(i).cidx = cidx(vidx,:,:);
    res(i).lifetime = lftData(i).lifetime_s(vidx);
    res(i).cidxAll = cidx;
end
    
%                 % split as a function of slave channel signal
%                 if isfield(lftData(i), 'significantMaster')
%                     sigIdx = lftData(i).significantMaster(:,ch)==1;
%                     res(i).sigIdx{c}(:,ch) = sigIdx(cidx);
%                 else
%                     res(i).sigIdx{c}(:,ch) = ones(numel(cidx),1);
%                 end
%                 if ch==1
%                     res(i).trackIdx{c} = lftData(i).index(cidx);
%                 end
%             else
%                 res(i).interpTracks{ch,c} = NaN(1,iLength(c));
%                 res(i).interpSigLevel{ch,c} = NaN(1,iLength(c));
%                 res(i).sigIdx{c}(:,ch) = NaN;
%                 if ch==1
%                     res(i).trackIdx{c} = NaN;
%                 end
%             end
%         end
%     end

