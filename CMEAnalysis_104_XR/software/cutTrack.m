%[tracks] = cutTrack(track, cutIdx) splits a track at the specified indexes
% The index positions are excluded from the resulting segments
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

function tracks = cutTrack(track, cutIdx)

np = numel(track.t);
ns = numel(cutIdx)+1;

cutGaps = cellfun(@(i) any(ismember(i, cutIdx)), track.gapIdx);

ub = [cellfun(@(i) i(1)-1, track.gapIdx(cutGaps)) np];
lb = [1 cellfun(@(i) i(end)+1, track.gapIdx(cutGaps))];

% fields with track-length content
fieldNames = fieldnames(track);
idx = structfun(@(i) size(i,2), track)==np;
fnames = fieldNames(idx);
bnames = {'x', 'y', 'A', 'c', 'A_pstd', 'sigma_r', 'SE_sigma_r', 'pval_Ar'};

framerate = track.t(2)-track.t(1);

tracks(1:ns) = struct(track); 
for s = 1:ns
    range = lb(s):ub(s);
    for f = 1:numel(fnames)
        tracks(s).(fnames{f}) = track.(fnames{f})(:,range);
    end
    tracks(s).start = tracks(s).f(1);
    tracks(s).end = tracks(s).f(end);
    tracks(s).lifetime_s = (tracks(s).end-tracks(s).start+1)*framerate;
    tracks(s).seqOfEvents = [tracks(s).start 1 1 NaN; tracks(s).end 2 1 NaN];

    if s==1
        tracks(s).startBuffer = track.startBuffer;
    else
        nb = numel(track.startBuffer.t);
        tracks(s).startBuffer.t = tracks(s).t(1) + (-nb:-1)*framerate;
        for f = 1:numel(bnames)
            % within track
            nbIn = min(cutIdx(s-1),nb);
            nbOut = nb-nbIn;
            tracks(s).startBuffer.(bnames{f})(:,nbOut+1:nb) = track.(bnames{f})(:,cutIdx(s-1)-(nbIn-1:-1:0));
            % within parent start buffer
            if nbOut>0
                tracks(s).startBuffer.(bnames{f})(:,1:nbOut) = track.startBuffer.(bnames{f})(:,end-nbOut+1:end);
            end
        end
    end
    
    if s==ns
        tracks(s).endBuffer = track.endBuffer;
    else
        nb = numel(track.endBuffer.t);
        tracks(s).endBuffer.t = tracks(s).t(end) + (1:nb)*framerate;
        for f = 1:numel(bnames)
            tracks(s).endBuffer.(bnames{f}) = track.(bnames{f})(:,cutIdx(s)+(0:nb-1));
        end
    end
    % gaps in current segment
    idx = cellfun(@(i) all(lb(s)<=i & i<=ub(s)), track.gapIdx);
    tracks(s).gapStatus = track.gapStatus(idx);
    tracks(s).gapIdx = track.gapIdx(idx);
end