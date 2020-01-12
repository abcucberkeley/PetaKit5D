%spath = getShortPath(data) returns the truncated path of a cell directory (3 levels) 
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

% Francois Aguet, 05/13/2011

function spath = getShortPath(data, level)
if nargin<2
    level = 3;
end

if numel(data.channels)>1
    mCh = strcmp(data.channels, data.source);
    sCh = setdiff(1:length(data.channels),mCh);
    spath = data.channels{sCh(1)};
    idx = regexp(spath, filesep);
    spath = spath(idx(end-level-1)+1:idx(end-1));
else
    spath = data.channels{1};
    idx = regexp(spath, filesep);
    spath = spath(idx(max(1,end-level))+1:idx(end)-1);
end
