function [dname, dpath] = getCellDir(data)

chParents = cellfun(@(c) getParentDir(c), data.channels, 'unif', 0);
sParent = getParentDir(data.source);

v = strcmp(sParent, chParents);

if numel(data.channels)==1 || ~all(v)
    dpath = data.source;
else % all channels at same level, below 'source'
    dpath = sParent;
end
dname = getDirFromPath(dpath);
dpath = getParentDir(dpath);
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
