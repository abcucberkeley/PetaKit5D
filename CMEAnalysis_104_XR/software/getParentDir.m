% [cpath] = getParentDir(cpath) Returns the parent directory's path from the input path
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

% Francois Aguet, 11/05/2010

function cpath = getParentDir(cpath, level)
if nargin<2
    level = 1;
end

fsIdx = regexp(cpath, filesep);

if strcmp(cpath(end), filesep)
    cpath = cpath(1:fsIdx(end-level));
else
    cpath = cpath(1:fsIdx(end-level+1));
end