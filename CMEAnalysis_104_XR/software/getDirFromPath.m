%[dirName] = getDirFromPath(dpath) returns the last directory contained in the input path
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

% Francois Aguet, November 2010

function [dirName, dirPath] = getDirFromPath(dpath)

idx = regexp(dpath, filesep);
if idx(end) == length(dpath)
    dirName = dpath(idx(end-1)+1:end-1);
    dirPath = dpath(1:idx(end-1));
else
    dirName = dpath(idx(end)+1:end);
    dirPath = dpath(1:idx(end));
end