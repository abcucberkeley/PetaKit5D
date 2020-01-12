%[fpath] = getExpDir(data) returns the parent directory of the movies in 'data'
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

function [fpath] = getExpDir(data)

nd = numel(data);
fpath = cell(1,nd);
for k = 1:nd
    [~,fpath{k}] = getCellDir(data(k));
end
fpath = unique(fpath);
if numel(fpath)>1
    fprintf(2, 'Data sets from different conditions were combined.\n');
end
fpath = fpath{1};