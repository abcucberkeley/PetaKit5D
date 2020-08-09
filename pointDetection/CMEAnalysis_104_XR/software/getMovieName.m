%[str] = getMovieName(data) returns the identifier string ' date movieName' for each movie in data
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

% Francois Aguet 08/2013

function str = getMovieName(data)

nd = numel(data);
str = cell(1,nd);
for i = 1:nd
    if isempty(data(i).date)
        str{i} = [' ' getCellDir(data(i))];
    else
        str{i} = [' ' num2str(data(i).date) filesep getCellDir(data(i))];
    end
end
if nd==1
    str = str{1};
end