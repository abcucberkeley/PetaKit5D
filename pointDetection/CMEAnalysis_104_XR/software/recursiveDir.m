%[dirList] = recursiveDir(path, maxlevel) recursively lists directories found under 'path'
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

% Francois Aguet, 062813

function p = recursiveDir(d, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('d');
ip.addOptional('maxlevel', [], @(x) isempty(x) || (isnumeric(x) && abs(round(x))==x));
ip.addOptional('returnQueryPath', true, @islogical);
ip.parse(d, varargin{:});
maxlevel = ip.Results.maxlevel;

if ~strcmpi(d(end), filesep)
    d = [d filesep];
end

isvisibledir = @(s) [s.isdir]' & arrayfun(@(i) ~strcmp(i.name(1), '.'), s);

if ip.Results.returnQueryPath
    p = {d}; % return query path
else
    p = {};
end

files = p;
iter = 0;
while ~isempty(files) && (isempty(maxlevel) || iter<maxlevel)

    querylevel = files;
    
    files = cellfun(@dir, querylevel, 'unif', 0);
    files(cellfun(@isempty, files)) = [];
    
    % retain visible dirs in each set
    files = cellfun(@(i) {i(isvisibledir(i)).name}', files, 'unif', 0);
    
    % generate full paths
    for k = 1:numel(files)
        files{k} = cellfun(@(i) [querylevel{k} i filesep], files{k}, 'unif', 0);
    end
    files = vertcat(files{:});
    
    p = [p; files]; %#ok<AGROW>
    iter = iter + 1;
end
