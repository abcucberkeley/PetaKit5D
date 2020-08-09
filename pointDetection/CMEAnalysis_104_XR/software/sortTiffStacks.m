%sortTiffStacks() 
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


% 1) Input directory contains TIFF stacks
%    Stacks are sorted into individual directories with the same file name.
% 2) List of TIFF stacks, including multiple channels. File names must contain a movie identifier.
%    Channels are sorted into their own directories. 
%    Example: cell1_1s_488.tif cell1_1s_561.tif cell2_1s_488.tif cell2_1s_561.tif,
% 3) 


% Francois Aguet, 09/27/2013

function sortTiffStacks(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('spath', [], @(x) ischar(x) || isempty(x) || iscell(x));
ip.addParamValue('ChannelNames', [], @iscell);
ip.addParamValue('MovieSelector', 'cell', @ischar);
ip.addParamValue('UseFileName', true, @islogical);
ip.parse(varargin{:});
stkpath = ip.Results.spath;
chSpec = ip.Results.ChannelNames;

if isempty(stkpath)
   stkpath = uigetdir('Select directory containing the STK files:'); 
   if (stkpath == 0)
       return;
   end
end

% Recursive call if input is cell array
if iscell(stkpath)
    cellfun(@(x)sortTiffStacks(x, chSpec), stkpath);
    return
end

stkpath = [stkpath filesep];
stkList = vdir(stkpath);

% if list of dirs: assume multi-channel, sort into sub-directories
if all([stkList.isdir])
    dirList = {stkList.name};
    for i = 1:numel(dirList)
        stkList = vdir([stkpath dirList{i}]);
        stkList = {stkList.name};
        stkList = stkList(~cellfun(@isempty, regexpi(stkList, '(\.tiff?|\.stk)$')));
        if ~isempty(chSpec)
            for c = 1:numel(stkList); % loop through channels
                [~,~,ext] = fileparts(stkList{c});
                destDir = [stkpath dirList{i} filesep chSpec{c} filesep];
                [~,~] = mkdir(destDir);
                [~] = movefile([stkpath dirList{i} filesep '*' chSpec{c} '*' ext], destDir);
            end
        else
            for c = 1:numel(stkList);
                [~,dname] = fileparts(stkList{c});
                destDir = [stkpath dirList{i} filesep dname filesep];
                [~,~] = mkdir(destDir);
                [~] = movefile([stkpath dirList{i} filesep stkList{c}], [destDir stkList{c}]);
            end
        end
        
    end
else % list of stacks
    stkList = {stkList(~[stkList.isdir]).name};
    idx = ~cellfun(@isempty, regexpi(stkList, '(\.tiff?|\.stk)$'));
    stkList = stkList(idx);
    
    N = length(stkList);
    if N==0
        fprintf('No TIFF files found in input directory.\n');
        return
    end
    
    nCh = numel(chSpec);
    
    %N = numel(movieID);
    for k = N:-1:1
        if ip.Results.UseFileName
            [~,dirName, ext] = fileparts(stkList{k});
        else
%             % stack index/identifier following selector (i.e., 12 for 'Cell12_2s')
%             movieID = str2double(regexpi(stkList{k}, ['(?<=' ip.Results.MovieSelector ')\d+'], 'match', 'once'));
%             % frame rate (if in file name)
%             framerate = str2double(regexpi(stkList{k}, '(?<=_)\d+(?=s)', 'match', 'once'));
%             [movieID,idx] = unique(movieID);
%             framerate = framerate(idx);
%             
%             if ~isnan(framerate(k))
%                 dirName = [ip.Results.MovieSelector num2str(movieID(k)) '_' num2str(framerate(k)) 's'];
%             else
%                 dirName = [ip.Results.MovieSelector num2str(movieID(k))];
%             end
        end
        
        [~,~] = mkdir([stkpath dirName]);
        
        if nCh>0
            for c = 1:nCh
                % create sub-directory for each channel
                [~,~] = mkdir([stkpath dirName filesep chSpec{c}]);
                % move matching files to sub-directory
                movefile([stkpath ip.Results.MovieSelector num2str(movieID(k)) '*' chSpec{c} '*.tif*'],...
                    [stkpath dirName filesep chSpec{c}]);
            end
        else
            movefile([stkpath dirName ext], [stkpath dirName]);
        end
    end
end


function d = vdir(path)
d = dir(path);
d = d(cellfun(@(i) ~strcmpi(i(1), '.'), {d.name}));
