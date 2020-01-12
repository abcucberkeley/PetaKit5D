%formatTickLabels(h, varargin) formats the x- and y-axis tick labels using format specifiers
% The default format uses the largest number of decimals found, except for '0'.
%
% Inputs:
%         h : figure handle
%
% Optional:
%   xformat : string specifier for the x-axis format
%   yformat : string specifier for the y-axis format
%
% Examples: formatTickLabels(h);
%           formatTickLabels(h, [], '%.3f');
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

% Francois Aguet, 2012

function formatTickLabels(varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('h', gca, @(x) all(arrayfun(@ishandle, x)));
ip.addOptional('XFormat', [], @(x) isempty(x) || ischar(x));
ip.addOptional('YFormat', [], @(x) isempty(x) || ischar(x));
ip.addParamValue('FormatXY', [true true], @(x) islogical(x) && numel(x)==2);
ip.addParamValue('MaxDigits', 2, @isscalar);
ip.parse(varargin{:});
h = ip.Results.h;

for i = 1:numel(h)
    xfmt = ip.Results.XFormat;
    yfmt = ip.Results.YFormat;
    
    if ip.Results.FormatXY(1)
        xticks = cellstr(get(h(i), 'XTickLabel'));
        ppos = regexpi(xticks, '\.');
        nchar = cellfun(@numel, xticks);
        idx = ~cellfun(@isempty, ppos);
        if ~all(idx==0) || ~isempty(xfmt)
            val = cellfun(@str2num, xticks);
            if isempty(xfmt)
                nd = max(nchar(idx)-[ppos{idx}]');
                xfmt = ['%.' num2str(nd) 'f'];
            end
            xticks = arrayfun(@(i) num2str(i, xfmt), val, 'unif', 0);
            idx = find(val==0);
            if ~isempty(idx)
                xticks{idx} = '0';
            end
            set(h(i), 'XTick', val, 'XTickLabel', xticks);
        end
    end
    if ip.Results.FormatXY(2)
        yticks = cellstr(get(h(i), 'YTickLabel'));
        nchar = cellfun(@numel, yticks);
        
        ppos = regexpi(yticks, '\.');
        idx = cellfun(@isempty, ppos);
        tmp = num2cell(nchar(idx)+1);
        [ppos{idx}] = tmp{:};
        ppos = [ppos{:}]';
        
        if any(ppos~=0) && ~all(idx)
            val = cellfun(@str2num, yticks);
            if isempty(yfmt)
                nd = max(max(nchar-ppos),0);
                if nd>ip.Results.MaxDigits % remove labels with more than 2 digits
                    nd = ip.Results.MaxDigits;
                    clearIdx = nchar-ppos>ip.Results.MaxDigits;
                else
                    clearIdx = [];
                end
                yfmt = ['%.' num2str(nd) 'f'];
            end
            yticks = arrayfun(@(i) num2str(i, yfmt), val, 'unif', 0);
            if ~isempty(clearIdx)
                [yticks{clearIdx}] = deal('');
            end
            idx = find(val==0);
            if ~isempty(idx)
                yticks{idx} = '0';
            end
            set(h(i), 'YTick', val, 'YTickLabel', yticks);
        end
    end
end
