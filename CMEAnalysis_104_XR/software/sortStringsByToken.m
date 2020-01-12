%[s, rm] = sortStringsByToken(s, token, mode) sorts input array by a number sequence preceding or following a specified token
%
% Inputs:
%      s : cell array of strings    
%  token : token to match
%   mode : 'pre' (sort by number sequence preceding token) or
%          'post' (by sequence following token)
%
% Outputs:
%      s : sorted array
%     rm : index of input that did not match token and was excluded from output
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

% Francois Aguet, 02/08/2013

function [s, rm] = sortStringsByToken(s, token, mode)

switch mode
    case 'post'
        q = ['(?<=' token ')\d+'];
    case 'pre'
        q = ['\d+(?=' token ')'];
end
idx = str2double(regexpi(s, q, 'match', 'once'));
rm = find(isnan(idx));
s(rm) = [];
idx(rm) = [];
[~,idx] = sort(idx);
if ~isempty(idx)
    s = s(idx);
end