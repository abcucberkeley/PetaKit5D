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