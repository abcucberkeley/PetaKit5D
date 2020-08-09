%[fi] = interpln(x, f, xi) implements lower-neighbor interpolation
% 
% Inputs:
%   x  : sample positions
%   f  : sample values
%   xi : interpolation positions
%
% Outputs:
%   fi : interpolated values
%
% Example:
%  interpln([0 0.5 2.7 3.2], [0 1 0 1], [1 3])
%  returns [1 0]

% Francois Aguet, 03/24/2014

function [fi] = interpln(x, f, xi)

% find lower bound index
% 1) subtract x from xi
idx = repmat(xi(:), [1 numel(x)]) - repmat(x(:)', [numel(xi) 1]);
% 2) position of last >=0 value is index
[~,idx] = max(idx(:,end:-1:1)>=0,[],2);
idx = numel(x)-idx+1;

fi = f(idx);