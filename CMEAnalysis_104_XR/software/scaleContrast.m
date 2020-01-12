%out = scaleContrast(in, rangeIn, rangeOut) adjusts the contrast of the input
%
% Inputs:
%         in : input signal
%    rangeIn : input range. If empty, [min(in(:)) max(in(:))]
%   rangeOut : output range
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

% Francois Aguet (Last modified: 03/22/2011)

function out = scaleContrast(in, rangeIn, rangeOut)

if nargin<2 || isempty(rangeIn)
    rangeIn = [min(in(:)) max(in(:))];
end
if nargin<3 || isempty(rangeOut)
    rangeOut = [0 255];
end

if rangeIn(2)-rangeIn(1) ~= 0
    out = (double(in)-rangeIn(1)) / diff(rangeIn) * diff(rangeOut) + rangeOut(1);
else
    out = zeros(size(in));
end