%[lengths, value] = binarySegmentLengths(signal) returns the lengths of the successive sequences of 0s and 1s in a binary signal
% Example: for the signal [0 1 1 0 0 0 1], lengths = [1 2 3 1], value = [0 1 0 1]
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

% Francois Aguet, 08/08/2012

function [lengths, value] = binarySegmentLengths(signal)

d = diff(signal);

changeIdx = find([0 d]~=0);

lengths = diff([1 changeIdx numel(signal)+1]);
value = signal([1 changeIdx]);