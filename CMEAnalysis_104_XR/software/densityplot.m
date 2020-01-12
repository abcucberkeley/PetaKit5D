%[ha, dRange] = densityplot(x, y, varargin) generates a 2D histogram / density plot
% from the input vectors
%
% Inputs:
%          x : sample vector for the 1st dimension
%          y : sample vector for the 2nd dimension
%
% Optional inputs:
%         xv : bin centers for the 1st dimension
%         yv : bin centers for the 2nd dimension
%
% Parameters ('specifier', value)
%           'Parent' : handle of the parent axes
%  'DisplayFunction' : intensity mapping, i.e., @sqrt
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

% Francois Aguet, 02/19/2013

function [ha, dRange] = densityplot(x, y, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('xv', linspace(min(x), max(x), 100));
ip.addOptional('yv', linspace(min(y), max(y), 100));
ip.addParamValue('Parent', gca, @ishandle);
ip.addParamValue('Div', 1, @isscalar);
ip.addParamValue('NormX', false, @islogical);
ip.addParamValue('DisplayFunction', @(x) x, @(x) isa(x, 'function_handle'));
ip.parse(varargin{:});
xv = ip.Results.xv(:)';
yv = ip.Results.yv(:)';
ha = ip.Results.Parent;

M = hist3([y(:) x(:)], {yv, xv});
if ip.Results.NormX
    M = M./repmat(sum(M,1), [size(M,1) 1]);
end

M = ip.Results.DisplayFunction(M/ip.Results.Div);

imagesc(xv, yv, M, 'Parent', ha);
axis(ha, 'xy', 'tight');
dRange = [min(M(:)) max(M(:))];
