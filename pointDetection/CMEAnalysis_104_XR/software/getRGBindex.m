function idxRGB = getRGBindex(markers)

hue = cellfun(@(x) rgb2hsv(wavelength2rgb(name2wavelength(x))), markers, 'UniformOutput', false);
hue = vertcat(hue{:});
hue = hue(:,1); % retain only 'h' from hsv

N = length(hue);

[hue, sortIdx] = sort(hue, 'ascend');

order = zeros(1,N);
order(sortIdx) = 1:N;


hueRef = [0 1/3 2/3]; % RGB
N = length(markers);
switch N
    case 1
        J = (hue-hueRef).^2;
        idxRGB = find(J==min(J));
    %case 2 needs to be implemented
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
    case {2,3}
        idxRGB = order;
    otherwise
        error('Max. 3 channels for RGB display.');
end