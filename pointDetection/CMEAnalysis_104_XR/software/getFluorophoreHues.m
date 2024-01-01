% colors = getHuesFromMarkers(markers) assigns primary hues to the markers
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

% Francois Aguet, February 2011

function colors = getFluorophoreHues(markers)

if ischar(markers)
    markers = {markers};
end

hue = arrayfun(@(x) rgb2hsv(wavelength2rgb(x)), name2wavelength(markers), 'UniformOutput', false);
hue = vertcat(hue{:});
hue = hue(:,1);

[hue, sortIdx] = sort(hue, 'descend');

N = length(markers);

switch N
    case 1
        hueV = [120 0 240]/360;
        D = (hue-hueV).^2;
        colors = hueV(D==min(D));
    case 2
        hueV = [120 0]/360; % green/red only
        D = (hue(1)-hueV).^2;
        [~,idx] = sort(D);
        colors = arrayfun(@(x) hueV(x), idx);
    case {3,4}
        colors = NaN(1,N);
        hueIdx = 1:N;
        % first assign green and red, then remaining hues
        D = (hue-1/3).^2;
        i = find(D==min(D));
        colors(hueIdx(i)) = 1/3; % green
        hueIdx(i) = [];
        hue(i) = [];
        
        D = hue.^2;
        i = find(D==min(D));
        colors(hueIdx(i)) = 0; % red
        hueIdx(i) = [];
        hue(i) = [];
                
        hueV = [240 180 60 30]/360; % blue cyan yellow orange
        D = (hue(1)-hueV).^2;
        colors(hueIdx(1)) = hueV(D==min(D));
        if N==4
            D = (hue(2)-hueV).^2;
            colors(hueIdx(2)) = hueV(D==min(D));
        end
    otherwise
        error('Unsupported number of colors.');
end        

colors = colors(sortIdx);