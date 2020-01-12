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
    case {2,3}
        idxRGB = order;
    otherwise
        error('Max. 3 channels for RGB display.');
end