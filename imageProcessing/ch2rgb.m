% Generates an RGB image from grayscale input channels. Channels are normalized to [0..255].
% Inputs can be empty.

% Francois Aguet, June 2010.

function imRGB = ch2rgb(R, G, B)

imRGB = cat(3, R, G, B);
[ny, nx, ~] = size(imRGB);

if isempty(R)
    R = zeros(ny,nx);
end
if isempty(G)
    G = zeros(ny,nx);
end
if isempty(B)
    B = zeros(ny,nx);
end
imRGB = uint8(cat(3, scaleContrast(R), scaleContrast(G), scaleContrast(B)));