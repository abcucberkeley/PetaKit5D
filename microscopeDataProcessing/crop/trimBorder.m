function [im] = trimBorder(im, borderSize, method)
% trim the border of an image
%   method: pre, post, both

if nargin < 3
    method = 'both';
end

if all(borderSize == 0)
    return;
end

sz = size(im);

switch method
    case 'pre'
        s = borderSize + 1;
        t = sz;
    case 'post'
        s = ones(1, 3);
        t = sz - borderSize;
    case 'both'
        s = borderSize + 1;
        t = sz - borderSize;
end

im = im(s(1) : t(1), s(2) : t(2), s(3) : t(3));

end