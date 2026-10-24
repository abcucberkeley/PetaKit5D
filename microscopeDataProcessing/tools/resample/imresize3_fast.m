function [out] = imresize3_fast(im, outSize, method, outBbox, outClass)
% imresize3(im, outSize, method) for 'linear' and 'cubic' with imresize3's
% weights and resize mex, but without its passes along dimensions that keep
% their size (copies) and without permuting the volume for the z pass (it
% resizes a [ny * nx, nz] reshape along its columns). Other methods use
% imresize3, and so does any error (e.g., changed MATLAB internals).
%
% outBbox (optional, [ymin, xmin, zmin, ymax, xmax, zmax]): only compute this
% region of the output, from the input region it depends on (the same values
% as cropping the whole output). outClass (optional): class to resize in
% (default: the class of im); only the input region that is used is converted.

if nargin < 3
    method = 'cubic';
end
if nargin < 4
    outBbox = [];
end
if nargin < 5 || isempty(outClass)
    outClass = class(im);
end

switch lower(method)
    case 'linear'
        kernel = @matlab.images.internal.resize.triangle;
        kernelWidth = 2;
    case 'cubic'
        kernel = @matlab.images.internal.resize.cubic;
        kernelWidth = 4;
    otherwise
        out = resize_with_imresize3(im, outSize, method, outBbox, outClass);
        return;
end

try
    sz = size(im, 1 : 3);
    outSize = ceil(outSize(:)');
    scale = outSize ./ sz;
    % imresize3 resizes along the dimensions in ascending order of scale
    [~, order] = sort(scale);
    if isempty(outBbox)
        out = im;
        if ~isa(out, outClass)
            out = cast(out, outClass);
        end
        for d = order
            if outSize(d) == sz(d)
                continue;
            end
            [weights, indices] = matlab.images.internal.resize.contributions(sz(d), outSize(d), ...
                scale(d), kernel, kernelWidth, true);
            out = resize_along_dim(out, d, outSize(d), weights, indices);
        end
    else
        % weights of the output region and the input region they use
        weights = cell(1, 3);
        indices = cell(1, 3);
        inRange = cell(1, 3);
        for d = 1 : 3
            if outSize(d) == sz(d)
                inRange{d} = outBbox(d) : outBbox(d + 3);
                continue;
            end
            [w, ind] = matlab.images.internal.resize.contributions(sz(d), outSize(d), ...
                scale(d), kernel, kernelWidth, true);
            weights{d} = w(outBbox(d) : outBbox(d + 3), :);
            ind = ind(outBbox(d) : outBbox(d + 3), :);
            s = min(ind(:));
            inRange{d} = s : max(ind(:));
            indices{d} = ind - s + 1;
        end
        out = im(inRange{1}, inRange{2}, inRange{3});
        if ~isa(out, outClass)
            out = cast(out, outClass);
        end
        for d = order
            if outSize(d) == sz(d)
                continue;
            end
            out = resize_along_dim(out, d, outBbox(d + 3) - outBbox(d) + 1, weights{d}, indices{d});
        end
    end
catch ME
    disp(ME);
    out = resize_with_imresize3(im, outSize, method, outBbox, outClass);
end

end


function [out] = resize_along_dim(in, d, outLen, weights, indices)
% resize along dimension d with imresize3's resize mex

if d == 3
    isz = size(in, 1 : 3);
    out = matlab.images.internal.resize.imresizemex(reshape(in, isz(1) * isz(2), isz(3)), ...
        weights', indices', 2);
    out = reshape(out, isz(1), isz(2), outLen);
else
    out = matlab.images.internal.resize.imresizemex(in, weights', indices', d);
end

end


function [out] = resize_with_imresize3(im, outSize, method, outBbox, outClass)
% imresize3 of the whole image, cropped to outBbox

if ~isa(im, outClass)
    im = cast(im, outClass);
end
out = imresize3(im, outSize, method);
if ~isempty(outBbox)
    out = out(outBbox(1) : outBbox(4), outBbox(2) : outBbox(5), outBbox(3) : outBbox(6));
end

end
