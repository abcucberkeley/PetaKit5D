function [out] = imresize3_fast(im, outSize, method)
% imresize3(im, outSize, method) for 'linear' and 'cubic' with imresize3's
% weights and resize mex, but without its passes along dimensions that keep
% their size (copies) and without permuting the volume for the z pass (it
% resizes a [ny * nx, nz] reshape along its columns). Other methods use
% imresize3, and so does any error (e.g., changed MATLAB internals).

if nargin < 3
    method = 'cubic';
end

switch lower(method)
    case 'linear'
        kernel = @matlab.images.internal.resize.triangle;
        kernelWidth = 2;
    case 'cubic'
        kernel = @matlab.images.internal.resize.cubic;
        kernelWidth = 4;
    otherwise
        out = imresize3(im, outSize, method);
        return;
end

try
    sz = size(im, 1 : 3);
    outSize = ceil(outSize(:)');
    scale = outSize ./ sz;
    % imresize3 resizes along the dimensions in ascending order of scale
    [~, order] = sort(scale);
    out = im;
    for d = order
        if outSize(d) == sz(d)
            continue;
        end
        [weights, indices] = matlab.images.internal.resize.contributions(sz(d), outSize(d), ...
            scale(d), kernel, kernelWidth, true);
        if d == 3
            osz = size(out, 1 : 3);
            out = matlab.images.internal.resize.imresizemex(reshape(out, osz(1) * osz(2), osz(3)), ...
                weights', indices', 2);
            out = reshape(out, osz(1), osz(2), outSize(3));
        else
            out = matlab.images.internal.resize.imresizemex(out, weights', indices', d);
        end
    end
catch ME
    disp(ME);
    out = imresize3(im, outSize, method);
end

end
