function [v] = fixed_order_stat(fun, x, opt)
% sum, mean, var or std of all elements of x (x(:)) in an order that doesn't
% depend on the number of threads: MATLAB's own order at one thread.
%
%   v = fixed_order_stat('sum', x)             sum(x(:))
%   v = fixed_order_stat('sum', x, 'double')   sum(x(:), 'double')
%   v = fixed_order_stat('mean', x)            mean(x(:))
%   v = fixed_order_stat('mean', x, 'double')  mean(x(:), 'double')
%   v = fixed_order_stat('var', x, w)          var(x(:), w), w 0 (default) or 1
%   v = fixed_order_stat('std', x, w)          std(x(:), w)
%
% MATLAB (R2024b, on Linux and ARM macOS) sums fewer than 89000 elements on
% one thread at any thread count, so its own functions are used there, and for
% other arrays (integer, logical, complex, sparse, gpuArray) and options. From
% 89000 elements on it sums blocks of 4096 elements in four lanes and adds the
% block sums in index order, but splits the blocks among its threads, which
% changes the last bits with the number of threads (see decon_sum_double).
% These sums are made here in its order at one thread: of single data in
% single (the MATLAB code below), of other data and with 'double' in double
% (decon_sum_double). As in MATLAB, the mean is the sum divided by the number
% of elements in double (rounded to single for single data), the variance the
% sum of abs(x - mean) .^ 2 in the class of x divided by n - 1 (w = 0) or n
% (w = 1) in double (rounded to the class), and std its square root.


if nargin < 3
    opt = [];
end
if ~any(strcmp(fun, {'sum', 'mean', 'var', 'std'}))
    error('fixed_order_stat: unknown function %s, use ''sum'', ''mean'', ''var'' or ''std''.', fun);
end
if any(strcmp(fun, {'sum', 'mean'}))
    supported = isempty(opt) || (ischar(opt) && strcmp(opt, 'double'));
else
    supported = isempty(opt) || (isnumeric(opt) && isscalar(opt) && (opt == 0 || opt == 1));
end
if numel(x) < 89000 || ~(isa(x, 'single') || isa(x, 'double')) || ~isreal(x) || issparse(x) || ~supported
    if isempty(opt)
        v = feval(fun, x(:));
    else
        v = feval(fun, x(:), opt);
    end
    return;
end

n = numel(x);
inDouble = ~isempty(opt) && ischar(opt);
switch fun
    case 'sum'
        v = class_sum(x, inDouble);
    case 'mean'
        v = class_mean(x, inDouble);
    otherwise
        mu = class_mean(x, false);
        if isa(x, 'single')
            s = single_sum(x, mu);
        else
            s = decon_sum_double(abs(x(:) - mu) .^ 2);
        end
        if isempty(opt) || opt == 0
            v = cast(double(s) / (n - 1), 'like', x);
        else
            v = cast(double(s) / n, 'like', x);
        end
        if strcmp(fun, 'std')
            v = sqrt(v);
        end
end

end


function [s] = class_sum(x, inDouble)
% the sum in single for single data, else (or with inDouble) in double

if isa(x, 'single') && ~inDouble
    s = single_sum(x);
else
    s = decon_sum_double(x);
end

end


function [m] = class_mean(x, inDouble)
% the sum divided by the number of elements in double, rounded to the class
% of x (unless inDouble)

m = double(class_sum(x, inDouble)) / numel(x);
if ~inDouble
    m = cast(m, 'like', x);
end

end


function [s] = single_sum(x, mu)
% sum(x(:)) of single data in MATLAB's order at one thread, with mu the sum of
% abs(x(:) - mu) .^ 2: the blocks in pieces of 256, as decon_sum_double's
% MATLAB code (sum along the second dimension of the 4 x 1024 elements of a
% block adds each lane in index order, at any thread count; the last block is
% padded with -0, as v + -0 is v), in single

B = 4096;
nPiece = 256;
x = x(:);
n = numel(x);
nb = ceil(n / B);
bs = zeros(nb, 1, 'single');
for b = 1 : nPiece : nb
    e = min(b + nPiece - 1, nb);
    v = x((b - 1) * B + 1 : min(e * B, n));
    if nargin > 1
        v = abs(v - mu) .^ 2;
    end
    v(end + 1 : (e - b + 1) * B) = -0;
    L = sum(reshape(v, 4, B / 4, []), 2);
    bs(b : e) = ((L(1, :) + L(2, :)) + L(3, :)) + L(4, :);
end
s = bs(1);
for b = 2 : nb
    s = s + bs(b);
end

end
