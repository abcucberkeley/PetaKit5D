function [s] = decon_sum_double(x)
% sum(x(:), 'double') of a real single or double array in an order that
% doesn't depend on the number of threads: MATLAB's own order at one thread.
%
% MATLAB (R2024b, R2025b) sums fewer than 89000 elements on one thread at any
% thread count (in index order below 44 elements, else in four interleaved
% lanes), so its sum is used there. From 89000 elements on it sums blocks of
% 4096 elements (each in four lanes, element t of a block in lane mod(t, 4),
% added as ((l1 + l2) + l3) + l4) and adds the block sums in index order, but
% splits the blocks among its threads, which changes the last bits of the sum
% with the number of threads. Those are summed here as on one thread, with
% decon_lucy_accel_mex in parallel or, without it, with the MATLAB code below
% (the same bits). Other arrays (gpuArray, integer, complex, sparse) are
% summed by MATLAB.


if numel(x) < 89000 || ~(isa(x, 'single') || isa(x, 'double')) || ~isreal(x) || issparse(x)
    s = sum(x(:), 'double');
    return;
end

if exist('decon_lucy_accel_mex', 'file') == 3
    try
        s = decon_lucy_accel_mex('sum', x, maxNumCompThreads);
        return;
    catch ME
        mex_fallback('decon_lucy_accel_mex', ME);
    end
else
    mex_fallback('decon_lucy_accel_mex');
end
s = block_sum(x);

end


function [s] = block_sum(x)
% the blocks in pieces of 256: sum along the second dimension of the 4 x 1024
% elements of a block adds each lane in index order (as cumsum, at any thread
% count); the last block is padded with -0 (v + -0 is v). The block sums are
% added in a loop, which starts from the first one as the mex (from -0).

B = 4096;
nPiece = 256;
x = x(:);
n = numel(x);
nb = ceil(n / B);
bs = zeros(nb, 1);
for b = 1 : nPiece : nb
    e = min(b + nPiece - 1, nb);
    v = x((b - 1) * B + 1 : min(e * B, n));
    v(end + 1 : (e - b + 1) * B) = -0;
    L = sum(reshape(v, 4, B / 4, []), 2, 'double');
    bs(b : e) = ((L(1, :) + L(2, :)) + L(3, :)) + L(4, :);
end
s = bs(1);
for b = 2 : nb
    s = s + bs(b);
end

end
