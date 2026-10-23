function [out] = decon_lucy_rfft(useGPU, varargin)
% run a step of the half-spectrum FFT decon ('otf', 'omw' or 'rlstep', see
% decon/mex/decon_lucy_rfft_mex.cpp) with decon_lucy_rfft_mex on the CPU,
% with the threads of maxNumCompThreads and the same FFTs as MATLAB's, or
% decon_lucy_rfft_gpu_mex on the GPU


persistent grain

if useGPU
    out = decon_lucy_rfft_gpu_mex(varargin{:});
    return;
end

if isempty(grain)
    grain = matlab_fft_grain();
end
out = decon_lucy_rfft_mex(varargin{:}, maxNumCompThreads, grain);

end


function [grain] = matlab_fft_grain()
% MATLAB's FFTs transform the lines of each dimension in groups of four and
% some lines alone, which for some lengths changes the last bits. When MATLAB
% runs its FFTs in parallel (it decides that at its first FFT, from the number
% of threads then), it splits the lines into pieces of at most 32 lines and
% runs the last (length % 4) lines of each piece alone, otherwise only the last
% lines of a dimension. Find out which split this MATLAB session uses from the
% FFT of test data, so that the mex runs the same lines alone (grain 32: the
% split, 0: none).

grains = [32, 0];
for n = [64, 128, 512]
    m = 2000;
    X = single(mod(reshape(0 : n * m - 1, n, m) * 0.6180339887, 1));
    Y = fft(X, [], 1);
    Y = Y(1 : n / 2 + 1, :);
    match = false(size(grains));
    for i = 1 : numel(grains)
        match(i) = isequal(decon_lucy_rfft_mex('fft1', X, maxNumCompThreads, grains(i)), Y);
    end
    if sum(match) == 1
        grain = grains(match);
        return;
    end
end
% the test data can't tell (or neither split matches): the parallel split
grain = grains(1);

end
