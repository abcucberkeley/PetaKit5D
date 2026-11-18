function bench_xcorr_k06(work, threads, repeats)
% normxcorr3_max_shift (its mex, normxcorr3_window_mex) on the pairs of xcorr_k06_setup.py, at each
% thread count; one warm-up call per pair and thread count

pairs = {'ovl300', [0, 0, 0; 40, 40, 20]; 'demo2p', [-98, -98, 0; 302, 302, 20]};
for p = 1 : size(pairs, 1)
    name = pairs{p, 1};
    T = readtiff(fullfile(work, [name, '_T.tif']));
    A = readtiff(fullfile(work, [name, '_A.tif']));
    for t = threads
        maxNumCompThreads(t);
        [off, mx] = normxcorr3_max_shift(T, A, pairs{p, 2});
        times = zeros(1, repeats);
        for r = 1 : repeats
            tic;
            [off, mx] = normxcorr3_max_shift(T, A, pairs{p, 2});
            times(r) = toc;
        end
        pk5d_bench_record(['max_shift_', name], t, times, 'offset', mat2str(off), 'peak', double(mx));
    end
end

end
