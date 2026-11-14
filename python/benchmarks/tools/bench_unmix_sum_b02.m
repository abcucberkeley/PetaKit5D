function bench_unmix_sum_b02(work, t)
% The unmixing kernel benchmark (B02), MATLAB side (Python: bench_unmix_sum_b02.py): this
% process pins itself to the CPUs of <work>/cpus_<t>.txt before any mex runs; then the sum of
% one demo tile pair in memory (unmix_channels_sum, the mex path: two uint16 channels, [1 -0.1],
% backgrounds [100 100], offset 1, uint16), 3 timed calls after one untimed (each output
% cleared after its toc).

cpus = strtrim(fileread(fullfile(work, sprintf('cpus_%d.txt', t))));
[s, msg] = system(sprintf('taskset -a -p -c %s %d', strrep(cpus, ' ', ','), feature('getpid')));
if s ~= 0
    error('taskset failed: %s', msg);
end
maxNumCompThreads(t);
tile = 'Scan_Iter_0000_0000_%s_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t';
ims = cellfun(@(c) readtiff([work, '/demo/', sprintf(tile, c), '.tif']), {'CamA', 'CamB'}, 'unif', 0);
pk5d_bench_mark('loaded');
times = zeros(1, 3);
for k = 0 : 3
    t0 = tic;
    out = unmix_channels_sum(ims, [1, -0.1], [100, 100], 1, 'uint16');
    dt = toc(t0);
    clear out;
    if k > 0
        times(k) = dt;
    end
end
pk5d_bench_record('sum', t, times);
end
