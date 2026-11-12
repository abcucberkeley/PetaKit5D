function n06_lsq_assign_ml(work, threads, repeats)
% N06 benchmark: stitch_global_grid_assignment on the 900-tile problems that
% n06_lsq_assign_setup.py wrote, per thread count
names = {'grid', 'free', 'tight'};
for k = 1 : numel(names)
    p = load(sprintf('%s/n06_p900_%s.mat', work, names{k}));
    args = {p.nF, p.mx, p.ash, p.ov, p.mas, p.thr, p.tileIdx, p.aw};
    for t = threads
        maxNumCompThreads(t);
        evalc('d = stitch_global_grid_assignment(args{:});');  % warm-up
        times = zeros(1, repeats);
        for r = 1 : repeats
            t0 = tic;
            evalc('d = stitch_global_grid_assignment(args{:});');
            times(r) = toc(t0);
        end
        pk5d_bench_record(sprintf('assign_900_%s', names{k}), t, times);
    end
end
end
