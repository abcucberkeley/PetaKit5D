function bench_registration_s04(work, which, threads, repeats)
% the stitching registration on the inputs of registration_s04_setup.py, at each thread count;
% one warm-up run each. which: 'pair300' (cross_correlation_registration_3d of the two tiles,
% MATLAB's zarr) or 'demo' (the demo's 12 registration tasks, cross_correlation_registration_wrapper,
% each into a new folder)

if strcmp(which, 'pair300')
    s = jsondecode(fileread(fullfile(work, 'pair300.json')));
    f1 = fullfile(work, 'matlab', 't1.zarr');
    f2 = fullfile(work, 'matlab', 't2.zarr');
    for t = threads
        maxNumCompThreads(t);
        [rel, mx] = cross_correlation_registration_3d(f1, f2, '', s.c1', s.c2', s.co', s.vs', ...
            s.dom', 'downSample', s.ds', 'MaxOffset', s.mo');
        times = zeros(1, repeats);
        for r = 1 : repeats
            t0 = tic;
            [rel, mx] = cross_correlation_registration_3d(f1, f2, '', s.c1', s.c2', s.co', s.vs', ...
                s.dom', 'downSample', s.ds', 'MaxOffset', s.mo');
            times(r) = toc(t0);
        end
        pk5d_bench_record('pair300', t, times, 'shift', mat2str(rel), 'peak', double(mx));
    end
    return;
end
calls = jsondecode(fileread(fullfile(work, 'demo_matlab.json')));
for t = threads
    maxNumCompThreads(t);
    times = zeros(1, repeats + 1);
    for r = 1 : repeats + 1
        d = tempname(work);
        mkdir(d);
        cd(d);
        t0 = tic;
        for k = 1 : numel(calls)
            if iscell(calls)
                c = calls{k};
            else
                c = calls(k);
            end
            img2 = c.img2;
            if ischar(img2)
                img2 = {img2};
            end
            cross_correlation_registration_wrapper(c.img1, img2(:)', c.out, c.pairs, c.c1(:)', ...
                c.c2, c.co, c.vs(:)', c.dom(:)', 'Stitch2D', false, 'downSample', c.ds(:)', ...
                'MaxOffset', c.mo(:)', 'largeFile', false, 'mipDirStr', '', 'poolSize', [], ...
                'parseCluster', false);
        end
        times(r) = toc(t0);
        cd(work);
    end
    pk5d_bench_record('demo16', t, times(2 : end), 'pairs', numel(calls));
end

end
