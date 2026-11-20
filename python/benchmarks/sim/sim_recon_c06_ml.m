function sim_recon_c06_ml(work, t)
% The SIM benchmark (C06), MATLAB side (Python: sim_recon_c06.py): this process pins itself (all
% its threads) to the CPUs of <work>/cpus_<t>.txt and sets maxNumCompThreads; then
% simReconFrame of <work>/data/Cell_5phase_small.tif with the PSF file (the OTF made in the
% call), the 2020 workflow's values (the defaults), DS, on the CPU, once untimed and then timed,
% the result folder removed after each (untimed).

cpus = strtrim(fileread(fullfile(work, sprintf('cpus_%d.txt', t))));
[s, msg] = system(sprintf('taskset -a -p -c %s %d', strrep(cpus, ' ', ','), feature('getpid')));
if s ~= 0
    error('taskset failed: %s', msg);
end
maxNumCompThreads(t);
pk5d_bench_mark('pinned');

data = [work, '/data/'];
for timed = [false, true]
    t0 = tic;
    simReconFrame([data, 'Cell_5phase_small.tif'], [data, 'PSF_5phase.tif'], 'useGPU', false, ...
        'DS', true);
    dt = toc(t0);
    rmdir([data, 'sim_recon'], 's');
    if timed
        pk5d_bench_record('frame', t, dt);
    end
end
end
