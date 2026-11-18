function chromatic_b03_ml(work, t)
% The chromatic shift benchmark (B03), MATLAB side (Python: chromatic_b03.py): this process pins
% itself (all its threads) to the CPUs of <work>/cpus_<t>.txt and sets maxNumCompThreads; then
% XR_chromatic_shift_correction_data_wrapper of <work>/ds/ as demo_useful_tools.m calls it
% (offsets [0 0 0; 1 -5 -2], 'valid' mode, tiff results, parseCluster false), then
% XR_chromatic_shift_estimation of the demo's two PSFs (maxOffset 20), each once untimed and
% then timed, the results removed after each (untimed).

cpus = strtrim(fileread(fullfile(work, sprintf('cpus_%d.txt', t))));
[s, msg] = system(sprintf('taskset -a -p -c %s %d', strrep(cpus, ' ', ','), feature('getpid')));
if s ~= 0
    error('taskset failed: %s', msg);
end
maxNumCompThreads(t);
pk5d_bench_mark('pinned');

ds = [work, '/ds/'];
psf = {[work, '/psf/560_c.tif'], [work, '/psf/488_2_c.tif']};
res = [work, '/ml_measurement.mat'];
here = pwd;
for timed = [false, true]
    t0 = tic;
    XR_chromatic_shift_correction_data_wrapper({ds}, 'chromaticOffset', [0, 0, 0; 1, -5, -2], ...
        'resultDirName', 'Chromatic_Shift_Corrected', 'mode', 'valid', 'padValue', 0, ...
        'newOrigin', [], 'channelPatterns', {'CamA', 'CamB'}, 'psfFullpaths', {}, ...
        'maxOffset', [20, 20, 20], 'cropLength', [0, 0, 0], 'zarrFile', false, ...
        'largeFile', false, 'saveZarr', false, 'batchSize', [1024, 1024, 1024], ...
        'blockSize', [256, 256, 256], 'parseCluster', false, 'masterCompute', true, ...
        'cpusPerTask', 3);
    dt = toc(t0);
    cd(here);
    rmdir([ds, 'Chromatic_Shift_Corrected'], 's');
    if timed
        pk5d_bench_record('dataset', t, dt);
    end
end
for timed = [false, true]
    t0 = tic;
    XR_chromatic_shift_estimation(psf, 'resultFullpath', res, 'maxOffset', [20, 20, 20]);
    dt = toc(t0);
    delete(res);
    if timed
        pk5d_bench_record('estimation', t, dt);
    end
end
end
