function fft_spectrum_b04_ml(work, t)
% The FFT spectrum benchmark (B04), MATLAB side (Python: fft_spectrum_b04.py): this process
% pins itself (all its threads) to the CPUs of <work>/cpus_<t>.txt and sets maxNumCompThreads;
% then XR_fftSpectrumComputingWrapper of <work>/matlab/ds/ as demo_useful_tools.m calls it
% (zarr input, outSize 1001^3, voxels of 0.1 * 0.488 / 1.33 um, no 3D stack, parseCluster
% false), once untimed and then timed, the result folder removed after each (untimed).

cpus = strtrim(fileread(fullfile(work, sprintf('cpus_%d.txt', t))));
[s, msg] = system(sprintf('taskset -a -p -c %s %d', strrep(cpus, ' ', ','), feature('getpid')));
if s ~= 0
    error('taskset failed: %s', msg);
end
maxNumCompThreads(t);
pk5d_bench_mark('pinned');

ds = [work, '/matlab/ds/'];
here = pwd;
for timed = [false, true]
    t0 = tic;
    XR_fftSpectrumComputingWrapper({ds}, 'resultDirName', 'FFT', 'xyPixelSize', 0.108, ...
        'dz', 0.3 * sind(32.45), 'zarrFile', true, 'outSize', [1001, 1001, 1001], ...
        'outPixelSize', 0.1 * 0.488 / 1.33, 'save3DStack', false, 'parseCluster', false, ...
        'masterCompute', true, 'cpusPerTask', 24);
    dt = toc(t0);
    % (the wrapper changes the current folder)
    cd(here);
    rmdir([ds, 'FFT'], 's');
    if timed
        pk5d_bench_record('dataset', t, dt);
    end
end
end
