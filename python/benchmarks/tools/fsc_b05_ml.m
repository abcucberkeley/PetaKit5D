function fsc_b05_ml(folder, threads, repeats)
% FSC benchmark (B05), the MATLAB side: XR_FSC_analysis_wrapper of folder/iters with the
% parameters of demo_useful_tools.m's FSC section (one-bit, dr 5, dtheta pi/6, halfSize 201,
% parseCluster false), with maxNumCompThreads set to each thread count; the results and
% figures removed before each run (untimed). The figures need MATLAB's JVM.

nt0 = maxNumCompThreads;
iters = [folder, '/iters/'];
for t = threads
    maxNumCompThreads(t);
    times = zeros(1, repeats);
    for r = 1 : repeats
        for d = {'FSCs', 'figures'}
            if exist([iters, d{1}], 'dir')
                rmdir([iters, d{1}], 's');
            end
        end
        clear XR_FSC XR_one_image_FSC_analysis
        t0 = tic;
        XR_FSC_analysis_wrapper({iters}, 'resultDirName', 'FSCs', 'xyPixelSize', 0.108, ...
            'dz', 0.3 * sind(32.45), 'dr', 5, 'dtheta', pi / 6, 'resThreshMethod', 'one-bit', ...
            'resThresh', 0.2, 'halfSize', [201, 201, 201], 'inputBbox', [], ...
            'ChannelPatterns', {'Iter'}, 'Channels', 488, 'suffix', 'decon', ...
            'iterInterval', 2, parseCluster=false, masterCompute=true, cpusPerTask=24, ...
            mccMode=false, configFile='');
        times(r) = toc(t0);
    end
    pk5d_bench_record('fsc_analysis_dataset', t, times);
end
maxNumCompThreads(nt0);
end
