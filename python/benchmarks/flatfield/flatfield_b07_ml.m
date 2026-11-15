function flatfield_b07_ml(folder, threads, repeats)
% Flat-field benchmark (B07), the MATLAB side, with maxNumCompThreads set to each thread
% count, the results removed before each run (untimed): XR_flatfield_background_wrapper of
% folder/dark (two cameras, parseCluster false) and XR_phase_image_flat_field_correction of
% each time point of folder/phase (the demo's loop, offset 6000). With B07_KEEP set, the last
% run's phase results are copied to $B07_KEEP/matlab.

nt0 = maxNumCompThreads;
dark = [folder, '/dark'];
phase = [folder, '/phase/'];
outPath = [phase, 'ff_corrected_1/'];
d = dir([phase, '*.tif']);
fns = cellfun(@(x) [phase, x], {d.name}', 'unif', 0);
for t = threads
    maxNumCompThreads(t);
    times = zeros(1, repeats);
    for r = 1 : repeats
        res = [dark, '/matlab_flat_field_estimation'];
        if exist(res, 'dir')
            rmdir(res, 's');
        end
        t0 = tic;
        XR_flatfield_background_wrapper(dark, 'channelPatterns', {'CamA_ch0', 'CamB_ch0'}, ...
            'parseCluster', false);
        times(r) = toc(t0);
    end
    pk5d_bench_record('flatfield_background_wrapper', t, times);
    times = zeros(1, repeats);
    for r = 1 : repeats
        if exist(outPath, 'dir')
            rmdir(outPath, 's');
        end
        mkdir(outPath);
        t0 = tic;
        for i = 1 : 2
            fns_i = fns(contains(fns, sprintf('Iter_%04d', i - 1)));
            XR_phase_image_flat_field_correction(fns_i, outPath, 6000);
        end
        times(r) = toc(t0);
    end
    pk5d_bench_record('phase_flat_field', t, times);
end
maxNumCompThreads(nt0);
keep = getenv('B07_KEEP');
if ~isempty(keep)
    mkdir([keep, '/matlab']);
    copyfile([outPath, '*'], [keep, '/matlab/']);
end
end
