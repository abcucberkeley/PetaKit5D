function psf_analysis_b06_ml(folder, what, threads, repeats)
% PSF analysis benchmark (B06), the MATLAB side, with maxNumCompThreads set to each thread
% count: 'analysis': XR_psf_analysis_wrapper of folder/PSF with the demo's parameters
% (demo_useful_tools' "PSF analysis"), the figures not visible, parseCluster false, the DS
% folder removed before each run (untimed; needs the JVM for the figures); 'detection':
% XR_psf_detection_and_cropping of folder/beads/beads.tif (crops of 64 x 64 x 31, isolation
% distance 64 x 64 x 40), its result folder removed before each run (untimed).

nt0 = maxNumCompThreads;
psf = [folder, '/PSF/'];
beads = [folder, '/beads/'];
for t = threads
    maxNumCompThreads(t);
    times = zeros(1, repeats);
    for r = 1 : repeats
        if strcmp(what, 'analysis')
            if exist([psf, 'DS'], 'dir')
                rmdir([psf, 'DS'], 's');
            end
            t0 = tic;
            XR_psf_analysis_wrapper({psf}, 'xyPixelSize', 0.108, 'dz', 0.3, 'skewAngle', 32.45, ...
                'ChannelPatterns', {'488'}, 'Channels', 488, 'Deskew', true, 'ObjectiveScan', false, ...
                'sourceStr', 'Demo PSF', 'RWFn', {[psf, 'RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif']}, ...
                'visible', false, 'parseCluster', false, 'masterCompute', true, 'cpusPerTask', 24);
        else
            if exist([beads, 'Cropped'], 'dir')
                rmdir([beads, 'Cropped'], 's');
            end
            t0 = tic;
            XR_psf_detection_and_cropping([beads, 'beads.tif'], [beads, 'Cropped/'], ...
                'cropSize', [64, 64, 31], 'distThresh', [64, 64, 40], 'prefix', 'beads');
        end
        times(r) = toc(t0);
        close all;
    end
    if strcmp(what, 'analysis')
        pk5d_bench_record('psf_analysis', t, times);
    else
        pk5d_bench_record('psf_detection', t, times);
    end
end
maxNumCompThreads(nt0);
end
