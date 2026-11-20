function detection_c02_ml(work, nthreads)
% The detection benchmark (C02), MATLAB side (Python: detection_c02.py): the background estimate
% and the lowSNR detector on one demo frame as XR_runDetection3D runs them, once untimed and
% then timed.
maxNumCompThreads(nthreads);
vol = double(readtiff(fullfile(work, 'frame.tif')));
sigma = [1.3758153893715306, 1.1049523723652759];
opts = {'Alpha', 0.05, 'Mask', [], 'RemoveRedundant', true, 'RedundancyRadius', 0.5, ...
    'RefineMaskLoG', false, 'WindowSize', [], 'Mode', 'xyzAc', 'FitMixtures', false, ...
    'MaxMixtures', 1, 'BackgroundCorrection', true, 'FilterByResidualSigma', true, ...
    'ResidualSigmaThresh', 3};
pk5d_bench_mark('ready');
run_once(vol, sigma, opts);
t0 = tic;
n = run_once(vol, sigma, opts);
pk5d_bench_record('detect_frame', nthreads, toc(t0));
fprintf('detections %d\n', n);
pk5d_bench_mark('done');
end

function n = run_once(vol, sigma, opts)
clear XR_pointSourceDetection3D_3
bi = XR_estimate_overall_background_information(vol, 'randseed', 1);
ps = XR_pointSourceDetection3D_3(vol, sigma, opts{:}, 'BackgroundInfo', [bi, 3]);
n = numel(ps.x);
end
