function detection_c02c_ml(work, nthreads)
% The frame-parallel detection benchmark (C02c), MATLAB side (Python: detection_c02c.py):
% XR_runDetection3D_frame_parallel on the setup's condition ('lowSNR' without the background
% correction, XR_runDetTrack3d's detection), the results removed before each run (untimed; the
% runner keeps existing files whatever Overwrite says); once untimed and then timed.
maxNumCompThreads(nthreads);
sigma = readmatrix(fullfile(work, 'sigma.txt'), 'FileType', 'text', 'Delimiter', ' ');
ex = [fullfile(work, 'Ex1_bench_z0p5') filesep];
data = XR_loadConditionData3D(ex, {'ch1', 'ch2'}, {'gfp', 'rfp'});
rp = [data(1).source 'Analysis' filesep];
opts = {'Sigma', sigma, 'Overwrite', false, 'WindowSize', [], 'ResultsPath', {rp}, ...
    'Mode', 'xyzAc', 'FitMixtures', false, 'MaxMixtures', 2, 'DetectionMethod', 'lowSNR', ...
    'BackgroundCorrection', false};
pk5d_bench_mark('ready');
run_once(data, rp, opts);
clean(rp);
t0 = tic;
run_once(data, rp, opts);
pk5d_bench_record('detect_movie', nthreads, toc(t0));
S = load([rp 'Detection3D.mat'], 'frameInfo');
fprintf('detections %d\n', sum(arrayfun(@(f) numel(f.x), S.frameInfo)));
pk5d_bench_mark('done');
end

function run_once(data, rp, opts)
clean(rp);
clear XR_pointSourceDetection3D_3
XR_runDetection3D_frame_parallel(data, opts{:});
end

function clean(rp)
if exist(rp, 'dir')
    rmdir(rp, 's');
end
end
