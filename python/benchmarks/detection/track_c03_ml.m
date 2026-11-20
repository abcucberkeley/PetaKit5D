function track_c03_ml(folder, repeats)
% Tracking benchmark (C03), the MATLAB side: XR_runTracking3D of folder/Analysis/Detection3D.mat
% with the demos' settings, the result removed before each run (untimed).

rp = [folder, '/Analysis/'];
S = load([rp, 'Detection3D.mat'], 'frameInfo');
data = struct('source', [folder, '/'], 'channels', {{[folder, '/']}}, 'movieLength', numel(S.frameInfo));
clear S
settings = loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2);
times = zeros(1, repeats);
for r = 1 : repeats
    fn = [rp, 'trackedFeatures.mat'];
    if exist(fn, 'file')
        delete(fn);
    end
    t0 = tic;
    XR_runTracking3D(data, settings, 'FileName', 'trackedFeatures.mat', 'Overwrite', true, ...
        'ResultsPath', {rp});
    times(r) = toc(t0);
end
pk5d_bench_record('track', maxNumCompThreads, times);
end
